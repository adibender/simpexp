#' Simulate event histories from piecewise exponential hazards
#'
#' `sim_pexp()` simulates single-event data, `sim_pexp_cr()` competing risks and
#' `sim_pexp_msm()` general multi-state processes, including recurrent states.
#' All three share one engine and one output format; they differ only in how the
#' process is specified.
#'
#' Hazards are piecewise constant on the grid `cut`, which makes the simulation
#' exact: waiting times are drawn by inverting the cumulative hazard in closed
#' form, and when several transitions compete, the realised one is drawn with
#' probabilities proportional to the transition-specific hazards within the
#' interval in which the event falls. The grid is therefore part of the model,
#' not a numerical tuning parameter: a hazard that varies quickly needs a fine
#' grid to be represented faithfully.
#'
#' @param formula One-sided formula for the log-hazard (single event).
#' @param formulas A list of one-sided formulas, one per competing cause. Names
#'   are used as the cause labels, otherwise causes are numbered.
#' @param transitions A list of [transition()] objects, or of
#'   `list(from, to, formula)`.
#' @param data A data frame of subjects, one row each. A column `id` is added if
#'   absent. Optional columns `from` (initial state) and `t_entry` (time of
#'   entry into that state) allow starting subjects in different states or at
#'   different times.
#' @param cut Time grid on which the hazards are piecewise constant. Must start
#'   at 0; its last value is the end of follow-up and therefore the time of
#'   administrative censoring.
#' @param origin Label of the initial state (single event and competing risks).
#' @param terminal_states States that end a trajectory. Defaults to all states
#'   that have no outgoing transition.
#' @param censoring_states States that end a trajectory with `status = 0`, i.e.
#'   that encode censoring as a competing event. This is how censoring that
#'   depends on covariates or on the state of the process is expressed; see
#'   [add_censoring()] for censoring applied after the fact.
#' @param max_transitions Guard against processes that never reach a terminal
#'   state; a warning is issued if it is hit.
#' @param keep_history Append the history features that the hazards were
#'   evaluated with (see [time_scales]) to the returned episodes. They are the
#'   ground truth of the data-generating process and are what a model would have
#'   to recover, so they are useful as oracle features in a simulation study.
#'
#' @return A data frame with one row per episode and the columns `id`,
#'   `episode`, `from`, `to`, `tstart`, `tstop`, `status`, `transition`, plus
#'   the columns of `data`. `status = 0` marks an episode that ended without a
#'   transition, i.e. censoring, and then `to` is `NA`.
#'
#' @seealso [time_scales] for what a hazard formula may depend on,
#'   [add_censoring()] for censoring that is not administrative.
#'
#' @examples
#' df <- data.frame(x = rnorm(50))
#'
#' # single event, hazard increasing in time
#' sim_pexp(~ log(0.1) + 0.05 * t + 0.3 * x, df, cut = seq(0, 10, by = 0.1)) |> head()
#'
#' # competing risks
#' sim_pexp_cr(list(relapse = ~ log(0.1) + 0.3 * x, death = ~ log(0.05)),
#'   df, cut = seq(0, 10, by = 0.1)) |> head()
#'
#' # illness-death with recovery; the hazard of dying while ill grows with the
#' # time already spent being ill, which makes the process non-Markov
#' sim_pexp_msm(
#'   list(
#'     transition("healthy", "ill",     ~ log(0.4) + 0.3 * x),
#'     transition("healthy", "dead",    ~ log(0.02)),
#'     transition("ill",     "healthy", ~ log(0.8) - 0.25 * t_sojourn),
#'     transition("ill",     "dead",    ~ log(0.05) + 0.15 * t_sojourn)
#'   ),
#'   df, cut = seq(0, 10, by = 0.1)) |> head()
#' @export
sim_pexp <- function(formula, data, cut, origin = "0", keep_history = FALSE) {
  assert_formula(formula)
  sim_event_history(as_transitions(formula, origin), data, cut, origin = origin,
    keep_history = keep_history)
}

#' @rdname sim_pexp
#' @export
sim_pexp_cr <- function(formulas, data, cut, origin = "0", keep_history = FALSE) {
  if (inherits(formulas, "formula")) formulas <- list(formulas)
  assert_list(formulas, min.len = 1L, types = "formula")
  sim_event_history(as_transitions(formulas, origin), data, cut, origin = origin,
    keep_history = keep_history)
}

#' @rdname sim_pexp
#' @export
sim_pexp_msm <- function(transitions, data, cut, terminal_states = NULL,
                         censoring_states = character(), max_transitions = 100L,
                         keep_history = FALSE) {
  sim_event_history(as_transitions(transitions), data, cut,
    terminal_states = terminal_states, censoring_states = censoring_states,
    max_transitions = max_transitions, keep_history = keep_history)
}

# ---- engine ----------------------------------------------------------------

sim_event_history <- function(trans, data, cut, origin = NULL,
                              terminal_states = NULL,
                              censoring_states = character(),
                              max_transitions = 100L,
                              keep_history = FALSE) {
  assert_data_frame(data, min.rows = 1L)
  cut   <- check_cut(cut)
  horizon <- cut[length(cut)]
  str   <- transition_structure(trans)
  assert_character(censoring_states, any.missing = FALSE)
  if (is.null(terminal_states)) terminal_states <- str$absorbing
  terminal_states <- as.character(terminal_states)
  assert_int(max_transitions, lower = 1L)
  assert_flag(keep_history)

  data <- as.data.frame(data)
  if (!("id" %in% names(data))) data$id <- seq_len(nrow(data))
  if (anyDuplicated(data$id)) stop("`id` must be unique.", call. = FALSE)
  covariates <- setdiff(names(data), c("id", "from", "t_entry"))

  state <- if ("from" %in% names(data)) as.character(data$from) else {
    if (is.null(origin)) origin <- str$transient[1L]
    rep.int(as.character(origin), nrow(data))
  }
  unknown <- setdiff(state, str$states)
  if (length(unknown)) {
    stop("No transitions defined for initial state(s): ",
      paste(unknown, collapse = ", "), ".", call. = FALSE)
  }
  entry <- if ("t_entry" %in% names(data)) as.numeric(data$t_entry) else rep.int(0, nrow(data))

  # history, carried along as one row per subject
  hist <- init_history(data$id, state, str$states)

  at_risk <- data.frame(row = seq_len(nrow(data)), state = state, entry = entry,
    episode = 1L, stringsAsFactors = FALSE)
  at_risk <- at_risk[!(at_risk$state %in% c(terminal_states, censoring_states)) &
      at_risk$entry < horizon, , drop = FALSE]

  episodes <- list()
  round <- 0L
  while (nrow(at_risk) > 0L) {
    round <- round + 1L
    if (round > max_transitions) {
      warning(nrow(at_risk), " subject(s) reached `max_transitions` without ",
        "reaching a terminal state; their trajectories are truncated.", call. = FALSE)
      break
    }
    next_risk <- list()
    for (s in unique(at_risk$state)) {
      idx <- at_risk[at_risk$state == s, , drop = FALSE]
      out <- sim_one_state(s, idx, data, hist, trans, str, cut, horizon, covariates,
        keep_history)
      episodes[[length(episodes) + 1L]] <- out$episodes   # append, never overwrite:
      hist <- out$history                                 # states may recur
      cont <- out$episodes$status == 1L &
        !(out$episodes$to %in% c(terminal_states, censoring_states)) &
        out$episodes$tstop < horizon
      if (any(cont)) {
        next_risk[[s]] <- data.frame(
          row     = out$rows[cont],
          state   = out$episodes$to[cont],
          entry   = out$episodes$tstop[cont],
          episode = out$episodes$episode[cont] + 1L,
          stringsAsFactors = FALSE)
      }
    }
    at_risk <- if (length(next_risk)) {
      do.call(rbind, c(next_risk, list(make.row.names = FALSE)))
    } else {
      at_risk[0L, , drop = FALSE]
    }
  }

  finalise(episodes, data, covariates, str)
}

# One round for the subjects that are currently in state `s`.
sim_one_state <- function(s, idx, data, hist, trans, str, cut, horizon, covariates,
                          keep_history = FALSE) {
  k        <- str$by_origin[[s]]
  targets  <- str$to[k]
  rows     <- idx$row

  base <- take_rows(data, rows, union("id", covariates))
  base$t_entry <- idx$entry
  features <- history_features(hist, rows)
  base <- cbind(base, features)
  exp_data <- expand_intervals(base, cut, entry = idx$entry)
  exp_data$t_sojourn <- exp_data$tend - exp_data$t_entry

  hazards <- vapply(k, function(i) {
    eta <- eval(trans[[i]]$formula[[2L]], exp_data, environment(trans[[i]]$formula))
    if (length(eta) == 1L) eta <- rep.int(eta, nrow(exp_data))
    assert_numeric(eta, any.missing = FALSE, finite = TRUE,
      .var.name = paste0("log-hazard of ", format(trans[[i]])))
    exp(eta)
  }, numeric(nrow(exp_data)))
  dim(hazards) <- c(nrow(exp_data), length(k))

  rate  <- rowSums(hazards)
  drawn <- draw_waiting_time(rate, exp_data$tstart, exp_data$tend, exp_data$.row)

  to     <- rep.int(NA_character_, length(rows))
  status <- as.integer(drawn$has_event)
  tstop  <- ifelse(drawn$has_event, drawn$time, horizon)
  if (any(drawn$has_event)) {
    hit <- drawn$interval[drawn$has_event]
    to[drawn$has_event] <- targets[draw_target(hazards, hit)]
  }

  ep <- data.frame(
    .row    = rows,
    id      = data$id[rows],
    episode = idx$episode,
    from    = s,
    to      = to,
    tstart  = idx$entry,
    tstop   = tstop,
    status  = status,
    stringsAsFactors = FALSE)
  if (keep_history) ep <- cbind(ep, t_entry = idx$entry, features)

  list(episodes = ep, rows = rows,
    history = update_history(hist, rows, s, to, status, tstop - idx$entry))
}

finalise <- function(episodes, data, covariates, str) {
  out <- do.call(rbind, c(episodes, list(make.row.names = FALSE)))
  out <- out[order(out$id, out$tstart), , drop = FALSE]
  out$transition <- factor(ifelse(is.na(out$to), NA_character_,
    paste0(out$from, "->", out$to)),
    levels = unique(paste0(str$from, "->", str$to)))
  out$from <- factor(out$from, levels = str$states)
  out$to   <- factor(out$to, levels = str$states)
  front <- c("id", "episode", "from", "to", "tstart", "tstop", "status", "transition")
  extra <- setdiff(names(out), c(front, ".row"))
  out <- cbind(out[c(front, extra)], take_rows(data, out$.row, covariates))
  rownames(out) <- NULL
  structure(out, class = c("sim_df", "data.frame"))
}
