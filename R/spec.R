#' Specify one transition of an event-history process
#'
#' A transition is the triple (origin state, target state, log-hazard). The
#' log-hazard is given as a one-sided formula whose right-hand side is evaluated
#' on the interval-expanded data, so it may use covariates of `data` as well as
#' the time scales and history features that the simulator provides (see
#' [time_scales]).
#'
#' @param from Origin state, a single string or number.
#' @param to Target state, a single string or number.
#' @param formula One-sided formula for the log-hazard of this transition, e.g.
#'   `~ log(0.1) + 0.3 * x + 0.05 * t_sojourn`.
#'
#' @return An object of class `pexp_transition`.
#' @seealso [sim_pexp_msm()], [time_scales]
#' @examples
#' transition(0, 1, ~ log(0.5) + 0.2 * x)
#' @export
transition <- function(from, to, formula) {
  assert_scalar(from)
  assert_scalar(to)
  assert_formula(formula)
  if (length(formula) != 2L) {
    stop("`formula` must be one-sided, e.g. ~ log(0.1) + 0.3 * x.", call. = FALSE)
  }
  from <- as.character(from)
  to   <- as.character(to)
  if (identical(from, to)) {
    stop("A transition must change the state: `from` and `to` are both '", from, "'.",
      call. = FALSE)
  }
  structure(list(from = from, to = to, formula = formula), class = "pexp_transition")
}

#' @export
format.pexp_transition <- function(x, ...) {
  paste0(x$from, " -> ", x$to, ": ", deparse1(x$formula[[2L]]))
}

#' @export
print.pexp_transition <- function(x, ...) {
  cat("<pexp_transition> ", format(x), "\n", sep = "")
  invisible(x)
}

# Bring the different ways of specifying a process into one internal form: a
# list of pexp_transition objects. Accepts a single formula (single event), a
# list of formulas (competing risks, names become target states) and a list of
# transitions (multi state).
as_transitions <- function(x, origin = "0") {
  if (inherits(x, "formula")) {
    return(list(transition(origin, "1", x)))
  }
  if (inherits(x, "pexp_transition")) {
    return(list(x))
  }
  assert_list(x, min.len = 1L, .var.name = "transitions")
  if (all(vapply(x, inherits, logical(1), "pexp_transition"))) {
    return(x)
  }
  if (all(vapply(x, inherits, logical(1), "formula"))) {
    targets <- names(x)
    if (is.null(targets) || any(!nzchar(targets))) {
      targets <- as.character(seq_along(x))
    }
    return(Map(transition, from = origin, to = targets, formula = x))
  }
  # list(list(from =, to =, formula =), ...), as used in some simulation scripts
  lapply(x, function(tr) {
    if (inherits(tr, "pexp_transition")) return(tr)
    assert_list(tr, min.len = 3L)
    from <- if (!is.null(tr$from)) tr$from else tr[[1L]]
    to   <- if (!is.null(tr$to)) tr$to else tr[[2L]]
    form <- if (!is.null(tr$formula)) tr$formula else tr[[3L]]
    transition(from, to, form)
  })
}

# States that can be left, states that cannot, and the transitions by origin.
transition_structure <- function(trans) {
  from <- vapply(trans, `[[`, character(1), "from")
  to   <- vapply(trans, `[[`, character(1), "to")
  dup  <- duplicated(paste0(from, "->", to))
  if (any(dup)) {
    stop("Duplicated transition(s): ",
      paste(unique(paste0(from, "->", to)[dup]), collapse = ", "), ".", call. = FALSE)
  }
  states <- unique(c(from, to))
  list(
    from = from,
    to = to,
    states = states,
    transient = unique(from),
    absorbing = setdiff(states, unique(from)),
    by_origin = split(seq_along(trans), factor(from, levels = unique(from)))
  )
}
