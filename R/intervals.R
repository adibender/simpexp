#' Time scales and history features available in a hazard formula
#'
#' The right-hand side of a transition formula is evaluated on the
#' interval-expanded data. Besides the columns of `data`, the following
#' variables are available.
#'
#' @section Time scales:
#' \describe{
#'   \item{`tstart`, `tend`}{Boundaries of the current interval on the
#'     process time scale. The hazard is constant on `(tstart, tend]`.}
#'   \item{`t`}{Alias of `tend`: time since the origin of the process
#'     ("clock forward"). The hazard of an interval is evaluated at its end,
#'     the same convention under which piecewise exponential models are fitted,
#'     so a hazard specified here and one estimated from the simulated data
#'     refer to the same time point.}
#'   \item{`t_entry`}{Time at which the current state was entered.}
#'   \item{`t_sojourn`}{Time spent in the current state, `tend - t_entry`
#'     ("clock reset").}
#' }
#'
#' @section History features:
#' \describe{
#'   \item{`n_transitions`}{Number of transitions made before the current one.}
#'   \item{`n_visits_<state>`}{Number of visits to `<state>` so far, including
#'     the current state.}
#'   \item{`t_total_<state>`}{Time spent in `<state>` before the current
#'     interval.}
#' }
#' A process whose hazards use only `t` is Markov; one that uses `t_sojourn`
#' is semi-Markov; one that uses the history features is neither.
#'
#' @name time_scales
NULL

# Validate a time grid: numeric, strictly increasing, starting at 0.
check_cut <- function(cut) {
  assert_numeric(cut, min.len = 2L, any.missing = FALSE, finite = TRUE,
    .var.name = "cut")
  if (is.unsorted(cut, strictly = TRUE)) {
    stop("`cut` must be strictly increasing.", call. = FALSE)
  }
  if (!isTRUE(all.equal(cut[1L], 0))) {
    stop("`cut` must start at 0.", call. = FALSE)
  }
  as.numeric(cut)
}

# Subset rows of a data frame without the row-name bookkeeping of
# `[.data.frame`. When `rows` repeats, that bookkeeping (make.unique on the row
# names, whose result is discarded) dominates the runtime of a simulation.
take_rows <- function(data, rows, cols = names(data)) {
  out <- lapply(data[cols], function(col) col[rows])
  attr(out, "row.names") <- .set_row_names(length(rows))
  class(out) <- "data.frame"
  out
}

# Expand one row per subject into one row per interval of `cut` that lies after
# the subject's entry time. The first interval of a subject starts at its entry
# time, so that subjects entering mid-interval are handled exactly.
#
# Returns the expanded data with columns tstart, tend, t, and the index of the
# originating row in `.row`.
expand_intervals <- function(data, cut, entry = rep(0, nrow(data))) {
  n_int <- length(cut) - 1L
  lower <- cut[-length(cut)]
  upper <- cut[-1L]
  # first interval of each subject: the one whose upper bound exceeds entry
  first <- findInterval(entry, upper, left.open = TRUE) + 1L
  first <- pmin(first, n_int)
  times <- n_int - first + 1L
  row   <- rep.int(seq_len(nrow(data)), times)
  int   <- sequence(times, from = first)
  out   <- take_rows(data, row)
  out$.row   <- row
  out$tstart <- pmax(lower[int], entry[row])
  out$tend   <- upper[int]
  out$t      <- out$tend
  out
}

# Cumulative sums within the blocks given by `group` (which must be sorted).
cumsum_by <- function(x, group) {
  cs    <- cumsum(x)
  start <- !duplicated(group)
  # value of the cumulative sum just before each block begins
  at     <- which(start)
  offset <- c(0, cs)[at]
  cs - rep.int(offset, diff(c(at, length(x) + 1L)))
}
