#' Censor simulated event histories
#'
#' The simulation itself censors administratively at the end of the time grid.
#' `add_censoring()` adds censoring that happens earlier: it cuts every
#' trajectory at its censoring time, shortens the episode that contains it and
#' drops everything after it.
#'
#' Censoring that depends on covariates or on the state of the process does not
#' need this function: give the process a transition into a censoring state and
#' name that state in `censoring_states` of [sim_pexp_msm()]. The engine then
#' draws the censoring time from its own hazard, in competition with the other
#' transitions, which is the natural way to express dependent censoring --
#' \emph{dependent} in the sense of covariates, not of the unobserved future.
#'
#' @param data Simulated event histories, as returned by [sim_pexp()],
#'   [sim_pexp_cr()] or [sim_pexp_msm()].
#' @param times Either a numeric vector of censoring times, named by `id` or in
#'   the order of the subjects, or a function of the number of subjects
#'   returning such a vector. Defaults to exponential censoring times, see
#'   `rate`.
#' @param rate Rate of the exponential censoring distribution, used when
#'   `times` is not given.
#' @param keep_censoring_time Append the censoring time of each subject as the
#'   column `censoring_time`.
#'
#' @return `data` with trajectories cut at the censoring times. Episodes that
#'   contain the censoring time get `status = 0`, `to = NA` and `tstop` equal to
#'   the censoring time; later episodes are removed. Subjects censored before
#'   their first episode begins drop out entirely.
#'
#' @examples
#' df  <- data.frame(x = rnorm(100))
#' sim <- sim_pexp(~ log(0.2) + 0.3 * x, df, cut = seq(0, 10, by = 0.1))
#' add_censoring(sim, rate = 0.1) |> head()
#'
#' # fixed end of follow-up per subject, e.g. staggered study entry
#' add_censoring(sim, times = runif(100, 2, 10)) |> head()
#' @export
add_censoring <- function(data, times = NULL, rate = 0.1,
                          keep_censoring_time = FALSE) {
  assert_data_frame(data, min.rows = 1L)
  ids <- unique(data$id)
  if (is.function(times)) times <- times(length(ids))
  if (is.null(times)) {
    assert_number(rate, lower = 0)
    times <- rexp(length(ids), rate)
  }
  assert_numeric(times, len = length(ids), any.missing = FALSE, lower = 0)
  if (!is.null(names(times))) {
    times <- times[match(as.character(ids), names(times))]
    if (anyNA(times)) stop("`times` must be named by every `id` in `data`.", call. = FALSE)
  }
  cens <- times[match(data$id, ids)]

  keep <- data$tstart < cens
  out  <- data[keep, , drop = FALSE]
  cens <- cens[keep]
  censored <- out$tstop > cens
  out$tstop[censored]  <- cens[censored]
  out$status[censored] <- 0L
  out$to[censored]     <- NA
  out$transition[censored] <- NA
  if (keep_censoring_time) out$censoring_time <- cens
  rownames(out) <- NULL
  out
}
