# Per-subject history of the process. Kept as matrices with one row per subject
# of `data`, so that lookups and updates are plain indexing: the engine walks
# through rounds, and in each round only the subjects currently in one state are
# touched.
init_history <- function(id, state, states) {
  n <- length(id)
  n_visits <- matrix(0L, nrow = n, ncol = length(states),
    dimnames = list(NULL, states))
  t_total  <- matrix(0, nrow = n, ncol = length(states),
    dimnames = list(NULL, states))
  n_visits[cbind(seq_len(n), match(state, states))] <- 1L
  list(n_visits = n_visits, t_total = t_total, n_transitions = integer(n),
    states = states)
}

# Feature columns handed to the hazard formulas. Names are made syntactic so
# that states such as "0" or "in hospital" can still be referred to in a
# formula; the mapping is state -> n_visits_<state> / t_total_<state>.
history_features <- function(hist, rows) {
  nm <- make.names(hist$states)
  out <- data.frame(n_transitions = hist$n_transitions[rows])
  for (j in seq_along(hist$states)) {
    out[[paste0("n_visits_", nm[j])]] <- hist$n_visits[rows, j]
    out[[paste0("t_total_", nm[j])]]  <- hist$t_total[rows, j]
  }
  out
}

# Record one round: add the time just spent in `from` to the time budget of that
# state and, for the subjects that actually moved, count the visit to `to` and
# the transition. Censored subjects contribute time but no transition.
update_history <- function(hist, rows, from, to, status, duration) {
  j_from <- match(from, hist$states)
  hist$t_total[cbind(rows, j_from)] <- hist$t_total[cbind(rows, j_from)] + duration
  moved <- status == 1L & !is.na(to)
  if (any(moved)) {
    j_to <- match(to[moved], hist$states)
    hist$n_visits[cbind(rows[moved], j_to)] <-
      hist$n_visits[cbind(rows[moved], j_to)] + 1L
    hist$n_transitions[rows[moved]] <- hist$n_transitions[rows[moved]] + 1L
  }
  hist
}
