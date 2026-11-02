# Draw waiting times from a hazard that is piecewise constant on the intervals
# of the expanded data.
#
# The cumulative hazard of a piecewise constant hazard is piecewise linear, so
# its inverse is available in closed form: draw r ~ Exp(1), find the interval in
# which the cumulative hazard first reaches r, and solve linearly within it.
# This is exact -- no quadrature and no interpolation of an inverse -- and it is
# the reason the hazards have to be specified on a grid.
#
# `rate` is the all-cause rate of the interval, `group` identifies the subject
# and must be sorted. Returns one row per group with the event time and whether
# an event happened before the end of the grid.
draw_waiting_time <- function(rate, tstart, tend, group) {
  if (any(rate < 0)) {
    stop("Hazards must be non-negative; check the transition formulas.", call. = FALSE)
  }
  grp   <- match(group, unique(group))
  n_grp <- max(grp)
  dH    <- rate * (tend - tstart)
  H     <- cumsum_by(dH, grp)
  r     <- rexp(n_grp)

  hit   <- H >= r[grp]
  # index of the first interval per subject in which the cumulative hazard
  # reaches the drawn quantile; NA if it never does
  first <- rep.int(NA_integer_, n_grp)
  if (any(hit)) {
    idx   <- which(hit)
    keep  <- idx[!duplicated(grp[idx])]
    first[grp[keep]] <- keep
  }
  has_event <- !is.na(first)

  time <- rep.int(NA_real_, n_grp)
  k    <- first[has_event]
  time[has_event] <- tstart[k] + (r[has_event] - (H[k] - dH[k])) / rate[k]

  list(
    time      = time,
    has_event = has_event,
    interval  = first
  )
}

# Draw which of the competing transitions is realised, given that the event
# happens in interval `interval`. Within an interval all hazards are constant,
# so the probabilities are proportional to the transition-specific hazards --
# no approximation is involved.
#
# `hazards` is a matrix with one column per competing transition, one row per
# interval of the expanded data.
draw_target <- function(hazards, interval) {
  h <- hazards[interval, , drop = FALSE]
  u <- runif(nrow(h)) * rowSums(h)
  # first column whose cumulative hazard share exceeds the uniform draw
  max.col(t(apply(h, 1L, cumsum)) >= u, ties.method = "first")
}
