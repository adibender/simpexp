# The point of simulating from piecewise exponential hazards is that the same
# model can be fitted back: expand the simulated episodes into intervals, count
# events per interval and fit a Poisson model with log exposure as offset. The
# coefficients must come back.
test_that("a piecewise exponential model recovers the simulated hazard", {
  skip_if_not_installed("stats")
  set.seed(31)
  n   <- 4000
  cut <- seq(0, 10, by = 0.5)
  df  <- data.frame(x = rbinom(n, 1, 0.5), z = rnorm(n))
  sim <- sim_pexp(~ log(0.2) + 0.7 * x - 0.3 * z, df, cut = cut)

  # interval expansion of the observed follow-up
  ped <- simpexp:::expand_intervals(sim, cut)
  ped <- ped[ped$tstart < sim$tstop[ped$.row], , drop = FALSE]
  ped$tend_obs <- pmin(ped$tend, sim$tstop[ped$.row])
  ped$dt       <- ped$tend_obs - ped$tstart
  ped$event    <- as.integer(sim$status[ped$.row] == 1L &
    sim$tstop[ped$.row] <= ped$tend & sim$tstop[ped$.row] > ped$tstart)

  fit <- glm(event ~ x + z + offset(log(dt)), data = ped, family = poisson())
  co  <- coef(fit)

  expect_equal(unname(co[["(Intercept)"]]), log(0.2), tolerance = 0.1)
  expect_equal(unname(co[["x"]]),  0.7, tolerance = 0.1)
  expect_equal(unname(co[["z"]]), -0.3, tolerance = 0.1)
})

test_that("a time-varying hazard is recovered interval by interval", {
  set.seed(32)
  n   <- 8000
  cut <- seq(0, 6, by = 1)
  sim <- sim_pexp(~ log(0.1) + 0.4 * t, data.frame(id = seq_len(n)), cut = cut)

  ped <- simpexp:::expand_intervals(sim, cut)
  ped <- ped[ped$tstart < sim$tstop[ped$.row], , drop = FALSE]
  ped$dt    <- pmin(ped$tend, sim$tstop[ped$.row]) - ped$tstart
  ped$event <- as.integer(sim$status[ped$.row] == 1L &
    sim$tstop[ped$.row] <= ped$tend & sim$tstop[ped$.row] > ped$tstart)

  fit <- glm(event ~ factor(tend) + offset(log(dt)), data = ped, family = poisson())
  # the hazard of interval (tend - 1, tend] was specified as 0.1 * exp(0.4 * tend)
  fitted_hazard <- exp(c(coef(fit)[1L], coef(fit)[1L] + coef(fit)[-1L]))
  expect_equal(unname(fitted_hazard), 0.1 * exp(0.4 * cut[-1L]), tolerance = 0.12)
})
