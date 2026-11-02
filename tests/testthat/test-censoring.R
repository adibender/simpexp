test_that("censoring cuts trajectories and marks the episode", {
  set.seed(21)
  sim  <- sim_pexp(~ log(0.2), data.frame(id = 1:2000), cut = seq(0, 20, by = 0.1))
  cens <- add_censoring(sim, rate = 0.2, keep_censoring_time = TRUE)

  expect_true(all(cens$tstop <= cens$censoring_time + 1e-9))
  expect_true(all(is.na(cens$to[cens$status == 0])))
  expect_lt(mean(cens$status), mean(sim$status))
  expect_equal(length(unique(cens$id)), length(unique(sim$id)))
})

test_that("censoring with given times is exact", {
  set.seed(22)
  sim  <- sim_pexp(~ log(0.3), data.frame(id = 1:500), cut = seq(0, 20, by = 0.1))
  cens <- add_censoring(sim, times = rep(1, 500))
  expect_true(all(cens$tstop <= 1 + 1e-9))
  expect_true(all(cens$status[cens$tstop == 1] == 0))
})

test_that("censoring removes episodes that start after the censoring time", {
  set.seed(23)
  sim <- sim_pexp_msm(
    list(transition("a", "b", ~ log(1)), transition("b", "a", ~ log(1)),
         transition("b", "c", ~ log(0.1))),
    data.frame(id = 1:300), cut = seq(0, 10, by = 0.05))
  cens <- add_censoring(sim, times = rep(2, 300))
  expect_true(all(cens$tstart < 2))
  expect_lt(nrow(cens), nrow(sim))
  expect_equal(length(unique(cens$id)), 300)
})

test_that("censoring can also be specified as a competing transition", {
  set.seed(24)
  n <- 20000
  sim <- sim_pexp_msm(
    list(transition("0", "event", ~ log(0.3)),
         transition("0", "lost",  ~ log(0.1) + 0.5 * x)),
    data.frame(id = seq_len(n), x = rbinom(n, 1, 0.5)),
    cut = seq(0, 40, by = 0.1), censoring_states = "lost")
  # the censoring hazard depends on x, the event hazard does not
  expect_gt(mean(sim$to[sim$x == 1] == "lost", na.rm = TRUE),
            mean(sim$to[sim$x == 0] == "lost", na.rm = TRUE))
  # with constant hazards the share of events is h_event / (h_event + h_cens)
  for (xx in c(0, 1)) {
    ev <- sim[sim$x == xx & !is.na(sim$to), ]
    expect_equal(mean(ev$to == "event"), 0.3 / (0.3 + 0.1 * exp(0.5 * xx)),
      tolerance = 0.02)
  }
})
