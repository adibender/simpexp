test_that("single-event simulation matches the specified distribution", {
  set.seed(3)
  n   <- 20000
  sim <- sim_pexp(~ log(0.3), data.frame(id = seq_len(n)), cut = seq(0, 60, by = 0.5))

  expect_s3_class(sim, "sim_df")
  expect_equal(nrow(sim), n)
  expect_setequal(names(sim),
    c("id", "episode", "from", "to", "tstart", "tstop", "status", "transition"))
  expect_true(all(sim$tstart == 0))
  expect_equal(mean(sim$tstop[sim$status == 1]), 1 / 0.3, tolerance = 0.02)
  expect_gt(ks.test(sim$tstop[sim$status == 1], "pexp", 0.3)$p.value, 0.01)
})

test_that("covariate effects come out as specified", {
  set.seed(4)
  n   <- 40000
  df  <- data.frame(x = rbinom(n, 1, 0.5))
  sim <- sim_pexp(~ log(0.2) + 0.7 * x, df, cut = seq(0, 80, by = 0.5))
  ev  <- sim[sim$status == 1, ]
  # mean waiting time ratio equals exp(-beta) for constant hazards
  ratio <- mean(ev$tstop[ev$x == 1]) / mean(ev$tstop[ev$x == 0])
  expect_equal(ratio, exp(-0.7), tolerance = 0.03)
})

test_that("administrative censoring happens at the end of the grid", {
  set.seed(5)
  sim <- sim_pexp(~ log(0.05), data.frame(id = 1:500), cut = seq(0, 5, by = 0.1))
  cens <- sim[sim$status == 0, ]
  expect_true(all(cens$tstop == 5))
  expect_true(all(is.na(cens$to)))
  expect_gt(nrow(cens), 0)
})
