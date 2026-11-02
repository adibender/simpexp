test_that("waiting times from a constant hazard are exponential", {
  set.seed(1)
  n   <- 20000
  dat <- data.frame(id = seq_len(n))
  e   <- simpexp:::expand_intervals(dat, cut = seq(0, 60, by = 0.5))
  res <- simpexp:::draw_waiting_time(rep(0.5, nrow(e)), e$tstart, e$tend, e$id)

  expect_true(all(res$has_event))
  expect_equal(mean(res$time), 2, tolerance = 0.02)
  expect_gt(ks.test(res$time, "pexp", 0.5)$p.value, 0.01)
})

test_that("a hazard specified on a fine grid reproduces its distribution", {
  set.seed(2)
  n   <- 20000
  dat <- data.frame(id = seq_len(n))
  e   <- simpexp:::expand_intervals(dat, cut = seq(0, 40, by = 0.01))
  h   <- (1.3 / 8) * (e$t / 8)^0.3          # Weibull(shape = 1.3, scale = 8)
  res <- simpexp:::draw_waiting_time(h, e$tstart, e$tend, e$id)

  expect_equal(mean(res$time[res$has_event]), 8 * gamma(1 + 1 / 1.3), tolerance = 0.01)
  expect_gt(ks.test(res$time[res$has_event], "pweibull", 1.3, 8)$p.value, 0.01)
})

test_that("subjects entering mid-interval start there, not at the interval bound", {
  e <- simpexp:::expand_intervals(data.frame(id = 1:2), cut = c(0, 1, 2, 3),
    entry = c(0, 1.5))
  expect_equal(e$tstart[e$id == 2], c(1.5, 2))
  expect_equal(e$tend[e$id == 2], c(2, 3))
})

test_that("negative hazards are refused", {
  expect_error(
    simpexp:::draw_waiting_time(c(-1, 1), c(0, 1), c(1, 2), c(1, 1)),
    "non-negative")
})
