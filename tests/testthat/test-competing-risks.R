test_that("competing risks reproduce both the time and the cause distribution", {
  set.seed(6)
  n   <- 20000
  sim <- sim_pexp_cr(list(a = ~ log(0.3), b = ~ log(0.7)),
    data.frame(id = seq_len(n)), cut = seq(0, 60, by = 0.5))
  ev  <- sim[sim$status == 1, ]

  # all-cause waiting time is Exp(0.3 + 0.7), cause a has probability 0.3
  expect_gt(ks.test(ev$tstop, "pexp", 1.0)$p.value, 0.01)
  expect_equal(mean(ev$to == "a"), 0.3, tolerance = 0.015)
})

test_that("the cause probabilities follow the hazards over time", {
  set.seed(7)
  n   <- 40000
  # hazard of b grows with time, so late events are more often of cause b
  sim <- sim_pexp_cr(list(a = ~ log(0.2), b = ~ log(0.02) + 0.3 * t),
    data.frame(id = seq_len(n)), cut = seq(0, 20, by = 0.05))
  ev  <- sim[sim$status == 1, ]
  share_early <- mean(ev$to[ev$tstop < 2] == "b")
  share_late  <- mean(ev$to[ev$tstop > 6] == "b")
  expect_lt(share_early, share_late)
  # at t the ratio of hazards is known: P(b | event at t) = h_b / (h_a + h_b)
  t_mid <- 8
  hb <- 0.02 * exp(0.3 * t_mid); ha <- 0.2
  near <- abs(ev$tstop - t_mid) < 0.5
  expect_equal(mean(ev$to[near] == "b"), hb / (ha + hb), tolerance = 0.05)
})

test_that("unnamed causes are numbered", {
  set.seed(8)
  sim <- sim_pexp_cr(list(~ log(0.3), ~ log(0.3)), data.frame(id = 1:200),
    cut = seq(0, 10, by = 0.1))
  expect_setequal(setdiff(levels(sim$to), "0"), c("1", "2"))
})
