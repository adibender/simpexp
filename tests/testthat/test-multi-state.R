illness_death <- function(recovery = TRUE, sojourn = 0) {
  tr <- list(
    transition("healthy", "ill",  ~ log(0.4)),
    transition("healthy", "dead", ~ log(0.02)),
    transition("ill", "dead",
      as.formula(sprintf("~ log(0.05) + %f * t_sojourn", sojourn)))
  )
  if (recovery) tr <- c(tr, list(transition("ill", "healthy", ~ log(0.8))))
  tr
}

test_that("recurrent states are kept, not overwritten", {
  set.seed(11)
  n   <- 500
  sim <- sim_pexp_msm(illness_death(recovery = TRUE), data.frame(id = seq_len(n)),
    cut = seq(0, 10, by = 0.1))

  expect_equal(length(unique(sim$id)), n)
  # with recovery, subjects pass through "healthy" more than once
  visits <- tapply(sim$from == "healthy", sim$id, sum)
  expect_gt(max(visits), 2)
  # episodes of a subject are contiguous and ordered
  by_id <- split(sim[c("episode", "tstart", "tstop", "from", "to")], sim$id)
  gaps <- vapply(by_id, function(d) {
    d <- d[order(d$episode), ]
    max(abs(d$tstart[-1] - d$tstop[-nrow(d)]), 0)
  }, numeric(1))
  expect_true(all(gaps < 1e-9))
  chained <- vapply(by_id, function(d) {
    d <- d[order(d$episode), ]
    nrow(d) == 1L || all(as.character(d$from[-1]) == as.character(d$to[-nrow(d)]))
  }, logical(1))
  expect_true(all(chained))
})

test_that("a progressive and a recurrent process both keep every subject", {
  set.seed(12)
  n <- 500
  prog <- sim_pexp_msm(illness_death(recovery = FALSE), data.frame(id = seq_len(n)),
    cut = seq(0, 10, by = 0.1))
  recu <- sim_pexp_msm(illness_death(recovery = TRUE), data.frame(id = seq_len(n)),
    cut = seq(0, 10, by = 0.1))
  expect_equal(length(unique(prog$id)), n)
  expect_equal(length(unique(recu$id)), n)
  expect_gt(nrow(recu), nrow(prog))
})

test_that("terminal states end the trajectory", {
  set.seed(13)
  sim <- sim_pexp_msm(illness_death(), data.frame(id = 1:300),
    cut = seq(0, 10, by = 0.1))
  last <- do.call(rbind, lapply(split(sim, sim$id), function(d) d[which.max(d$episode), ]))
  expect_true(all(is.na(last$to) | last$to == "dead" | last$tstop == 10))
  expect_false(any(sim$from == "dead"))
})

test_that("the sojourn time enters the hazard as specified", {
  set.seed(14)
  n <- 4000
  grid <- seq(0, 20, by = 0.05)
  flat  <- sim_pexp_msm(illness_death(sojourn = 0),   data.frame(id = seq_len(n)), cut = grid)
  steep <- sim_pexp_msm(illness_death(sojourn = 0.4), data.frame(id = seq_len(n)), cut = grid)
  # with a constant recovery hazard, a death hazard that grows with the sojourn
  # lets a larger share of illness episodes end in death
  share_dead <- function(x) {
    ill <- x[which(x$from == "ill"), ]
    mean(!is.na(ill$to) & ill$to == "dead")
  }
  expect_gt(share_dead(steep), share_dead(flat))
})

test_that("history features are available and are counted correctly", {
  set.seed(15)
  n <- 4000
  sim <- sim_pexp_msm(
    list(
      transition("healthy", "ill",     ~ log(0.5)),
      transition("ill",     "healthy", ~ log(0.9)),
      # the hazard of dying rises with the number of previous illness episodes
      transition("ill",     "dead",    ~ log(0.02) + 0.5 * n_visits_ill)
    ),
    data.frame(id = seq_len(n)), cut = seq(0, 10, by = 0.1), keep_history = TRUE)

  # the empirical hazard of dying while ill has to follow the specification:
  # h = 0.02 * exp(0.5 * n_visits_ill), estimated as events per person-time
  ill <- sim[which(sim$from == "ill"), ]
  ill$dead <- !is.na(ill$to) & ill$to == "dead"
  ill$time <- ill$tstop - ill$tstart
  rate <- tapply(ill$dead, ill$n_visits_ill, sum) /
    tapply(ill$time, ill$n_visits_ill, sum)
  k <- as.integer(names(rate))
  expect_equal(as.vector(rate[k <= 3]), 0.02 * exp(0.5 * k[k <= 3]), tolerance = 0.25)
  expect_true(all(diff(rate[k <= 3]) > 0))
})

test_that("subjects can start in different states and at different times", {
  set.seed(16)
  df <- data.frame(id = 1:200, from = rep(c("healthy", "ill"), each = 100),
    t_entry = rep(c(0, 2), each = 100))
  sim <- sim_pexp_msm(illness_death(), df, cut = seq(0, 10, by = 0.1))
  first <- do.call(rbind, lapply(split(sim, sim$id), function(d) d[which.min(d$episode), ]))
  expect_equal(as.character(first$from[first$id <= 100]), rep("healthy", 100))
  expect_true(all(first$tstart[first$id > 100] == 2))
})

test_that("malformed specifications are refused", {
  expect_error(transition("a", "a", ~ log(1)), "must change the state")
  expect_error(transition("a", "b", y ~ x), "one-sided")
  expect_error(
    sim_pexp_msm(list(transition("a", "b", ~ log(1)), transition("a", "b", ~ log(2))),
      data.frame(id = 1), cut = c(0, 1)),
    "Duplicated transition")
  expect_error(sim_pexp(~ log(0.1), data.frame(id = 1), cut = c(1, 2)), "start at 0")
  expect_error(sim_pexp(~ log(0.1), data.frame(id = 1), cut = c(0, 2, 1)), "increasing")
  expect_error(sim_pexp(~ log(0.1), data.frame(id = c(1, 1)), cut = c(0, 1)), "unique")
})
