# Getting started

``` r

library(simpexp)
```

## The idea

A hazard that is constant on the intervals of a grid determines
everything else in closed form. Its cumulative hazard is piecewise
linear, so the waiting time can be drawn by inverting it exactly; and
when several transitions compete, the hazards are constant within the
interval in which the event falls, so the transition that is realised
can be drawn with probabilities proportional to them. No quadrature, no
root finding, no rejection step.

The grid is therefore part of the model, not a numerical tuning
parameter. A hazard that varies quickly needs a fine grid to be
represented faithfully – and the same grid is the one on which piecewise
exponential additive models are fitted, which is what makes simulation
and estimation speak about the same quantities.

## Three interfaces, one engine

Everything is specified by a one-sided formula for the **log**-hazard.

### A single event

``` r

df  <- data.frame(x = rnorm(500), trt = rbinom(500, 1, 0.5))
cut <- seq(0, 10, by = 0.1)

sim <- sim_pexp(~ log(0.1) + 0.05 * t + 0.7 * trt - 0.3 * x, df, cut = cut)
head(sim)
#>   id episode from   to tstart      tstop status transition           x trt
#> 1  1       1    0 <NA>      0 10.0000000      0       <NA> -0.91804941   1
#> 2  2       1    0    1      0  0.1656644      1       0->1  0.73987621   1
#> 3  3       1    0    1      0  8.7447251      1       0->1  0.87236005   0
#> 4  4       1    0    1      0  6.6065062      1       0->1 -0.03432841   0
#> 5  5       1    0    1      0  1.3145779      1       0->1  1.27357501   1
#> 6  6       1    0    1      0  0.5269736      1       0->1 -1.09479437   1
```

The output is always in the same long format: one row per episode, with
`tstart`, `tstop`, `from`, `to`, `status` and the covariates of the
input data. `status = 0` marks an episode that ended without a
transition, and then `to` is `NA`. For a single event that is all there
is to it – one episode per subject.

``` r

mean(sim$status)
#> [1] 0.806
```

A hazard given as one value per interval needs no closed form at all;
the variable `interval` is the index of the current interval of `cut`:

``` r

rates <- c(0.05, 0.2, 0.35, 0.2, 0.1)
sim_pexp(~ log(rates[interval]) + 0.5 * trt, df, cut = c(0, 1, 2, 4, 7, 10)) |>
  head(3)
#>   id episode from to tstart    tstop status transition          x trt
#> 1  1       1    0  1      0 1.296862      1       0->1 -0.9180494   1
#> 2  2       1    0  1      0 2.556285      1       0->1  0.7398762   1
#> 3  3       1    0  1      0 3.170560      1       0->1  0.8723600   0
```

### Competing risks

One formula per cause; the names become the cause labels.

``` r

cr <- sim_pexp_cr(
  list(relapse = ~ log(0.1) + 0.5 * trt,
       death   = ~ log(0.03) - 0.2 * x),
  df, cut = cut)
table(cr$to, useNA = "ifany")
#> 
#>       0 relapse   death    <NA> 
#>       0     323      71     106
```

`to` is `NA` exactly for the administratively censored subjects, which
is `status == 0`.

### Multi-state processes

A transition is the triple (origin, target, log-hazard). States may
recur: the engine keeps drawing until every subject has reached a
terminal state, been censored, or run out of follow-up.

``` r

msm_sim <- sim_pexp_msm(
  list(
    transition("healthy", "ill",     ~ log(0.4) + 0.3 * x),
    transition("healthy", "dead",    ~ log(0.02)),
    transition("ill",     "healthy", ~ log(0.8)),
    transition("ill",     "dead",    ~ log(0.05) + 0.15 * t_sojourn)
  ),
  df, cut = cut)

head(msm_sim, 8)
#>   id episode    from      to   tstart     tstop status   transition          x
#> 1  1       1 healthy     ill 0.000000  1.386584      1 healthy->ill -0.9180494
#> 2  1       2     ill healthy 1.386584  1.947929      1 ill->healthy -0.9180494
#> 3  1       3 healthy     ill 1.947929  3.464626      1 healthy->ill -0.9180494
#> 4  1       4     ill healthy 3.464626  4.577484      1 ill->healthy -0.9180494
#> 5  1       5 healthy    <NA> 4.577484 10.000000      0         <NA> -0.9180494
#> 6  2       1 healthy     ill 0.000000  4.165867      1 healthy->ill  0.7398762
#> 7  2       2     ill healthy 4.165867  5.460125      1 ill->healthy  0.7398762
#> 8  2       3 healthy     ill 5.460125  5.856368      1 healthy->ill  0.7398762
#>   trt
#> 1   1
#> 2   1
#> 3   1
#> 4   1
#> 5   1
#> 6   1
#> 7   1
#> 8   1
table(msm_sim$transition, useNA = "ifany")
#> 
#>  healthy->ill healthy->dead  ill->healthy     ill->dead          <NA> 
#>          1219            67          1028            78           355
```

States that have no outgoing transition are terminal by default (`dead`
here); `terminal_states` overrides that.

``` r

# how often subjects pass through `ill`
table(table(msm_sim$id[msm_sim$to == "ill" & !is.na(msm_sim$to)]))
#> 
#>   1   2   3   4   5   6 
#> 101 138 115  68  33  10
```

## Starting elsewhere, starting later

Two optional columns of `data` control where a subject starts: `from`
gives the initial state and `t_entry` the time at which it was entered.
The first interval of such a subject starts at its entry time, so
entering in the middle of an interval is handled exactly rather than
rounded to the grid.

``` r

staggered <- data.frame(
  id      = 1:6,
  x       = 0,
  from    = c("healthy", "healthy", "ill", "ill", "ill", "healthy"),
  t_entry = c(0, 0, 1.37, 2.5, 4, 0))

sim_pexp_msm(
  list(
    transition("healthy", "ill",  ~ log(0.4)),
    transition("ill",     "dead", ~ log(0.2))),
  staggered, cut = seq(0, 10, by = 0.5)) |>
  subset(select = c(id, episode, from, to, tstart, tstop, status))
#>   id episode    from   to     tstart      tstop status
#> 1  1       1 healthy  ill 0.00000000 0.45591307      1
#> 2  1       2     ill dead 0.45591307 0.85770620      1
#> 3  2       1 healthy  ill 0.00000000 0.05610374      1
#> 4  2       2     ill dead 0.05610374 5.90193646      1
#> 5  3       1     ill dead 1.37000000 6.22890148      1
#> 6  4       1     ill dead 2.50000000 4.08344092      1
#> 7  5       1     ill dead 4.00000000 5.14905646      1
#> 8  6       1 healthy  ill 0.00000000 0.79673510      1
#> 9  6       2     ill dead 0.79673510 1.45679491      1
```

## Fitting what was simulated

The point of simulating on a grid is that the same grid carries over to
estimation. Expand the episodes into intervals, fit a Poisson model with
the log of the interval length as offset, and the coefficients are the
ones that were put in.

``` r

truth <- c(intercept = log(0.2), trt = 0.7, x = -0.3)
dat   <- data.frame(x = rnorm(3000), trt = rbinom(3000, 1, 0.5))
sim   <- sim_pexp(~ log(0.2) + 0.7 * trt - 0.3 * x, dat, cut = seq(0, 10, by = 0.5))

# one row per subject and interval at risk
expand <- function(sim, cut) {
  lo <- cut[-length(cut)]; hi <- cut[-1L]
  do.call(rbind, lapply(seq_along(lo), function(j) {
    at_risk <- sim$tstop > lo[j]
    d <- sim[at_risk, , drop = FALSE]
    d$offset <- log(pmin(d$tstop, hi[j]) - lo[j])
    d$event  <- as.integer(d$status == 1L & d$tstop <= hi[j])
    d$j      <- j
    d
  }))
}
ped <- expand(sim, seq(0, 10, by = 0.5))

fit <- glm(event ~ trt + x, offset = offset, family = poisson(), data = ped)
cbind(truth = truth, estimate = coef(fit), confint.default(fit))
#>               truth   estimate      2.5 %     97.5 %
#> intercept -1.609438 -1.6229434 -1.6779484 -1.5679383
#> trt        0.700000  0.6806815  0.6052453  0.7561178
#> x         -0.300000 -0.2903343 -0.3291241 -0.2515444
```

A constant baseline hazard was simulated, so no baseline term is needed
here. With a time-varying baseline, `j` (or a smooth function of
interval midpoints) goes into the model – that is precisely a piecewise
exponential additive model.

## Where to go next

- [`vignette("time-scales")`](https://adibender.github.io/simpexp/articles/time-scales.md)
  – which time scale a hazard may depend on, and how that decides
  whether the process is Markov, semi-Markov or non-Markov.
- [`vignette("censoring")`](https://adibender.github.io/simpexp/articles/censoring.md)
  – administrative, random and covariate-dependent censoring.
- [`vignette("validation")`](https://adibender.github.io/simpexp/articles/validation.md)
  – the same data-generating process simulated with **simpexp** and with
  other packages, and the same model fitted to both.
