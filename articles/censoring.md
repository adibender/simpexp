# Censoring

``` r

library(simpexp)
```

There are three kinds of censoring, and they are three different things,
not three settings of one option.

## Administrative censoring happens by itself

The last value of `cut` is the end of follow-up. Whoever has not had an
event by then gets `status = 0` and `to = NA`.

``` r

df  <- data.frame(x = rnorm(2000))
sim <- sim_pexp(~ log(0.05) + 0.3 * x, df, cut = seq(0, 10, by = 0.1))
table(status = sim$status)
#> status
#>    0    1 
#> 1191  809
range(sim$tstop[sim$status == 0])
#> [1] 10 10
```

## Random censoring: cut the trajectories afterwards

If the censoring time is independent of everything in the process,
drawing it separately and cutting the trajectories is correct and
simple. That is what
[`add_censoring()`](https://adibender.github.io/simpexp/reference/add_censoring.md)
does: it shortens the episode containing the censoring time and drops
everything after it.

``` r

cens <- add_censoring(sim, rate = 0.15, keep_censoring_time = TRUE)
head(cens[c("id", "tstart", "tstop", "status", "censoring_time")], 5)
#>   id tstart        tstop status censoring_time
#> 1  1      0 2.9104729711      0   2.910473e+00
#> 2  2      0 0.9966957476      0   9.966957e-01
#> 3  3      0 0.0004771062      0   4.771062e-04
#> 4  4      0 5.6387400624      0   5.638740e+00
#> 5  5      0 4.7712454998      1   1.257624e+01
mean(cens$status)
#> [1] 0.2175
```

Any distribution will do – pass the times directly, for instance for a
staggered end of follow-up:

``` r

add_censoring(sim, times = runif(nrow(df), 2, 10)) |> (\(d) mean(d$status))()
#> [1] 0.266
```

On a multi-state process the cut falls wherever the censoring time
lands, and the later episodes of that subject disappear:

``` r

msm_sim <- sim_pexp_msm(
  list(
    transition("healthy", "ill",  ~ log(0.4)),
    transition("ill", "healthy",  ~ log(0.6)),
    transition("ill", "dead",     ~ log(0.1))),
  data.frame(id = 1:500), cut = seq(0, 20, by = 0.1))

cens_msm <- add_censoring(msm_sim, rate = 0.1)
n_episodes <- function(d) {
  tb <- table(factor(table(d$id), levels = 1:6))
  setNames(as.integer(tb), paste0("n=", 1:6))
}
rbind(uncensored = n_episodes(msm_sim), censored = n_episodes(cens_msm))
#>            n=1 n=2 n=3 n=4 n=5 n=6
#> uncensored   0  56   1  42  12  62
#> censored    96  91  65  43  37  43
```

## Dependent censoring belongs in the process

If the censoring hazard depends on covariates or on the state of the
process, cutting afterwards is the wrong construction: the censoring
time would have to be drawn in competition with the other transitions,
not independently of them. Express it as what it is – a transition into
a state that ends the trajectory without an event:

``` r

dep <- sim_pexp_msm(
  list(
    transition("healthy", "ill",     ~ log(0.3)),
    transition("healthy", "dropout", ~ log(0.08) + 0.5 * x),
    transition("ill",     "dead",    ~ log(0.15)),
    transition("ill",     "dropout", ~ log(0.20) + 0.5 * x)),
  data.frame(id = 1:3000, x = rnorm(3000)),
  cut = seq(0, 10, by = 0.1),
  censoring_states = "dropout")

table(dep$to, dep$status, useNA = "ifany")
#>          
#>              0    1
#>   healthy    0    0
#>   ill        0 2295
#>   dropout    0 1813
#>   dead       0  865
#>   <NA>     322    0
```

Subjects who reach `dropout` get `status = 0`, like any other censored
subject, but their censoring time came out of the process: it is three
times as likely while ill, and it depends on `x`. That dependence is
visible in the data:

``` r

drop_rate <- function(d, state) {
  e <- d[d$from == state, ]
  sum(e$to == "dropout", na.rm = TRUE) / sum(e$tstop - e$tstart)
}
c(healthy = drop_rate(dep, "healthy"), ill = drop_rate(dep, "ill"))
#>    healthy        ill 
#> 0.08710871 0.19410835
```

The two rates are the covariate-averaged versions of
`0.08 * exp(0.5 * x)` and `0.20 * exp(0.5 * x)`, i.e. about `exp(0.125)`
times the baseline for standard normal `x`.

``` r

c(healthy = 0.08, ill = 0.20) * exp(0.5^2 / 2)
#>    healthy        ill 
#> 0.09065188 0.22662969
```

## Which one to use

- **administrative** – set the end of `cut`; nothing else to do.
- **independent random** –
  [`add_censoring()`](https://adibender.github.io/simpexp/reference/add_censoring.md).
- **depends on covariates or state** – a transition into a
  `censoring_states` state.

Only the third kind changes what a correct estimator must do, which is
exactly why it should not be faked by the first two.
