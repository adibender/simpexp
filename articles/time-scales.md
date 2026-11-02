# Time scales: Markov, semi-Markov, non-Markov

``` r

library(simpexp)
```

## What a formula may depend on

The right-hand side of a transition formula is evaluated on the
interval-expanded data. Besides the covariates of `data`, the simulator
provides the variables listed in
[`?time_scales`](https://adibender.github.io/simpexp/reference/time_scales.md).
Which of them a hazard uses decides the probabilistic structure of the
process – nothing else does:

| variable | meaning | resulting process |
|----|----|----|
| `t` (`= tend`) | time since the origin of the process | Markov |
| `t_sojourn` | time already spent in the current state | semi-Markov |
| `t_entry` | time at which the current state was entered | non-Markov |
| `n_visits_<state>` | visits to `<state>` so far | non-Markov |
| `t_total_<state>` | time spent in `<state>` so far | non-Markov |
| `n_transitions` | transitions made so far | non-Markov |
| `interval` | index of the current interval of `cut` | – |

There is no switch to set: writing `t_sojourn` into a formula *is* the
specification of a semi-Markov process.

## Markov: the clock never resets

``` r

df <- data.frame(id = 1:2000)

markov <- sim_pexp_msm(
  list(
    transition("healthy", "ill",  ~ log(0.3)),
    transition("ill", "healthy",  ~ log(0.5)),
    transition("ill", "dead",     ~ log(0.05) + 0.1 * t)),
  df, cut = seq(0, 20, by = 0.1))
```

The hazard of dying while ill grows with calendar time, identically for
a subject who has just fallen ill and for one who has been ill for
years.

## Semi-Markov: the clock resets on entry

``` r

semi <- sim_pexp_msm(
  list(
    transition("healthy", "ill",  ~ log(0.3)),
    transition("ill", "healthy",  ~ log(0.5)),
    transition("ill", "dead",     ~ log(0.05) + 0.4 * t_sojourn)),
  df, cut = seq(0, 20, by = 0.1))
```

Now the risk depends on how long *this* episode of illness has lasted.
The two processes have the same transition structure and differ only in
that one variable; their sojourn-time distributions are visibly
different.

``` r

soj <- function(s) {
  e <- s[s$from == "ill", ]
  e$tstop - e$tstart
}
rbind(markov = summary(soj(markov)), `semi-Markov` = summary(soj(semi)))
#>                     Min.   1st Qu.   Median     Mean  3rd Qu.     Max.
#> markov      0.0001061092 0.4324054 1.066274 1.513517 2.132567 12.16489
#> semi-Markov 0.0003750519 0.4750897 1.151195 1.533443 2.207250  8.73657
```

Because the semi-Markov hazard rises steeply within an episode, long
illness episodes are cut short – the sojourn distribution has a much
lighter tail.

## Non-Markov: history dependence

A process whose hazards look further back than the current state and the
time spent in it is non-Markov. The history features make that a
one-line change.

``` r

relapsing <- sim_pexp_msm(
  list(
    transition("healthy", "ill",  ~ log(0.3) + 0.25 * n_visits_ill),
    transition("ill", "healthy",  ~ log(0.5)),
    transition("ill", "dead",     ~ log(0.03) + 0.1 * t_total_ill)),
  df, cut = seq(0, 20, by = 0.1), keep_history = TRUE)
```

Each episode of illness makes the next one more likely, and the
accumulated time spent ill raises the risk of dying.
`keep_history = TRUE` returns the feature values the hazards were
evaluated with, which is the ground truth of the data-generating
process:

``` r

head(relapsing[c("id", "episode", "from", "to", "tstart", "tstop",
                 "n_visits_ill", "t_total_ill")], 8)
#>   id episode    from      to    tstart     tstop n_visits_ill t_total_ill
#> 1  1       1 healthy     ill  0.000000  1.678785            0   0.0000000
#> 2  1       2     ill healthy  1.678785  2.326540            1   0.0000000
#> 3  1       3 healthy     ill  2.326540  5.291553            1   0.6477546
#> 4  1       4     ill healthy  5.291553  8.104056            2   0.6477546
#> 5  1       5 healthy     ill  8.104056  8.479224            2   3.4602572
#> 6  1       6     ill healthy  8.479224  8.578599            3   3.4602572
#> 7  1       7 healthy     ill  8.578599 10.319620            3   3.5596318
#> 8  1       8     ill healthy 10.319620 13.654427            4   3.5596318
```

That the mechanism really bites can be read off the data: the rate of
becoming ill again, by the number of previous illness episodes.

``` r

healthy <- relapsing[relapsing$from == "healthy", ]
healthy$k <- healthy$n_visits_ill
tab <- do.call(rbind, lapply(sort(unique(healthy$k)), function(k) {
  d <- healthy[healthy$k == k, ]
  data.frame(n_visits_ill = k, episodes = nrow(d),
    rate = sum(d$to == "ill", na.rm = TRUE) / sum(d$tstop - d$tstart),
    implied = 0.3 * exp(0.25 * k))
}))
print(tab[tab$episodes > 30, ], digits = 3, row.names = FALSE)
#>  n_visits_ill episodes  rate implied
#>             0     2000 0.300   0.300
#>             1     1856 0.387   0.385
#>             2     1667 0.506   0.495
#>             3     1387 0.653   0.635
#>             4     1076 0.771   0.815
#>             5      750 1.120   1.047
#>             6      487 1.375   1.345
#>             7      274 2.078   1.726
#>             8      140 2.286   2.217
#>             9       64 3.148   2.846
#>            10       32 3.573   3.655
```

The empirical occurrence/exposure rates track `0.3 * exp(0.25 * k)`,
which is what the formula asked for.

## Why oracle features matter

In a simulation study the interesting question is usually whether a
method can *recover* a dependence of this kind from the observed
trajectories. Having the true feature values alongside the data
separates two failure modes: a model that cannot represent the
dependence, and a model that can but does not find it. Fitting with the
oracle features gives the ceiling against which a data-driven model is
measured.
