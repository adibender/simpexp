# simpexp

Simulate single-event, competing-risks and multi-state event histories from
hazards that are **piecewise constant on a time grid**.

Each transition is specified by a formula for its log-hazard. Because the
hazards are piecewise constant by construction, waiting times are drawn by
inverting the cumulative hazard in closed form, and the realised transition is
drawn with probabilities proportional to the transition-specific hazards inside
the interval in which the event falls. No quadrature, no interpolated inverse,
no rejection sampling — the grid is part of the model, not a tuning parameter.

```r
# install.packages("pak")
pak::pak("adibender/simpexp")
```

## Three interfaces, one engine

```r
library(simpexp)
df  <- data.frame(x = rnorm(500))
cut <- seq(0, 10, by = 0.1)

# single event, hazard rising over time
sim_pexp(~ log(0.1) + 0.05 * t + 0.3 * x, df, cut)

# competing risks: one formula per cause
sim_pexp_cr(list(relapse = ~ log(0.1) + 0.3 * x, death = ~ log(0.05)), df, cut)

# multi state, here illness-death with recovery
sim_pexp_msm(
  list(
    transition("healthy", "ill",     ~ log(0.4) + 0.3 * x),
    transition("healthy", "dead",    ~ log(0.02)),
    transition("ill",     "healthy", ~ log(0.8) - 0.25 * t_sojourn),
    transition("ill",     "dead",    ~ log(0.05) + 0.15 * t_sojourn)
  ),
  df, cut)
```

All three return the same long format, one row per episode:

| id | episode | from | to | tstart | tstop | status | transition | … |
|---|---|---|---|---|---|---|---|---|

`status = 0` marks an episode that ended without a transition, with `to` missing.

## Markov, semi-Markov, neither

What a hazard may depend on is what distinguishes the three cases, and it is
simply a matter of which variable the formula uses:

| variable | meaning | process |
|---|---|---|
| `t` (= `tend`) | time since the origin | Markov |
| `t_sojourn` | time in the current state | semi-Markov |
| `t_entry` | time the current state was entered | neither |
| `n_visits_<state>` | visits to a state so far | neither |
| `t_total_<state>` | time spent in a state so far | neither |
| `n_transitions` | transitions so far | neither |

```r
# the hazard of dying rises with every previous illness episode
transition("ill", "dead", ~ log(0.02) + 0.5 * n_visits_ill)
```

`keep_history = TRUE` returns these features with the episodes. They are the
ground truth of the process, which makes them the oracle features of a
simulation study.

## Censoring

Follow-up ends at the last value of `cut`, which censors administratively.
Anything earlier is added afterwards:

```r
sim <- sim_pexp(~ log(0.2) + 0.3 * x, df, cut)
add_censoring(sim, rate = 0.1)              # exponential censoring times
add_censoring(sim, times = runif(500, 2, 10))  # given times, e.g. staggered entry
```

Censoring that depends on covariates or on the state of the process is better
expressed as part of the process itself — a transition into a censoring state,
competing with the others:

```r
sim_pexp_msm(
  list(transition("0", "event", ~ log(0.3)),
       transition("0", "lost",  ~ log(0.1) + 0.5 * x)),
  df, cut, censoring_states = "lost")
```

## Fitting back what was simulated

The hazards are specified under the same convention under which piecewise
exponential additive models are fitted — constant on `(tstart, tend]`,
evaluated at `tend` — so a model fitted to the simulated data estimates the
quantity that was specified. `tests/testthat/test-roundtrip.R` does exactly
that: expand the episodes into intervals, fit a Poisson model with log exposure
as offset, and recover the coefficients.

## Relation to pammtools

The simulation functions of
[pammtools](https://github.com/adibender/pammtools) are the origin of this
package. This package is developed independently of them: it shares no code and
no dependency in either direction, so simulation and model fitting are no longer
tied to each other. pammtools keeps its own `sim_pexp()` for now.
