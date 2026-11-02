# Simulate event histories from piecewise exponential hazards

`sim_pexp()` simulates single-event data, `sim_pexp_cr()` competing
risks and `sim_pexp_msm()` general multi-state processes, including
recurrent states. All three share one engine and one output format; they
differ only in how the process is specified.

## Usage

``` r
sim_pexp(formula, data, cut, origin = "0", keep_history = FALSE)

sim_pexp_cr(formulas, data, cut, origin = "0", keep_history = FALSE)

sim_pexp_msm(
  transitions,
  data,
  cut,
  terminal_states = NULL,
  censoring_states = character(),
  max_transitions = 100L,
  keep_history = FALSE
)
```

## Arguments

- formula:

  One-sided formula for the log-hazard (single event).

- data:

  A data frame of subjects, one row each. A column `id` is added if
  absent. Optional columns `from` (initial state) and `t_entry` (time of
  entry into that state) allow starting subjects in different states or
  at different times.

- cut:

  Time grid on which the hazards are piecewise constant. Must start at
  0; its last value is the end of follow-up and therefore the time of
  administrative censoring.

- origin:

  Label of the initial state (single event and competing risks).

- keep_history:

  Append the history features that the hazards were evaluated with (see
  [time_scales](https://adibender.github.io/simpexp/reference/time_scales.md))
  to the returned episodes. They are the ground truth of the
  data-generating process and are what a model would have to recover, so
  they are useful as oracle features in a simulation study.

- formulas:

  A list of one-sided formulas, one per competing cause. Names are used
  as the cause labels, otherwise causes are numbered.

- transitions:

  A list of
  [`transition()`](https://adibender.github.io/simpexp/reference/transition.md)
  objects, or of `list(from, to, formula)`.

- terminal_states:

  States that end a trajectory. Defaults to all states that have no
  outgoing transition.

- censoring_states:

  States that end a trajectory with `status = 0`, i.e. that encode
  censoring as a competing event. This is how censoring that depends on
  covariates or on the state of the process is expressed; see
  [`add_censoring()`](https://adibender.github.io/simpexp/reference/add_censoring.md)
  for censoring applied after the fact.

- max_transitions:

  Guard against processes that never reach a terminal state; a warning
  is issued if it is hit.

## Value

A data frame with one row per episode and the columns `id`, `episode`,
`from`, `to`, `tstart`, `tstop`, `status`, `transition`, plus the
columns of `data`. `status = 0` marks an episode that ended without a
transition, i.e. censoring, and then `to` is `NA`.

## Details

Hazards are piecewise constant on the grid `cut`, which makes the
simulation exact: waiting times are drawn by inverting the cumulative
hazard in closed form, and when several transitions compete, the
realised one is drawn with probabilities proportional to the
transition-specific hazards within the interval in which the event
falls. The grid is therefore part of the model, not a numerical tuning
parameter: a hazard that varies quickly needs a fine grid to be
represented faithfully.

## See also

[time_scales](https://adibender.github.io/simpexp/reference/time_scales.md)
for what a hazard formula may depend on,
[`add_censoring()`](https://adibender.github.io/simpexp/reference/add_censoring.md)
for censoring that is not administrative.

## Examples

``` r
df <- data.frame(x = rnorm(50))

# single event, hazard increasing in time
sim_pexp(~ log(0.1) + 0.05 * t + 0.3 * x, df, cut = seq(0, 10, by = 0.1)) |> head()
#>   id episode from   to tstart     tstop status transition          x
#> 1  1       1    0    1      0  2.451034      1       0->1 -1.4726915
#> 2  2       1    0    1      0  2.306163      1       0->1  0.1236927
#> 3  3       1    0 <NA>      0 10.000000      0       <NA> -1.0955924
#> 4  4       1    0    1      0  5.914207      1       0->1  0.3883423
#> 5  5       1    0 <NA>      0 10.000000      0       <NA> -0.7550303
#> 6  6       1    0    1      0  2.445891      1       0->1 -1.2316026

# competing risks
sim_pexp_cr(list(relapse = ~ log(0.1) + 0.3 * x, death = ~ log(0.05)),
  df, cut = seq(0, 10, by = 0.1)) |> head()
#>   id episode from      to tstart     tstop status transition          x
#> 1  1       1    0   death      0  1.220442      1   0->death -1.4726915
#> 2  2       1    0   death      0  2.555442      1   0->death  0.1236927
#> 3  3       1    0   death      0  3.291608      1   0->death -1.0955924
#> 4  4       1    0    <NA>      0 10.000000      0       <NA>  0.3883423
#> 5  5       1    0 relapse      0  8.370190      1 0->relapse -0.7550303
#> 6  6       1    0 relapse      0  1.481785      1 0->relapse -1.2316026

# illness-death with recovery; the hazard of dying while ill grows with the
# time already spent being ill, which makes the process non-Markov
sim_pexp_msm(
  list(
    transition("healthy", "ill",     ~ log(0.4) + 0.3 * x),
    transition("healthy", "dead",    ~ log(0.02)),
    transition("ill",     "healthy", ~ log(0.8) - 0.25 * t_sojourn),
    transition("ill",     "dead",    ~ log(0.05) + 0.15 * t_sojourn)
  ),
  df, cut = seq(0, 10, by = 0.1)) |> head()
#>   id episode    from      to   tstart     tstop status    transition          x
#> 1  1       1 healthy    dead 0.000000  2.001132      1 healthy->dead -1.4726915
#> 2  2       1 healthy     ill 0.000000  5.834358      1  healthy->ill  0.1236927
#> 3  2       2     ill healthy 5.834358  8.166126      1  ill->healthy  0.1236927
#> 4  2       3 healthy    <NA> 8.166126 10.000000      0          <NA>  0.1236927
#> 5  3       1 healthy     ill 0.000000  1.575736      1  healthy->ill -1.0955924
#> 6  3       2     ill healthy 1.575736  2.576441      1  ill->healthy -1.0955924
```
