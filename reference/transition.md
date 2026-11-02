# Specify one transition of an event-history process

A transition is the triple (origin state, target state, log-hazard). The
log-hazard is given as a one-sided formula whose right-hand side is
evaluated on the interval-expanded data, so it may use covariates of
`data` as well as the time scales and history features that the
simulator provides (see
[time_scales](https://adibender.github.io/simpexp/reference/time_scales.md)).

## Usage

``` r
transition(from, to, formula)
```

## Arguments

- from:

  Origin state, a single string or number.

- to:

  Target state, a single string or number.

- formula:

  One-sided formula for the log-hazard of this transition, e.g.
  `~ log(0.1) + 0.3 * x + 0.05 * t_sojourn`.

## Value

An object of class `pexp_transition`.

## See also

[`sim_pexp_msm()`](https://adibender.github.io/simpexp/reference/sim_pexp.md),
[time_scales](https://adibender.github.io/simpexp/reference/time_scales.md)

## Examples

``` r
transition(0, 1, ~ log(0.5) + 0.2 * x)
#> <pexp_transition> 0 -> 1: log(0.5) + 0.2 * x
```
