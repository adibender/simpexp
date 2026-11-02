# Package index

## Simulate

The three interfaces to the simulation engine. They differ only in how
the process is specified; the output format is the same.

- [`sim_pexp()`](https://adibender.github.io/simpexp/reference/sim_pexp.md)
  [`sim_pexp_cr()`](https://adibender.github.io/simpexp/reference/sim_pexp.md)
  [`sim_pexp_msm()`](https://adibender.github.io/simpexp/reference/sim_pexp.md)
  : Simulate event histories from piecewise exponential hazards

## Specify

Building blocks of a process specification: a transition, and the
variables a hazard formula may depend on.

- [`transition()`](https://adibender.github.io/simpexp/reference/transition.md)
  : Specify one transition of an event-history process
- [`time_scales`](https://adibender.github.io/simpexp/reference/time_scales.md)
  : Time scales and history features available in a hazard formula

## Censor

- [`add_censoring()`](https://adibender.github.io/simpexp/reference/add_censoring.md)
  : Censor simulated event histories

## Package

- [`simpexp`](https://adibender.github.io/simpexp/reference/simpexp-package.md)
  [`simpexp-package`](https://adibender.github.io/simpexp/reference/simpexp-package.md)
  : simpexp: Simulate Event Histories from Piecewise Exponential Hazards
