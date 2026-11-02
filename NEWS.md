# simpexp 0.0.0.9000

* First version. `sim_pexp()`, `sim_pexp_cr()` and `sim_pexp_msm()` simulate
  single-event, competing-risks and multi-state histories from hazards that are
  piecewise constant on a time grid, `add_censoring()` adds censoring on top.
* The hazard of each transition is given as a formula and may use the time
  scales and history features the simulator maintains, so Markov, semi-Markov
  and history-dependent processes are all expressed the same way.
