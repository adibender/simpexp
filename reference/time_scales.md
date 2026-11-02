# Time scales and history features available in a hazard formula

The right-hand side of a transition formula is evaluated on the
interval-expanded data. Besides the columns of `data`, the following
variables are available.

## Time scales

- `tstart`, `tend`:

  Boundaries of the current interval on the process time scale. The
  hazard is constant on `(tstart, tend]`.

- `t`:

  Alias of `tend`: time since the origin of the process ("clock
  forward"). The hazard of an interval is evaluated at its end, the same
  convention under which piecewise exponential models are fitted, so a
  hazard specified here and one estimated from the simulated data refer
  to the same time point.

- `t_entry`:

  Time at which the current state was entered.

- `t_sojourn`:

  Time spent in the current state, `tend - t_entry` ("clock reset").

- `interval`:

  Index of the current interval in `cut`, so that a hazard given as one
  value per interval can be used directly:
  `~ log(rates[interval]) + 0.3 * x`.

## History features

- `n_transitions`:

  Number of transitions made before the current one.

- `n_visits_<state>`:

  Number of visits to `<state>` so far, including the current state.

- `t_total_<state>`:

  Time spent in `<state>` before the current interval.

A process whose hazards use only `t` is Markov; one that uses
`t_sojourn` is semi-Markov; one that uses the history features is
neither.
