# simpexp: Simulate Event Histories from Piecewise Exponential Hazards

Simulates single-event, competing-risks and multi-state event histories
from hazards that are piecewise constant on a user-supplied time grid.
Each transition is specified by a formula for its log-hazard, which may
use the time scales and history features that the simulator maintains,
so that non-Markov processes can be specified directly. Because the
hazards are piecewise constant by construction, waiting times are drawn
by exact inversion of the cumulative hazard rather than by numerical
approximation.

## See also

Useful links:

- <https://adibender.github.io/simpexp/>

- <https://github.com/adibender/simpexp>

- Report bugs at <https://github.com/adibender/simpexp/issues>

## Author

**Maintainer**: Andreas Bender <bender.at.R@gmail.com>
([ORCID](https://orcid.org/0000-0001-5628-8611))

Authors:

- Andreas Bender <bender.at.R@gmail.com>
  ([ORCID](https://orcid.org/0000-0001-5628-8611))
