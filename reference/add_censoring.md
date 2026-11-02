# Censor simulated event histories

The simulation itself censors administratively at the end of the time
grid. `add_censoring()` adds censoring that happens earlier: it cuts
every trajectory at its censoring time, shortens the episode that
contains it and drops everything after it.

## Usage

``` r
add_censoring(data, times = NULL, rate = 0.1, keep_censoring_time = FALSE)
```

## Arguments

- data:

  Simulated event histories, as returned by
  [`sim_pexp()`](https://adibender.github.io/simpexp/reference/sim_pexp.md),
  [`sim_pexp_cr()`](https://adibender.github.io/simpexp/reference/sim_pexp.md)
  or
  [`sim_pexp_msm()`](https://adibender.github.io/simpexp/reference/sim_pexp.md).

- times:

  Either a numeric vector of censoring times, named by `id` or in the
  order of the subjects, or a function of the number of subjects
  returning such a vector. Defaults to exponential censoring times, see
  `rate`.

- rate:

  Rate of the exponential censoring distribution, used when `times` is
  not given.

- keep_censoring_time:

  Append the censoring time of each subject as the column
  `censoring_time`.

## Value

`data` with trajectories cut at the censoring times. Episodes that
contain the censoring time get `status = 0`, `to = NA` and `tstop` equal
to the censoring time; later episodes are removed. Subjects censored
before their first episode begins drop out entirely.

## Details

Censoring that depends on covariates or on the state of the process does
not need this function: give the process a transition into a censoring
state and name that state in `censoring_states` of
[`sim_pexp_msm()`](https://adibender.github.io/simpexp/reference/sim_pexp.md).
The engine then draws the censoring time from its own hazard, in
competition with the other transitions, which is the natural way to
express dependent censoring – *dependent* in the sense of covariates,
not of the unobserved future.

## Examples

``` r
df  <- data.frame(x = rnorm(100))
sim <- sim_pexp(~ log(0.2) + 0.3 * x, df, cut = seq(0, 10, by = 0.1))
add_censoring(sim, rate = 0.1) |> head()
#>   id episode from   to tstart     tstop status transition            x
#> 1  1       1    0    1      0 8.2974118      1       0->1 -1.400043517
#> 2  2       1    0    1      0 0.1030237      1       0->1  0.255317055
#> 3  3       1    0 <NA>      0 1.0312832      0       <NA> -2.437263611
#> 4  4       1    0    1      0 0.7301007      1       0->1 -0.005571287
#> 5  5       1    0 <NA>      0 0.3253287      0       <NA>  0.621552721
#> 6  6       1    0    1      0 0.6014634      1       0->1  1.148411606

# fixed end of follow-up per subject, e.g. staggered study entry
add_censoring(sim, times = runif(100, 2, 10)) |> head()
#>   id episode from   to tstart     tstop status transition            x
#> 1  1       1    0 <NA>      0 7.9441459      0       <NA> -1.400043517
#> 2  2       1    0    1      0 0.1030237      1       0->1  0.255317055
#> 3  3       1    0 <NA>      0 5.9746020      0       <NA> -2.437263611
#> 4  4       1    0    1      0 0.7301007      1       0->1 -0.005571287
#> 5  5       1    0 <NA>      0 9.6633530      0       <NA>  0.621552721
#> 6  6       1    0    1      0 0.6014634      1       0->1  1.148411606
```
