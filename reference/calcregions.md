# Fit segmented regression models for all combinations of breakpoints

`calcregions()` enumerates all possible combinations of breakpoints to
fit multivariate segmented regression models. `addregions()` adds models
with additional numbers of regions to the resulting output object.
`ncombos()` computes an upper bound on the number of breakpoint
combinations that will be tested.

## Usage

``` r
calcregions(
  pco,
  scores,
  noregions,
  minvert = 3,
  cont = TRUE,
  exhaus = TRUE,
  includebp = NULL,
  omitbp = NULL,
  ncombos_file_trigger = 1e+07,
  temp_file_dir = tempdir(TRUE),
  cl = NULL,
  verbose = interactive()
)

addregions(
  regions_results,
  noregions,
  exhaus = TRUE,
  ncombos_file_trigger = 1e+07,
  temp_file_dir = tempdir(TRUE),
  cl = NULL,
  verbose = interactive()
)

# S3 method for class 'regions_results'
summary(object, ...)

ncombos(pco, noregions, minvert = 3, includebp = NULL, omitbp = NULL)
```

## Arguments

- pco:

  a `regions_pco` object; the output of a call to
  [`svdPCO()`](https://aagillet.github.io/MorphoRegions/reference/svdPCO.md).

- scores:

  `numeric`; the indices of the PCO scores to use as the outcomes in
  fitting the models (e.g., `1:4` to use the first four scores). Can
  also be the output of a call to
  [`PCOselect()`](https://aagillet.github.io/MorphoRegions/reference/PCOselect.md).

- noregions:

  `numeric`; for `calcregions()`, the maximum number of regions for
  which models are fit (e.g, 4 to request models with 1 to 4 regions);
  for `addregions()`, a vector containing the numbers of regions to add
  (e.g., 5:6 to request models with 5 and 6 regions); for `ncombos()`, a
  vector containing the numbers of regions to check.

- minvert:

  `numeric`; the minimum number of vertebrae allowed in each region.
  Default is 3.

- cont:

  `logical`; whether to fit models that are continuous (`TRUE`) or
  discontinuous (`FALSE`) at the breakpoints. Default is `TRUE`.

- exhaus:

  `logical`; whether to fit all possible models (`TRUE`) or use
  heuristics to reduce the number of models fit (`FALSE`). Default is
  `TRUE`. See Details. Setting to `FALSE` can reduce the size of the
  resulting object.

- includebp:

  an optional vector of vertebrae that must be included in any tested
  set of breakpoints, e.g., if it is known that two regions are divided
  at that vertebra. `includebp` does not have to obey the `minvert`
  rules, but a warning will be thrown if it doesn't.

- omitbp:

  an optional vector of vertebrae to be omitted from the list of
  possible breakpoints, e.g., if it is known that two adjacent vertebrae
  belong to the same region.

- ncombos_file_trigger:

  `numeric`; when the number of eligible combinations of breakpoints
  exceeds this number, the problem will be split into smaller problems,
  with the results of each stored in its own temporary file in the
  directory supplied to `temp_file_dir` before being re-read into
  memory. The primary purpose of this is to preserve memory when
  `exhaus = FALSE` by delegating storage of the results to disk instead
  of RAM.

- temp_file_dir:

  string; the directory where the temporary files will be saved (and
  then deleted) when the number of breakpoint combinations exceeds
  `ncombos_file_trigger`. Default is the directory produced by
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html), but it is much
  safer to provide your own directory, which must already exist on your
  machine. See Details.

- cl:

  a cluster object created by
  [`parallel::makeCluster()`](https://rdrr.io/r/parallel/makeCluster.html),
  an integer to indicate number of child-processes (integer values are
  ignored on Windows) for parallel evaluations, or `"future"` to use a
  future backend. `NULL` (the default) refers to sequential evaluation
  (no parallelization). See
  [`pbapply::pbapply()`](https://peter.solymos.org/pbapply/reference/pbapply.html)
  for details.

- verbose:

  `logical`; whether to print information about the fitting process,
  including a progress bar. Default is `TRUE` when running interactively
  and `FALSE` otherwise.

- regions_results, object:

  a `regions_results` object; the output of a call to `calcregions()` or
  `addregions()`.

- ...:

  ignored.

## Value

A `regions_results` object with the following components:

- `results` - the results of the fitting process for each combination of
  breakpoints

- `stats` - statistics summarizing the fitting process. Use
  [`summary()`](https://rdrr.io/r/base/summary.html) to view this
  information in a clean format.

`ncombos()` returns a numeric vector with the number of breakpoint
combinations for each number of regions (which are stored as the names).

## Details

`calcregions()` enumerates all possible combinations of breakpoints that
satisfy the constraint imposed by `minvert` (i.e., that breakpoints need
to be at least `minvert` vertebrae apart) and fits the segmented
regression models implied by each combination. These are multivariate
regression models with the PCO scores specified by `scores` as the
outcomes. When `cont = TRUE`, these regression models are continuous;
i.e., the regression lines for each region connect at the breakpoints.
Otherwise, the models are discontinuous so that each region has its own
intercept and slope. The models are fit using
[`.lm.fit()`](https://rdrr.io/r/stats/lmfit.html), which efficiently
implements ordinary least squares regression.

When `exhaus = FALSE`, heuristics are used to reduce the number of
models to fit, which can be useful for keeping the size of the resulting
object down by avoiding fitting models corresponding to breakpoint
combinations that yield a poor fit to the data. Only breakpoint
combinations that correspond to the breakpoints of the best fitting
model with a smaller number of regions +/- 3 vertebrae are used, and
only models that have an RSS smaller than half a standard deviation more
the smallest RSS are kept.

`addregions()` should be used on an existing `regions_results` object to
add models with more regions. Internally, it works just the same as
`calcregions()`.

`ncomobs()` computes an upper bound on the number of possible breakpoint
combinations. When `exhaus = FALSE` or `includebp` is specified, the
actual number of combinations will be smaller than that produced by
`ncombos()`.

When the number of possible combinations of breakpoints for a given
number of regions (as computed by `ncombos()`) is larger than
`ncombos_file_trigger`, the problem will be split into smaller problems,
with the results of each stored in temporary files that are deleted when
the function completes. These temporary files will be stored in the
directory supplied to `temp_file_dir`. By default, this is the temporary
directory produced by
[`tempdir()`](https://rdrr.io/r/base/tempfile.html). However, this
directory can be deleted by R at any time without warning, which will
cause the function to crash, so it is a good idea to supply your own
directory that will be preserved. You can use `ncombos()` to check to
see if the number of breakpoint combinations exceeds
`ncombos_file_trigger`.

## See also

[`calcmodel()`](https://aagillet.github.io/MorphoRegions/reference/calcmodel.md)
to fit a segmented regression model for a single set of breakpoints;
[`modelselect()`](https://aagillet.github.io/MorphoRegions/reference/modelselect.md)
to select the best model for each number of regions based on RSS;
[`modelsupport()`](https://aagillet.github.io/MorphoRegions/reference/modelsupport.md)
to compute statistics the describe the support of the best models;
[`calcBPvar()`](https://aagillet.github.io/MorphoRegions/reference/calcBPvar.md)
to compute the variability in the optimal breakpoints.

## Examples

``` r
data("alligator")

alligator_data <- process_measurements(alligator,
                                       pos = "Vertebra")

# Compute PCOs
alligator_PCO <- svdPCO(alligator_data)

# Fit segmented regression models for 1 to 5 regions
# using PCOs 1 to 4 and a continuous model with a
# non-exhaustive search
regionresults <- calcregions(alligator_PCO,
                             scores = 1:4,
                             noregions = 5,
                             minvert = 3,
                             cont = TRUE,
                             exhaus = FALSE)

regionresults
#> A `regions_results` object
#>  - number of PCOs used: 4 
#>  - number of regions: 1, 2, 3, 4, 5 
#>  - model type: continuous 
#>  - min vertebrae per region: 3 
#>  - total models saved: 85 
#> Use `summary()` to examine summaries of the fitting process.

# View model fitting summary
summary(regionresults)
#>  Regions Possible Tested Saved Comp. method Saving method
#>        1        1      1     1   Exhaustive           All
#>        2       17     17    17   Exhaustive           All
#>        3      105    105    17   Non-exhaus          SD/2
#>        4      286    237    24   Non-exhaus          SD/2
#>        5      330    323    26   Non-exhaus          SD/2

# Add additional regions to existing results,
# exhaustive search this time
regionresults <- addregions(regionresults,
                            noregions = 6:7,
                            exhaus = TRUE)

regionresults
#> A `regions_results` object
#>  - number of PCOs used: 4 
#>  - number of regions: 1, 2, 3, 4, 5, 6, 7 
#>  - model type: continuous 
#>  - min vertebrae per region: 3 
#>  - total models saved: 218 
#> Use `summary()` to examine summaries of the fitting process.

summary(regionresults)
#>  Regions Possible Tested Saved Comp. method Saving method
#>        1        1      1     1   Exhaustive           All
#>        2       17     17    17   Exhaustive           All
#>        3      105    105    17   Non-exhaus          SD/2
#>        4      286    237    24   Non-exhaus          SD/2
#>        5      330    323    26   Non-exhaus          SD/2
#>        6      126    126   126   Exhaustive           All
#>        7        7      7     7   Exhaustive           All

# Fit segmented regression models for 1 to 5 regions
# using PCOs 1 to 4 and a discontinuous model with a
# exhaustive search, excluding breakpoints at vertebrae
# 10 and 15
regionresults <- calcregions(alligator_PCO,
                             scores = 1:4,
                             noregions = 5,
                             minvert = 3,
                             cont = FALSE,
                             omitbp = c(10, 15))

regionresults
#> A `regions_results` object
#>  - number of PCOs used: 4 
#>  - number of regions: 1, 2, 3, 4, 5 
#>  - model type: discontinuous 
#>  - min vertebrae per region: 3 
#>  - omitted breakpoints: 10, 15 
#>  - total models saved: 385 
#> Use `summary()` to examine summaries of the fitting process.

summary(regionresults)
#>  Regions Possible Tested Saved Comp. method Saving method
#>        1        1      1     1   Exhaustive           All
#>        2       15     15    15   Exhaustive           All
#>        3       78     78    78   Exhaustive           All
#>        4      165    165   165   Exhaustive           All
#>        5      126    126   126   Exhaustive           All

# Compute the number of breakpoint combinations for given
# specification using `ncombos()`; if any number exceeds
# the value supplied to `ncombos_file_trigger`, results
# will temporary be stored in files before being read in and
# deleted.
ncombos(alligator_PCO,
        noregions = 1:8,
        minvert = 3)
#>   1   2   3   4   5   6   7   8 
#>   1  17 105 286 330 126   7   0 
```
