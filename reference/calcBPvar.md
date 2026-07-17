# Calculate variability of breakpoints

`calcBPvar()` computes an estimate of the variability of the breakpoints
for a given number of regions. This involves computing the weighted mean
and standard deviation of each breakpoint using Akaike weights.

## Usage

``` r
calcBPvar(regions_results, noregions, pct = 0.05, criterion = "aic")
```

## Arguments

- regions_results:

  a `regions_results` object; the output of a call to
  [`calcregions()`](https://aagillet.github.io/MorphoRegions/reference/calcregions.md)
  or
  [`addregions()`](https://aagillet.github.io/MorphoRegions/reference/calcregions.md).

- noregions:

  the number of regions for which the weighted mean and standard
  deviation are to be computed.

- pct:

  the proportion of best model to keep from the original total number of
  possible models

- criterion:

  string; the criterion used to compute the weights. Allowable options
  include `"aic"` and `"bic"`. Abbreviations allowed.

## Value

A `regions_BPvar` object, which has two components:

- `WeightedBP` is a matrix containing the weighted mean and standard
  deviation of each breakpoint

- `BestModels` is a data frame containing the models used to compute the
  weighted breakpoint statistics and the weights each one is given.

## See also

[`calcregions()`](https://aagillet.github.io/MorphoRegions/reference/calcregions.md)
for fitting segmented regression models to all combinations of
breakpoints.

## Examples

``` r
data("alligator")

alligator_data <- process_measurements(alligator,
                                       pos = "Vertebra")

# Compute PCOs
alligator_PCO <- svdPCO(alligator_data)

# Fit segmented regression models for 1 to 7 regions
# using PCOs 1 to 4 and a continuous model with a
# exhaustive search
regionresults <- calcregions(alligator_PCO,
                             scores = 1:4,
                             noregions = 7,
                             minvert = 3,
                             cont = TRUE,
                             exhaus = TRUE,
                             verbose = FALSE)

# Compute Akaike-weighted location and SD of optimal
# breakpoints using top 10% of models with 4 regions
calcBPvar(regionresults, noregions = 4,
          pct = .1, criterion = "aic")
#>        BP 1   BP 2   BP 3
#> wMean 8.904 12.808 19.506
#> wSD   0.441  0.596  0.976
#> - Computed using top 10.14% of models
```
