# Evaluate model support

`modelsupport()` computes measures of the relative support of each of
the best models identified by
[`modelselect()`](https://aagillet.github.io/MorphoRegions/reference/modelselect.md)
to facilitate selecting the optimal number and position of regions.
These measures are in the form of information criteria (AICc and BIC).

## Usage

``` r
modelsupport(models)
```

## Arguments

- models:

  a `regions_modelselect` object; the output of a call to
  [`modelselect()`](https://aagillet.github.io/MorphoRegions/reference/modelselect.md).

## Value

A `regions_modelsupport` object, which contains the best model for each
number of regions as determined by the AICc and BIC. The computed
statistics are `AICc`/`BIC`–the value of the information criterion (IC)
for each model, `deltaAIC`/`deltaBIC`–the difference between the IC for
the corresponding model and that of the model with the lowest IC value,
`model_lik`–the likelihood ratio of the model against the model with the
lowest IC value, and `Ak_weight`/`BIC_weight`–the Akaike weights for
each model used to compute the region score. The region score is a
weighted average of the numbers of regions, weighted by the Akaike
weights to represent the variability around the optimal number of
regions.

## See also

[`modelselect()`](https://aagillet.github.io/MorphoRegions/reference/modelselect.md),
[`calcregions()`](https://aagillet.github.io/MorphoRegions/reference/calcregions.md),
[`calcBPvar()`](https://aagillet.github.io/MorphoRegions/reference/calcBPvar.md),
[`modelperf()`](https://aagillet.github.io/MorphoRegions/reference/modelperf.md),
[`plotsegreg()`](https://aagillet.github.io/MorphoRegions/reference/plotsegreg.md)

## Examples

``` r
data("alligator")

alligator_data <- process_measurements(alligator,
                                       pos = "Vertebra")

# Compute PCOs
alligator_PCO <- svdPCO(alligator_data)

# Fit segmented regression models for 1 to 7 regions
# using PCOs 1 to 4 and a continuous model with a
# non-exhaustive search
regionresults <- calcregions(alligator_PCO,
                             scores = 1:4,
                             noregions = 7,
                             minvert = 3,
                             cont = TRUE,
                             exhaus = FALSE,
                             verbose = FALSE)

regionresults
#> A `regions_results` object
#>  - number of PCOs used: 4 
#>  - number of regions: 1, 2, 3, 4, 5, 6, 7 
#>  - model type: continuous 
#>  - min vertebrae per region: 3 
#>  - total models saved: 112 
#> Use `summary()` to examine summaries of the fitting process.

# For each number of regions, identify best
# model based on minimizing RSS
bestresults <- modelselect(regionresults)
bestresults
#>  Regions BP 1 BP 2 BP 3 BP 4 BP 5 BP 6 sumRSS RSS.1 RSS.2 RSS.3 RSS.4
#>        1    .    .    .    .    .    .  0.725 0.221 0.327 0.135 0.041
#>        2   12    .    .    .    .    .  0.356 0.077 0.127 0.112 0.040
#>        3    9   14    .    .    .    .  0.154 0.019 0.039 0.063 0.033
#>        4    9   13   19    .    .    .  0.098 0.008 0.010 0.045 0.035
#>        5    6    9   13   19    .    .  0.054 0.008 0.010 0.017 0.020
#>        6    6    9   12   15   19    .  0.042 0.007 0.009 0.017 0.009
#>        7    6    9   12   15   18   21  0.042 0.006 0.011 0.014 0.010

# Evaluate support for each model and rank models
supp <- modelsupport(bestresults)
supp
#> - Model support (AICc)
#>  Regions BP 1 BP 2 BP 3 BP 4 BP 5 BP 6 sumRSS     AICc deltaAIC model_lik
#>        5    6    9   13   19    .    .  0.054 -566.603    0.000     1.000
#>        6    6    9   12   15   19    .  0.042 -565.883    0.720     0.698
#>        7    6    9   12   15   18   21  0.042 -537.487   29.116     0.000
#>        4    9   13   19    .    .    .  0.098 -534.934   31.668     0.000
#>        3    9   14    .    .    .    .  0.154 -512.839   53.764     0.000
#>        2   12    .    .    .    .    .  0.356 -454.040  112.563     0.000
#>        1    .    .    .    .    .    .  0.725 -404.530  162.073     0.000
#>  Ak_weight
#>      0.589
#>      0.411
#>      0.000
#>      0.000
#>      0.000
#>      0.000
#>      0.000
#> Region score: 5.41 
#> 
#> - Model support (BIC)
#>  Regions BP 1 BP 2 BP 3 BP 4 BP 5 BP 6 sumRSS      BIC deltaBIC model_lik
#>        6    6    9   12   15   19    .  0.042 -525.687    0.000      1.00
#>        5    6    9   13   19    .    .  0.054 -524.763    0.924      0.63
#>        7    6    9   12   15   18   21  0.042 -503.838   21.849      0.00
#>        4    9   13   19    .    .    .  0.098 -495.206   30.481      0.00
#>        3    9   14    .    .    .    .  0.154 -478.160   47.527      0.00
#>        2   12    .    .    .    .    .  0.356 -426.753   98.933      0.00
#>        1    .    .    .    .    .    .  0.725 -386.534  139.152      0.00
#>  BIC_weight
#>       0.613
#>       0.387
#>       0.000
#>       0.000
#>       0.000
#>       0.000
#>       0.000
#> Region score: 5.61 

# 5 regions best based on AICc; 6 regions based on BIC
```
