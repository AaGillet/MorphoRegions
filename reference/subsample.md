# Subsample a dataset

`subsample()` creates a smaller version of the original dataset by
sampling its rows. Because PCOs should be computed on the full dataset
and most other functions take in `regions_pco` objects, `subsample()`
requires a `regions_pco` object as its input.

## Usage

``` r
subsample(pco, sample = NULL, type = "seq")
```

## Arguments

- pco:

  a `regions_pco` object; the output of a call to
  [`svdPCO()`](https://aagillet.github.io/MorphoRegions/reference/svdPCO.md).

- sample:

  `numeric`; either the number or proportion of vertebrae to remain the
  sampled dataset. If `NULL`, the original dataset is returned.

- type:

  string; the type of subsampling to do, either `"seq"` for sampling in
  sequence or `"random"` for random sampling. Default is `"seq"`.
  Abbreviations allowed.

## Value

A `regions_pco` object, a subset of the original supplied to `pco`. The
original dataset is stored as an attribute, which itself contains the
subsampling indices.

## See also

[`svdPCO()`](https://aagillet.github.io/MorphoRegions/reference/svdPCO.md),
[`process_measurements()`](https://aagillet.github.io/MorphoRegions/reference/process_measurements.md),
[`plotvertmap()`](https://aagillet.github.io/MorphoRegions/reference/plotvertmap.md)
to visualize the vertebral map after subsampling.

## Examples

``` r
data("alligator")

alligator_data <- process_measurements(alligator,
                                       pos = "Vertebra")

# Compute PCOs
alligator_PCO <- svdPCO(alligator_data)

# Plot vertebrae before subsampling
plotvertmap(alligator_PCO, dropNA = FALSE,
            text = TRUE)


# Subsample data after estimating PCOs; subsample down
# to 15 vertebrae
alligator_PCO_sample <- subsample(alligator_PCO,
                                  sample = 15)

# Plot vertebrae after subsampling; gray vertebrae
# have been dropped
plotvertmap(alligator_PCO_sample, dropNA = FALSE,
            text = TRUE)
```
