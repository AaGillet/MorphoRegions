# Process vertebral measurements and PC scores of traditional morphometric data

`process_PC()` initializes the analysis workflow by processing a dataset
of vertebral measurements and already computed PC scores into an object
usable by MorphoRegions. Such processing includes identifying the
vertebra indices, measurements, and PC scores.

## Usage

``` r
process_PC(data, pos = 1L, pcscores, eigenvals, posPC)
```

## Arguments

- data:

  a data frame containing a column of vertebra indices and measurements
  for each vertebra, or a named list thereof for multiple specimens with
  the names of the list corresponding to the unique identifiers of each
  specimen.

- pos:

  the name or index of the variable in `data` containing the vertebra
  indices. Default is to use the first column.

- pcscores:

  a matrix or data frame containing PC scores of each vertebra. Note
  that for multiple specimens, an ordination method should have been
  performed on the concatenated data across all specimens so that
  `pcscores` is a single data frame or matrix containing scores for all
  specimens.

- eigenvals:

  a numeric vector containing the eigenvalues of each PC axis.

- posPC:

  a vector containing the positional information of vertebrae in
  `pcscores`, or a named list thereof for multiple specimens with the
  names of the list corresponding to the unique identifiers of each
  specimen. If not provided, will assume the row order of `pcscores`
  matches the vertebra order in `pos`.

## Value

A `regions_pco` object, which contains user-provided eigenvectors in the
`scores` component and eigenvalues in the `eigen.val` component. The
original dataset, including positional information, is stored in the
`data` attribute.

## Details

Unlike
[`process_measurements()`](https://aagillet.github.io/MorphoRegions/reference/process_measurements.md),
`process_PC()` does not fill in missing values and these should have
been removed or replaced by numeric values before running
`process_PC()`.

## See also

- [`process_gmPC()`](https://aagillet.github.io/MorphoRegions/reference/process_gmPC.md)
  for processing PC scores from 2D or 3D geometric morphometric datasets

- [`svdPCO()`](https://aagillet.github.io/MorphoRegions/reference/svdPCO.md)
  for computing principal coordinate axes from processed vertebra data.

- [`plot.regions_pco()`](https://aagillet.github.io/MorphoRegions/reference/plot.regions_pco.md)
  for plotting PCO axes

## Examples

``` r
# Load dataset; vertebral index in first column
# ('Vertebra' column)
data("dolphin")

# Compute PC scores with prcomp:
PCA <- prcomp(dolphin[-1], scale = TRUE)

# Extract PC scores and eigenvalues:
PCA_scores <- as.data.frame(PCA$x)
PCA_eigenval <- PCA$sdev^2


# Process PC scores and morphological data:
pco_dolphin <- process_PC(data = dolphin,
                          pos = "Vertebra",
                          pcscores = PCA_scores,
                          eigenvals = PCA_eigenval,
                          posPC = dolphin[[1]])
```
