# Process vertebral measurements and PC scores of 2D or 3D geometric morphometric data

`process_gmPC()` initializes the analysis workflow by processing
vertebral measurements and already computed PC scores from 2D or 3D
geometric morphometric datasets into an object usable by MorphoRegions.
Such processing includes identifying the vertebra indices, measurements,
and PC scores.

## Usage

``` r
process_gmPC(data, pos = "names", pcscores, eigenvals, specimens = NULL)
```

## Arguments

- data:

  a 3D array (\\p \times k \times n\\) of Procrustes aligned landmarks
  where \\p\\ is the number of of landmark points per vertebra, \\k\\ is
  the number of landmark dimensions (2D or 3D), and \\n\\ is the number
  of vertebrae sampled for the specimen.

- pos:

  either character string `"names"` in which case the function will
  extract vertebral position info form the names of the array (using
  `dimnames(data)[[3]]`), or a vector containing numeric values of
  vertebral position. Default is `"names"`.

- pcscores:

  a matrix or data frame containing PC scores of each vertebra.

- eigenvals:

  a numeric vector containing the eigenvalues of each PC axis.

- specimens:

  if `data` contains multiple specimens, a vector of specimen names of
  length equal to the number of vertebrae in `data`.

## Value

A `regions_pco` object, which contains user-provided eigenvectors in the
`scores` component and eigenvalues in the `eigen.val` component. The
original dataset, including landmark coordinates and positional
information, is stored in the `data` attribute.

## Details

Unlike
[`process_measurements()`](https://aagillet.github.io/MorphoRegions/reference/process_measurements.md),
`process_gmPC()` does not fill in missing values and these should have
been removed or replaced by numeric values before running
`process_gmPC()`.

## See also

- [`process_PC()`](https://aagillet.github.io/MorphoRegions/reference/process_PC.md)
  for processing PC scores from traditional morphometric datasets

- [`svdPCO()`](https://aagillet.github.io/MorphoRegions/reference/svdPCO.md)
  for computing principal coordinate axes from processed vertebra data.

- [`plot.regions_pco()`](https://aagillet.github.io/MorphoRegions/reference/plot.regions_pco.md)
  for plotting PCO axes

## Examples

``` r
# Load coordinates of aligned landmarks, PC scores, and
# eigenvalues; vertebra index as names of array elements
# of aligned landmarks
data("seal_gmm")

# Process PC scores and morphological data:
seal_pca <- process_gmPC(data = seal_gmm$coord_gpa,
                         pcscores = seal_gmm$scores,
                         eigenvals = seal_gmm$eigenvals)

# Process multiple datasets; specimen names in 'specimens'
# vector, vertebrae index in 'vertebrae' vector
data("wolves_gmm")

# Process PC scores and geometric morphometric dataset:
wolves_pca <- process_gmPC(data = wolves_gmm$coord_gpa,
                           pos = wolves_gmm$vertebrae,
                           pcscores = wolves_gmm$scores,
                           eigenvals = wolves_gmm$eigenvals,
                           specimens = wolves_gmm$specimens)
#> Warning: `pos` supplied as numeric vector; assuming its order matches the order of
#> `data`.
#> Warning: The row names of `pcscores` do not match `pos`; assuming row order matches
#> `pos`.
```
