# Geometric morphometric data and PC scores of the vertebral column two wolves

Procrustes aligned 3d landmarks (`coord_gpa`) and PC scores (`scores`)
from 10 subsampled vertebrae along the backbone of two specimens of
Canis lupus (MU 071 and MU 073).

## Usage

``` r
data("wolves_gmm")
```

## Format

A list containing: a 3D array of Procrustes aligned 3D landmarks
(`coord_gpa`) of 10 subsampled vertebrae per specimen, and 30 landmarks
per vertebra, a matrix of PC scores from a PCA performed on Procrustes
aligned landmarks (`scores`), a numeric vector with the eigenvalues from
the PCA (`eigenvals`), a numeric vector with the position information of
vertebrae (`vertebrae`), and a vector with the names of the specimens
(`specimens`).

## References

Schwab, J., Figueirido, B., & Jones, K. E. (2026). Ecological inference
from isolated vertebrae: Evaluating functional signal across the
carnivoran spine. Journal of Morphology, 287(1), e70109.
