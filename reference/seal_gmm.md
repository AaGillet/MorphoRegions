# Geometric morphometric data and PC scores of the vertebral column a seal

Procrustes aligned 3d landmarks (`coord_gpa`) and PC scores (`scores`)
of the backbone of Phoca vitulina (ZMUC 160).

## Usage

``` r
data("seal_gmm")
```

## Format

A list containing: a 3D array of Procrustes aligned 3D landmarks
(`coord_gpa`) of 25 presacral vertebrae with 40 landmarks per vertebra,
a matrix of PC scores from a PCA performed on Procrustes aligned
landmarks (`scores`), and a numeric vector with the eigenvalues from the
PCA (`eigenvals`).

## References

Esteban, J. M., Martín-Serra, A., Pérez-Ramos, A., Mulot, B., Jones, K.
E., & Figueirido, B. (2023). The impact of the land-to-sea transition on
evolutionary integration and modularity of the pinniped backbone.
Communications Biology, 6, 1141.
