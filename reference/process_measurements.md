# Process vertebra measurements

`process_measurements()` initializes the analysis workflow by processing
a dataset of vertebra measurements into an object usable by
MorphoRegions. Such processing includes identifying the vertebra indices
and the measurements and filling in missing values.

## Usage

``` r
process_measurements(data, pos = 1L, measurements, fillNA = TRUE)
```

## Arguments

- data:

  a data.frame containing a column of vertebra indices and measurements
  for each vertebra, or a list thereof for multiple specimens.

- pos:

  the name or index of the variable in `data` containing the vertebra
  indices. Default is to use the first column.

- measurements:

  the names or indices of the variables in `data` containing the
  relevant vertebra measurements. If unspecified, will use all variables
  other than that specified in `pos`.

- fillNA:

  `logical`; whether to fill in missing values using a simple linear
  imputation. Default is `TRUE`. See Details.

## Value

A `regions_data` object, which is a list of data.frames (one for each
specimen) with attributes containing metadata.

## Details

Any rows with missing values for all measurements will be removed. When
missing values in non-removed rows are present and `fillNA` is set to
`TRUE`, `process_measurements()` fills them in if the sequence of
missing values is no greater than 2 in length. For numeric variables, it
uses a linear interpolation, and for categorical variables, it fills in
the missing values with the surrounding non-missing values if they are
identical and leaves them missing otherwise. Otherwise, missing values
are left as they are.

When a list of data frames is supplied to `data`, only the variables
named in `measurements` that are common across datasets will be stored
as measurement variables.

## See also

[`svdPCO()`](https://aagillet.github.io/MorphoRegions/reference/svdPCO.md)
for computing principal coordinate axes from processed vertebra data.

## Examples

``` r
# Process dataset; vertebra index in "Vertebra" column
data("alligator")

alligator_data <- process_measurements(alligator,
                                       pos = "Vertebra")

# Process multiple datasets; vertebra index in first column
data("porpoise")

porpoise_data <- process_measurements(list(porpoise1,
                                           porpoise2,
                                           porpoise3),
                                      pos = 1)
```
