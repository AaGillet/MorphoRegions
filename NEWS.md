# *MorphoRegions* News and Updates

## Version 0.2.0.9000 (development)

### New features
- Added support for geometric morphometric (GMM) data
- Allowed direct input of PCA scores computed outside MorphoRegions

### Changes
- Added new function `process_gmPC()` to allow direct input of PC scores from 2D or 3D geometric morphometric data
- Added new function `process_PC()` to allow direct input of PC scores from traditional morphometric data
- Added new data to demonstrate new GMM features
- Added a `NEWS.md` file to track changes to the package
- The `verbose` argument is now set to `FALSE` when not in an interactive session

### Notes
- This is a development version; functionality may change before CRAN release