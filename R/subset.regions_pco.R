#' Subset a `regions_pco` object
#'
#' `subset()` allows one to subset a `regions_pco` object, the output of a call to [svdPCO()]. This function is useful only to include other observations when computing the PCOs but to focus the analysis on a core set of observations.
#'
#' @param x a `regions_pco` object, the output of a call to [svdPCO()].
#' @param subset a `logical` vector with length equal to the number of rows of the original dataset passed to `process_measurements()` before `svdPCO()`. The vector is evaluated in the environment of the original dataset.
#' @param drop `logical`. When `TRUE`  (the default), the observations not in the subset are removed from the dataset. When `FALSE`, all observations remain, but the indices of the subsetted observations are stored. This only affects how `plotvertmap()` displays observations when subsetting changes the range of the vertebra indices. See Examples.
#' @param dots ignored.
#'
#' @returns A `regions_pco` object subsetted according to the argument passed to `subset`.
#'
#' @seealso
#' [subsample()] for randomly subsetting a `regions_pco` object.
#'
#' @example man/examples/example-subset.R
#'
#' @exportS3Method subset regions_pco
subset.regions_pco <- function(x, subset, drop = TRUE, ...) {
  chkDots(...)
  if (missing(subset)) return(x)

  dat <- attr(x, "data")

  e <- substitute(subset)
  r <- eval(e, dat, parent.frame())
  chk::chk_logical(r, "`subset`")

  r <- r & !is.na(r)

  if (all(r)) return(x)

  if (nrow(dat) != length(attr(dat, "subset"))) {
    chk::err("`subset()` cannot be used on a `regions_pco` object after using `subsample()` on it.")
  }

  chk::chk_flag(drop)

  x$scores <- x$scores[r, , drop = FALSE]
  x$eigen.val <- x$eigen.val[r]

  if (drop) {
    attr(x, "data") <- dat[r, , drop = FALSE]
    attr(attr(x, "data"), "pos_ind") <- attr(dat, "pos_ind")
    attr(attr(x, "data"), "subset") <- seq_len(sum(r))
  }
  else {
    attr(attr(x, "data"), "subset") <- which(r)
  }

  x
}