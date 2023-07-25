#' Subsample a dataset
#'
#' `subsample()` creates a smaller version of the original dataset by sampling its rows. Because PCOs should be computed on the full dataset and most other functions take in `regions_pco` objects, `subsample()` requires a `regions_pco` object as its input.
#'
#' @param pco a `regions_pco` object; the output of a call to [svdPCO()].
#' @param sample `numeric`; either the number or proportion of vertebrae to remain the sampled dataset. If `NULL`, the original dataset is returned.
#' @param type string; the type of subsampling to do, either `"seq"` for sampling in sequence or `"random"` for random sampling. Default is `"seq"`. Abbreviations allowed.
#'
#' @returns A `regions_pco` object, a subset of the original supplied to `pco`. The original dataset is stored as an attribute, which itself contains the subsampling indices.
#'
#' @seealso [svdPCO()], [process_measurements()], [plotvertmap()] to visualize the vertebral map after subsampling.
#'
#' @example man/examples/example-subsample.R
#'

#' @export
subsample <- function(pco, sample = NULL, type = "seq") {

  chk::chk_is(pco, "regions_pco")

  if (nrow(attr(pco, "data")) != length(attr(attr(pco, "data"), "subset"))) {
    chk::err("`subsample()` cannot be used on a `regions_pco` object after using `subsample()` or `subset(., drop = FALSE)` on it.")
  }

  dat <- attr(pco, "data")
  pos <- dat[[attr(dat, "pos_ind")]]

  chk::chk_number(sample)
  chk::chk_gt(sample, 0)
  chk::chk_lte(sample, length(pos))

  chk::chk_string(type)
  type <- tolower(type)
  type <- .match_arg(type, c("seq", "random"))

  if (sample <= 1) {
    sample <- round(sample * length(pos))
  }

  sampled_pos_ind <- switch(type,
                            "seq" = round(seq(1, length(pos), length.out = sample)),
                            "random" = sample(seq_along(pos), sample))

  pco$scores <- pco$scores[sampled_pos_ind,, drop = FALSE]
  pco$eigen.val <- pco$eigen.val[sampled_pos_ind]

  attr(attr(pco, "data"), "subset") <- sampled_pos_ind

  pco
}
