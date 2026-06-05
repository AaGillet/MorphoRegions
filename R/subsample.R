#' Subsample a dataset
#'
#' `subsample()` creates a smaller version of the original dataset by sampling its rows. Because PCOs should be computed on the full dataset and most other functions take in `regions_pco` objects, `subsample()` requires a `regions_pco` object as its input.
#'
#' @param pco a `regions_pco` object; the output of a call to [svdPCO()].
#' @param sample `numeric`; either the number or proportion of vertebrae to remain the sampled dataset. If `NULL`, the original dataset is returned.
#' @param type string; the type of subsampling to do, either `"seq"` for sampling in sequence or `"random"` for random sampling. Default is `"seq"`. Abbreviations allowed.
#'
#' @return
#' A `regions_pco` object, a subset of the original supplied to `pco`. The original dataset is stored as an attribute, which itself contains the subsampling indices.
#'
#' @seealso [svdPCO()], [process_measurements()], [plotvertmap()] to visualize the vertebral map after subsampling.
#'
#' @example man/examples/example-subsample.R
#'

#' @export
subsample <- function(pco, sample = NULL, type = "seq") {

  arg::arg_supplied(pco)
  arg::arg_is(pco, "regions_pco")

  pos <- .get_pos(pco)

  if (!identical(pos, .get_pos(pco, subset = FALSE))) {
    arg::err("{.fun subsample} cannot be used on a {.cls regions_pco} object after using {.fun subsample} or {.code subset(., drop = FALSE)} on it")
  }

  eligible_vertebrae <- .get_eligible_vertebrae(pco)

  arg::arg_number(sample)
  arg::arg_gt(sample, 0)
  arg::arg_lte(sample, length(eligible_vertebrae))

  type <- arg::match_arg(type, c("seq", "random"))

  if (sample <= 1) {
    sample <- round(sample * length(eligible_vertebrae))
  }

  sampled_pos_ind <- switch(type,
                            "seq" = round(seq(1, length(eligible_vertebrae), length.out = sample)),
                            "random" = sample(seq_along(eligible_vertebrae), sample))

  pco$scores <- pco$scores[pos %in% eligible_vertebrae[sampled_pos_ind],, drop = FALSE]
  pco$eigen.val <- pco$eigen.val[pos %in% eligible_vertebrae[sampled_pos_ind]]

  attr(attr(pco, "data"), "eligible_vertebrae") <- eligible_vertebrae[sampled_pos_ind]

  pco
}
