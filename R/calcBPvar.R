#' Calculate variability of breakpoints
#'
#' `calcBPvar()` computes an estimate of the variability of the breakpoints for a given number of regions. This involves computing the weighted mean and standard deviation of each breakpoint using Akaike weights.
#'
#' @param regions_results a `regions_results` object; the output of a call to [calcregions()] or [addregions()].
#' @param noregions the number of regions for which the weighted mean and standard deviation are to be computed.
#' @param pct the proportion of best model to keep from the original total number of possible models
#' @param criterion string; the criterion used to compute the weights. Allowable options include `"aic"` and `"bic"`. Abbreviations allowed.
#'
#' @returns A `regions_BPvar` object, which has two components:
#' * `WeightedBP` is a matrix containing the weighted mean and standard deviation of each breakpoint
#' * `BestModels` is a data frame containing the models used to compute the weighted breakpoint statistics and the weights each one is given.
#'
#' @seealso [calcregions()] for fitting segmented regression models to all combinations of breakpoints.
#'
#' @example man/examples/example-calcBPvar.R
#'

#' @export
calcBPvar <- function(regions_results, noregions, pct = .05, criterion = "aic") {
  arg::arg_is(regions_results, "regions_results")

  arg::arg_supplied(noregions)
  arg::arg_count(noregions)
  arg::arg_between(noregions, range(regions_results$stats$Nregions))

  arg::arg_number(pct)
  arg::arg_between(pct, c(0, 1), inclusive = c(FALSE, TRUE))

  criterion <- arg::match_arg(criterion, c("aic", "bic"))

  nmodel <- regions_results$stats$Nmodel_possible[regions_results$stats$Nregions == noregions]

  nvert <- length(attr(regions_results, "pos"))
  cont <- attr(regions_results, "cont")

  # Define number of models to keep to correspond to given pct of models wants to keep
  ntop <- ceiling(nmodel * pct)

  # Get number of PCs on which analysis was run
  noPC <- sum(startsWith(colnames(regions_results$results), "RSS."))

  # Extract models corresponding to the given number of regions
  dat <- regions_results$results[regions_results$results$regions == noregions, , drop = FALSE]

  # Order by increasing sumRSS (from best to worst model)
  dat <- dat[order(dat$sumRSS), , drop = FALSE]

  if (nrow(dat) >= ntop) {
    dat <- dat[seq_len(ntop), , drop = FALSE]
  }
  else {
    pct_calc <- round(nrow(dat) / nmodel * 100, 2L)
    arg::wrn("number of models provided lower than percentage requested. Weighted means and SD calculated on {pct_calc}% of total number of models")
  }

  # Calculate probability (with AICc/BIC) and weight of each model:
  IC <- vapply(dat$sumRSS, function(r) {
    ic <- .AICcalc(r, noPC = noPC, nvert = nvert, noregions = noregions, cont = cont)
    ic[switch(criterion, "aic" = "AICc", "bic" = "BIC")]
  }, numeric(1L))

  deltaIC <- IC - min(IC)	#Calculate AIC difference
  model_lik <- exp(-deltaIC / 2)	#Likelihood of the model

  IC_weight <- model_lik / sum(model_lik)	#Akaike Weights

  # Calculate weighted mean and weighted SD of each BP position using
  # Akaike weights: (formulae from Symonds & Moussalli, Behav Ecol Sociobiol 2011)
  bps <- as.matrix(dat[paste0("breakpoint", seq_len(noregions - 1L))])

  # Weighted mean is the sum of bps values multiplied by the corresponding Akaike weight
  wMean <- colSums(bps * IC_weight)
  wSD <- vapply(seq_len(ncol(bps)), function(i) {
    sqrt(sum(IC_weight * (bps[, i] - wMean[i])^2))
  }, numeric(1L))

  WeightedBP <- rbind(wMean, wSD)

  wt_dat <- data.frame(IC_weight, CumWeight = cumsum(IC_weight))
  names(wt_dat)[1L] <- switch(criterion, "aic" = "AICweight", "bic" = "BICweight")

  out <- list(WeightedBP = WeightedBP, BestModels = cbind(dat, wt_dat))

  attr(out, "pct") <- nrow(dat) / nmodel

  class(out) <- "regions_BPvar"

  out
  # WeightedBP = matrix with row1 = weighted means and row2 = weighted SD
  # BestModels = top pct % models ordered from best to worst with corresponding AIC weight and cumulative weight
}

#' @exportS3Method print regions_BPvar
print.regions_BPvar <- function(x, digits = 3, ...) {
  x0 <- x
  colnames(x$WeightedBP) <- sub("breakpoint", "BP ", colnames(x$WeightedBP), fixed = TRUE)
  x$WeightedBP <- round(x$WeightedBP, digits)

  print(x$WeightedBP)
  cat("- Computed using top ", round(100 * attr(x, "pct"), 2), "% of models\n", sep = "")

  invisible(x0)
}
