#' Calculate results of a single segmented regression model
#'
#' `calcmodel()` fits a multivariate segmented regression model using the supplied PCOs and breakpoints.
#'
#' @param x a `regions_pco` object; the output of a call to `svdPCO()`.
#' @param scores `numeric`; the indices of the PCO scores to use as the outcomes in fitting the model (e.g., `1:4` to use the first four scores).
#' @param bps `numeric`; the indices of the breakpoints to use in fitting the model. To request a model with no breakpoints, set to `NA`.
#' @param cont `logical`; whether to fit a model that is continuous (`TRUE`) or discontinuous (`FALSE`) at the breakpoints. Default is `TRUE`. Ignored when `bps` is `NA`.
#'
#' @returns A `regions_results_single` objects, which contains the results of the model (breakpoints and RSS of each PCO and overall) and model support statistics.
#'
#' @seealso [calcregions()] and [addregions()] for computing all possible models instead of just a single one; [plotsegreg()], for which the `plot` method is an alias, for plotting the fitted regression lines; [modelsupport()] for interpreting the model support statistics.
#'
#' @example man/examples/example-calcmodel.R

#' @export
calcmodel <- function(x, scores, bps, cont = TRUE) {
  if (inherits(x, "regions_sim")) {
    if (missing(scores)) {
      scores <- seq_len(ncol(x[["Yvar"]]))
    }
    Xvar <- x[["Xvar"]]
    Yvar <- x[["Yvar"]][, scores, drop = FALSE]
  }
  else if (inherits(x, "regions_pco")) {
    chk::chk_not_missing(scores, "`scores`")
    chk::chk_whole_numeric(scores)
    chk::chk_range(scores, c(1, ncol(x[["scores"]])))

    subset <- attr(attr(x, "data"), "subset")
    Xvar <- attr(attr(x, "data"), "pos")[subset]
    Yvar <- x[["scores"]][subset, scores, drop = FALSE]
  }
  else {
    chk::err("`x` must be a `regions_pco` or `regions_sim` object")
  }

  chk::chk_not_missing(bps, "`bps`")


  if (chk::vld_atomic(bps) && all(is.na(bps))) {
    bps <- NA_real_
    cont <- TRUE
  }
  else {
    chk::chk_numeric(bps)
    chk::chk_range(bps, range(Xvar))
    chk::chk_flag(cont)
  }

  BPs <- sort(.drop_na(bps))
  nbp <- length(BPs)

  names(BPs) <- paste0("breakpoint", seq_len(nbp))

  noPC <- ncol(Yvar)

  colhead <- c("regions", paste0("breakpoint", seq_len(max(1, nbp))),
               "sumRSS", paste("RSS", 1:noPC, sep = "."))

  #Fit the model
  x <- .design_matrix(Xvar, BPs, cont)

  lines <- .lm.fit(x = x, y = Yvar)

  RSS <- sum(lines$residuals^2)
  if (noPC > 1) {
    rsq <- colSums(lines$residuals^2)
  }
  else {
    rsq <- RSS
  }

  res <- as.data.frame(matrix(c(nbp + 1, BPs, RSS, rsq), nrow = 1,
                              dimnames = list(NULL, colhead)))

  supp <- AICcalc(RSS, noPC = noPC, nvert = length(Xvar),
                  noregions = nbp + 1, cont = cont)

  out <- list(results = res,
              support = supp)

  attr(out, "scores") <- Yvar
  attr(out, "cont") <- cont
  attr(out, "nvert") <- length(Xvar)
  attr(out, "pos") <- Xvar

  class(out) <- "regions_results_single"

  return(out)
}

#' @exportS3Method print regions_results_single
print.regions_results_single <- function(x, digits = 3, ...) {
  x0 <- x
  colnames(x[["results"]])[1] <- "Regions"
  colnames(x[["results"]]) <- sub("breakpoint", "BP ", colnames(x[["results"]]),
                                  fixed = TRUE)
  print(round(x[["results"]], digits), row.names = FALSE)
  cat("\n- Support:\n")
  print(x[["support"]])

  invisible(x0)
}

#' @exportS3Method plot regions_results_single
plot.regions_results_single <- function(x, scores, ...) {
  chk::chk_not_missing(scores, "`scores`")
  plotsegreg.regions_results_single(x, scores = scores, ...)
}
