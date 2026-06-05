#' Calculate results of a single segmented regression model
#'
#' `calcmodel()` fits a multivariate segmented regression model using the supplied PCOs and breakpoints.
#'
#' @param x a `regions_pco` object; the output of a call to `svdPCO()`.
#' @param scores `numeric`; the indices of the PCO scores to use as the outcomes in fitting the model (e.g., `1:4` to use the first four scores).
#' @param bps `numeric`; the indices of the breakpoints to use in fitting the model. To request a model with no breakpoints, set to `NA`.
#' @param cont `logical`; whether to fit a model that is continuous (`TRUE`) or discontinuous (`FALSE`) at the breakpoints. Default is `TRUE`. Ignored when `bps` is `NA`.
#'
#' @returns A `regions_results_single` object, which contains the results of the model (breakpoints and RSS of each PCO and overall) and model support statistics.
#'
#' @seealso [calcregions()] and [addregions()] for computing all possible models instead of just a single one; [plotsegreg()], for which the `plot` method is an alias, for plotting the fitted regression lines; [modelsupport()] for interpreting the model support statistics.
#'
#' @example man/examples/example-calcmodel.R

#' @export
calcmodel <- function(x, scores, bps, cont = TRUE) {
  arg::arg_is(x, c("regions_pco", "regions_sim"))

  if (inherits(x, "regions_pco")) {
    arg::arg_supplied(scores)
    arg::arg_whole_numeric(scores)
    arg::arg_between(scores, c(1L, ncol(x[["scores"]])))

    Xvar <- .get_pos(x)
    Yvar <- x[["scores"]][, scores, drop = FALSE]
    eligible_vertebrae <- .get_eligible_vertebrae(x)
  }
  else {
    if (missing(scores)) {
      scores <- seq_len(ncol(x[["Yvar"]]))
    }
    Xvar <- x[["Xvar"]]
    Yvar <- x[["Yvar"]][, scores, drop = FALSE]
    eligible_vertebrae <- sort(unique(Xvar))
  }

  arg::arg_supplied(bps)

  if (is.atomic(bps) && all(is.na(bps))) {
    bps <- NA_real_
    cont <- TRUE
  }
  else {
    arg::arg_numeric(bps)
    arg::arg_between(bps, range(eligible_vertebrae))
    arg::arg_flag(cont)
  }

  # Re-order data
  o <- order(Xvar)
  Xvar <- Xvar[o]
  Yvar <- Yvar[o, , drop = FALSE]

  BPs <- sort(.drop_na(bps))
  nbp <- length(BPs)

  names(BPs) <- paste0("breakpoint", seq_len(nbp))

  noPC <- ncol(Yvar)

  colhead <- c("regions", paste0("breakpoint", seq_len(max(1L, nbp))),
               "sumRSS", paste("RSS", 1:noPC, sep = "."))

  #Calculate weights to ensure each vertebra counts equally
  vert_tab <- tabulate(Xvar)
  w <- 1 / vert_tab[Xvar]

  #Fit the model
  x <- .design_matrix(Xvar, if (!anyNA(BPs)) BPs, cont)

  lines <- .fast_lm(x = x, y = Yvar, w = w)

  RSS <- sum(w * lines$residuals^2)

  rsq <- {
    if (noPC > 1) colSums(w * lines$residuals^2)
    else RSS
  }

  res <- as.data.frame(matrix(c(nbp + 1, BPs, RSS, rsq), nrow = 1L,
                              dimnames = list(NULL, colhead)))

  supp <- .AICcalc(RSS, noPC = noPC, nvert = length(Xvar),
                  noregions = nbp + 1, cont = cont)

  out <- list(results = res,
              support = supp)

  attr(out, "scores") <- Yvar
  attr(out, "cont") <- cont
  attr(out, "nvert") <- length(Xvar)
  attr(out, "pos") <- Xvar

  class(out) <- "regions_results_single"

  out
}

#' @exportS3Method print regions_results_single
print.regions_results_single <- function(x, digits = 3, ...) {
  x0 <- x
  colnames(x[["results"]])[1L] <- "Regions"
  colnames(x[["results"]]) <- sub("breakpoint", "BP ", colnames(x[["results"]]),
                                  fixed = TRUE)
  print(round(x[["results"]], digits), row.names = FALSE)
  cat("\n- Support:\n")
  print(x[["support"]])

  invisible(x0)
}

#' @exportS3Method plot regions_results_single
plot.regions_results_single <- function(x, scores, ...) {
  arg::arg_supplied(scores)
  plotsegreg.regions_results_single(x, scores = scores, ...)
}
