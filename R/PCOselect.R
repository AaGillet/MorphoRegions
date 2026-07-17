#' Select PCO scores
#'
#' `PCOselect()` provides several methods to select the number of principal coordinates (PCOs) analysis scores to be used in subsequent analyses.
#'
#' @param pco a `regions_pco` object; the output of a call to [svdPCO()].
#' @param method string; the method used to select the number of PCOs. Allowable options include `"manual"`, `"boot"`, `"variance"`, and `"max"`. Default is `"manual"`. Abbreviations allowed. See Details.
#' @param scores when `method = "manual"`, the number of PCO scores to use.
#' @param nreps when `method = "boot"`, the number of bootstrap replications to use.
#' @param cutoff when `method = "variance"`, the cutoff for the variance explained by each PCO score.
#' @param results when `method = "max"`, a `regions_results` object, the output of a call to [calcregions()].
#' @param criterion when `method = "max"`, which criterion should be used to select the number of scores. Allowable options include `"aic"` and `"bic"`. Abbreviations allowed.
#' @param verbose when `method = "boot"`, whether to display a progress bar. Default is `TRUE` when running interactively and `FALSE` otherwise.
#' @param x for `plot.regions_pco_select()`, a `regions_pco_select` object, the output of a call to `PCOselect()` with `method = "boot"` or `"max"`.
#' @param object a `regions_pco_select` object, the output of a call to `PCOselect()` with `method = "max"`.
#' @param \dots ignored.
#'
#' @return For `PCOselect()`, a `regions_pco_select` object, which is a numeric vector containing the indices of the chosen PCs, with attributes containing information about the PC scores chosen by the specified method. When `method = "boot"`, the bootstrap results are stored in the `"boot"` attribute. When `method = "max"`, the `regions_results` object passed to `regions` and other information about the quality of fit for each number of PCs are stored in the `"pcomax"` attribute.
#'
#' The `plot()` methods each return a `ggplot` object that can manipulated using \pkg{ggplot2} syntax. The `summary()` method returns a data.frame of results.
#'
#' @details
#' Each method provides an alternate way to select the number of scores. These are described below.
#'
#' ## `method = "manual"`:
#'
#' This simply returns the number supplied to `scores` after running some checks to ensure it is valid.
#'
#' ## `method = "boot"`
#'
#' Bootstrapping works by comparing the eigenvalue distributions of PCOs to those with randomized data in order to extract PCO axes with significant signal, which are defined as those with eigenvalues greater than those from randomized data. The returned PCO cutoff is the largest PCO axis whose eigenvalues fall below the mean eigenvalue for that axis from the randomized data. Data are randomly sampled by row. Bootstrapping is sensitive to unequal variances of columns, so `scale = TRUE` should be set in the call to [svdPCO()], which is the default; the data are scaled in the same way prior to bootstrapping. The `plot()` method displays the eigenvalues of the true PCOs and boxplots summarizing the distribution of the bootstrapped eigenvalues for each PCO. The `boot` method is unavailable for user-supplied PC scores (`regions_pco` object obtained from `process_PC` or `process_gmPC` functions).
#'
#' ## `method = "variance"`
#'
#' This method works by computing the ratio of each eigenvalue to the sum of the eigenvalues (i.e., to compute the proportion of variance explained by each PC score) and select the number of scores with ratios greater than the cutoff value supplied to `cutoff`.
#'
#' ## `method = "max"`
#'
#' This method works by selecting the smallest number of PCs that gives a region score within .001 of the maximum possible region score for the segmented models fit in the object supplied to `results`. Which criterion is maximized (AIC or BIC) is determined by the value supplied to `criterion`. The `summary()` method displays the region score (estimated number of regions) for each PC (`RSind`) and for PCs cumulatively (`RScum`) selected using the AICc or BIC as well as the cumulative proportion of variance explained by the PCs. The `plot()` method displays this information graphically, with the left y-axis displaying the region score for the PCs individually (pale blue triangles) and cumulatively (orange circles) using each of the two criteria, and the right y-axis displaying the cumulative percentage of variance explained by the PCs.
#'
#' @example man/examples/example-PCOselect.R
#'

#' @export
PCOselect <- function(pco, method = "manual", scores = NULL, cutoff = .05, nreps = 500,
                      results = NULL, criterion = "aic", verbose = interactive()) {
  arg::arg_is(pco, "regions_pco")

  method <- arg::match_arg(method, c("manual", "boot", "variance", "max"))

  boot <- pcomax <- NULL
  info <- list(method = method)
  if (method == "manual") {
    arg::arg_non_null(scores)
    arg::arg_count(scores,
                   .msg = "{.arg scores} must be the number of PC scores to select")
    arg::arg_between(scores, c(1, ncol(pco$scores)),
                     .msg = sprintf("{.arg scores} must be between 1 and the total number of principal coordinates (i.e., %s)",
                                    ncol(pco$scores)))
  }
  else if (method == "boot") {
    arg::arg_count(nreps)
    arg::arg_gte(nreps, 1)
    arg::arg_flag(verbose)

    boot <- .PCOboot(pco, nreps, verbose)
    scores <- boot$sigpco

    info$nreps <- nreps
  }
  else if (method == "variance") {
    arg::arg_number(cutoff)
    arg::arg_between(cutoff, c(0, 1), inclusive = FALSE)

    scores <- sum(pco$eigen.val / sum(pco$eigen.val) > cutoff)
    info$cutoff <- cutoff
  }
  else {
    arg::arg_non_null(results)
    arg::arg_is(results, "regions_results")

    criterion <- arg::match_arg(criterion, c("aic", "bic"))

    pcomax <- .PCOmax(results)

    scores <- pcomax[[switch(criterion,
                             aic = "pco.max.AICc",
                             bic = "pco.max.BIC")]]

    pcomax$pco <- pco
    pcomax$results <- results

    info$criterion <- criterion
  }

  out <- seq_len(scores)

  attr(out, "info") <- info
  attr(out, "boot") <- boot
  attr(out, "pcomax") <- pcomax

  class(out) <- c("regions_pco_select", class(out))

  out
}

#' @exportS3Method print regions_pco_select
print.regions_pco_select <- function(x, ...) {
  info <- attr(x, "info")
  cat("A `regions_pco_select` object\n")
  cat(sprintf("- PC scores selected: %s\n", toString(x[])))
  cat(sprintf("- Method: %s\n",
              switch(info$method,
                     "boot" = sprintf("boot (%s replications)", info$nreps),
                     "variance" = sprintf("variance (cutoff: %s)", info$cutoff),
                     "max" = sprintf("max (criterion: %s)", toupper(info$criterion)),
                     info$method)))
  invisible(x)
}

#' @exportS3Method plot regions_pco_select
#' @rdname PCOselect
plot.regions_pco_select <- function(x, ...) {
  if (!identical(attr(x, "info")$method, "boot") &&
      !identical(attr(x, "info")$method, "max")) {
    arg::err("{.fun plot} can only be used on {.cls regions_pco_select} objects when {.arg method} is {.val boot} or {.val max}")
  }

  if (identical(attr(x, "info")$method, "boot")) {
    boot <- attr(x, "boot")
    ind <- seq_along(boot$eigen.true)

    eigen.boot <- as.vector(boot$eigen.boot)
    ind.boot <- as.vector(row(boot$eigen.boot))

    p <- ggplot() +
      geom_point(aes(x = ind, y = boot$eigen.true), shape = "circle filled") +
      geom_line(aes(x = ind, y = boot$eigen.true)) +
      geom_boxplot(aes(y = eigen.boot, x = ind.boot, group = factor(ind.boot)),
                   outlier.shape = NA) +
      labs(x = "PC axis", y = "Eigenvalue", title = "Eigenvalue cutoff") +
      theme_bw()
  }
  else {
    s <- summary(x, plot = FALSE)

    p <- .plot_summary_regions_pco_select(s)
  }

  p
}

#' @exportS3Method summary regions_pco_select
#' @rdname PCOselect
summary.regions_pco_select <- function(object, ...) {
  if (!identical(attr(object, "info")$method, "max") ||
      is.null(attr(object, "pcomax"))) {
    arg::err("{.fun summary} can only be used on {.cls regions_pco_select} objects when {.arg method} is {.val max}")
  }

  results <- attr(object, "pcomax")$results
  eigenvals <- attr(object, "pcomax")$pco$eigen.val

  noregions <- max(results$stats$Nregions)

  nPCO <- sum(startsWith(colnames(results$results), "RSS."))

  # Calculate variance of each PCO:
  var.exp <- cumsum(eigenvals / sum(eigenvals))[seq_len(nPCO)]

  pco.no.test <- data.frame(PCO = seq_len(nPCO),
                            RSind.AICc = NA_real_,
                            RScum.AICc = NA_real_,
                            RSind.BIC = NA_real_,
                            RScum.BIC = NA_real_,
                            CumulVar = var.exp)

  for (a in seq_len(nPCO)) {
    #Run for individual PCs
    models.ind <- modelselect(results, scores = a)
    support.ind <- modelsupport(models.ind)

    #Run for cumulative PCs
    models.cum <- modelselect(results, scores = 1:a)
    support.cum <- modelsupport(models.cum)

    pco.no.test[["RSind.AICc"]][a] <- support.ind$Region_score
    pco.no.test[["RScum.AICc"]][a] <- support.cum$Region_score
    pco.no.test[["RSind.BIC"]][a] <- support.ind$Region_score_BIC
    pco.no.test[["RScum.BIC"]][a] <- support.cum$Region_score_BIC
  }

  attr(pco.no.test, "noregions") <- noregions
  class(pco.no.test) <- c("summary.regions_pco_select", class(pco.no.test))

  pco.no.test
}

#' @exportS3Method print summary.regions_pco_select
print.summary.regions_pco_select <- function(x, digits = 3L, ...) {
  arg::arg_whole_number(digits)

  x0 <- x
  for (i in seq_len(ncol(x))[-1]) {
    x[[i]] <- round(x[[i]], digits)
  }
  print.data.frame(x, row.names = FALSE, ...)
  invisible(x0)
}

.PCOboot <- function(pco, nreps = 500, verbose = FALSE) {

  if (verbose) {
    cat("Bootstrapping...\n")
  }
  else {
    opb <- pbapply::pboptions(type = "none")
    on.exit(pbapply::pboptions(opb))
  }

  #### Edited code to add support for GMM data: ##
  metric <- attr(pco, "metric")
  if (identical(metric, "custom")) {
    arg::err("method {.val boot} is not allowed for user-supplied PC scores")
  }

  #calculate 'true' eigenvalues as percentage variance
  eigen.true <- prop.table(pco$eigen.val)

  data <- .get_data_without_pos(pco)

  randdata <- data

  eigen.boot <- do.call("cbind", pbapply::pblapply(seq_len(nreps), function(i) {
    #Shuffle each row of the dataset
    for (j in seq_len(nrow(data))) {
      randdata[j, ] <- sample(data[j, ])
    }

    dist <- cluster::daisy(randdata, metric = metric, stand = attr(pco, "scale"))

    if (anyNA(dist)) {
      return(NULL)
    }

    #calculate bootstrapped eigenvalues
    pco <- .svdPCO_internal(dist, val.only = TRUE)

    #calculate as percentage variance
    prop.table(pco$eigen.val)
  }))

  eigen.mean <- rowMeans(eigen.boot)  #calculate mean and SD of bootstrapped values
  eigen.sd <- apply(eigen.boot, 1L, sd)
  diff <- eigen.true - eigen.mean  #figure out which PCOs have greater eigenvalues for the 'true' dataset
  diff[diff < 0] <- 0

  #split the dataset at the zeros, and calculate number of pcos in the first string
  sigpco <- length(split(diff, cumsum(diff == 0))[["0"]])

  list(eigen.true = eigen.true,
       eigen.mean = eigen.mean,
       eigen.sd = eigen.sd,
       sigpco = sigpco,
       eigen.boot = eigen.boot)
}

.PCOmax <- function(results, tol = .001) {

  regiondata <- results$results

  nvar <- sum(startsWith(colnames(regiondata), "RSS."))

  pco.no.test <- data.frame(PCO = seq_len(nvar),
                            RS_AICc = NA_real_,
                            RS_BIC = NA_real_)

  for (i in seq_len(nvar)) {
    #Run for cumulative PCs
    models.cum <- modelselect(results, scores = seq_len(i))
    support.cum <- modelsupport(models.cum)

    pco.no.test[["RS_AICc"]][i] <- support.cum$Region_score
    pco.no.test[["RS_BIC"]][i] <- support.cum$Region_score_BIC
  }

  pco.max.AICc <- which(.equiv(pco.no.test[["RS_AICc"]], max(pco.no.test[["RS_AICc"]]), tol = tol))[1L]
  pco.max.BIC <- which(.equiv(pco.no.test[["RS_BIC"]], max(pco.no.test[["RS_BIC"]]), tol = tol))[1L]

  list(pco.max.AICc = pco.max.AICc,
       pco.max.BIC = pco.max.BIC,
       pco.dist = pco.no.test)
}

.plot_summary_regions_pco_select <- function(x, ...) {
  arg::arg_is(x, "summary.regions_pco_select")

  pco.no.test.long <- reshape(x, direction = "long", idvar = "PCO", varying = list(c(2, 3), c(4, 5)),
                              timevar = "PCOtype", times = c("Single PCO", "Cumulated PCOs"),
                              v.names = c("RS.AICc", "RS.BIC"))
  pco.no.test.long <- reshape(pco.no.test.long, direction = "long", varying = 4:5,
                              timevar = "Testtype", times = c("AICc", "BIC"),
                              v.names = "value")

  noregions <- attr(x, "noregions")

  p <- ggplot(pco.no.test.long) +
    geom_point(data = pco.no.test.long,
               aes(x = .data$PCO, y = .data$value,
                   color = .data$PCOtype, shape = .data$PCOtype)) +
    scale_color_manual(values = c("#fc8d62", "#8da0cb")) +
    geom_line(data = x, aes(x = .data$PCO, y = (.data$CumulVar * (noregions - 1) + 1)),
              color = "darkgrey", linewidth = 1) +
    scale_y_continuous(name = "Region score",
                       sec.axis = sec_axis(~ (. - 1) / (noregions - 1), labels = scales::percent,
                                           # limits = c(0, 1),
                                           name = "Cumulated variance explained")) +
    scale_x_continuous(breaks = scales::breaks_extended(Q = c(0:5))) +
    facet_wrap(~Testtype) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom") +
    labs(color = "Region score", shape = "Region score")

  p
}
