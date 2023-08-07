#' Fit segmented regression models for all combinations of breakpoints
#'
#' `calcregions()` enumerates all possible combinations of breakpoints to fit multivariate segmented regression models. `addregions()` adds models with additional numbers of regions to the resulting output object.
#'
#' @param pco a `regions_pco` object; the output of a call to [svdPCO()].
#' @param scores `numeric`; the indices of the PCO scores to use as the outcomes in fitting the models (e.g., `1:4` to use the first four scores). Can also be the ouput of a call to [PCOselect()].
#' @param noregions `numeric`; for `calcregions()`, the maximum number of regions for which models are fit (e.g, 4 to request models with 1 to 4 regions); for `addregions()`, a vector containing the number of regions to add (e.g., 5:6 to request models with 5 and 6 regions).
#' @param minvert `numeric`; the minimum number of vertebrae allowed in each region. Default is 3.
#' @param cont `logical`; whether to fit models that are continuous (`TRUE`) or discontinuous (`FALSE`) at the breakpoints. Default is `TRUE`.
#' @param exhaus `logical`; whether to fit all possible models (`TRUE`) or use heuristics to reduce the number of models fit (`FALSE`). Default is `TRUE`. See Details. Setting to `FALSE` can reduce the size of the resulting object.
#' @param omitbp an optional vector of vertebrae to be omitted from the list of possible breakpoints, e.g., if it is known that two adjacent vertebrae belong to the same region.
#' @param cl a cluster object created by [parallel::makeCluster()], an integer to indicate number of child-processes (integer values are ignored on Windows) for parallel evaluations, or `"future"` to use a future backend. `NULL` (the default) refers to sequential evaluation (no parallelization). See [pbapply::pbapply()] for details.
#' @param verbose `logical`; whether to print information about the fitting process, including a progress bar. Default is `TRUE`.
#' @param regions_results,object a `regions_results` object; the output of a call to `calcregions()` or `addregions()`.
#' @param \dots ignored.
#'
#' @returns A `regions_results` object with the following components:
#' * `results` - the results of the fitting process for each combination of breakpoints
#' * `stats` - statistics summarizing the fitting process. Use `summary()` to view this information in a clean format.
#'
#' @details `calcregions()` enumerates all possible combinations of breakpoints that satisfy the constraint imposed by `minvert` (i.e., that breakpoints need to be at least `minvert` vertebrae apart) and fits the segmented regression models implied by each combination. These are multivariate regression models with the PCO scores specified by `scores` as the outcome. When `cont = TRUE`, these regression models are continuous; i.e., the regression lines for each region connect at the breakpoints. Otherwise, the models are discontinuous so that each region has its own intercept and slope. The models are fit using [.lm.fit()], which efficiently implements ordinary least squares regression.
#'
#' When `exhaus = FALSE`, heuristics are used to reduce the number of models to fit, which can be useful for keeping the size of the resulting object down by avoiding fitting models corresponding to breakpoint combinations that yield a poor fit to the data. Only breakpoint combinations that correspond to the breakpoints of the best fitting model with a smaller number of regions +/- 3 vertebrae are used, and only models that have an RSS smaller than half a standard deviation more the smallest RSS are kept.
#'
#' `addregions()` should be used on an existing `regions_results` object to add models with more regions. Internally, it works just the same as `calcregions()`.
#'
#' @seealso [calcmodel()] to fit a segmented regression model for a single set of breakpoints; [modelselect()] to select the best model for each number of regions based on RSS; [modelsupport()] to compute statistics the describe the support of the best models; [calcBPvar()] to compute the variability in the optimal breakpoints.
#'
#' @example man/examples/example-calcregions.R

#' @export
calcregions <- function(pco, scores, noregions, minvert = 3, cont = TRUE,
                        exhaus = TRUE, omitbp = NULL, cl = NULL, verbose = TRUE) {
  # Argument checks
  if (inherits(pco, "regions_pco")) {
    chk::chk_not_missing(scores, "`scores`")
    chk::chk_whole_numeric(scores)
    chk::chk_range(scores, c(1, ncol(pco[["scores"]])))

    # subset <- attr(attr(pco, "data"), "subset")
    # Xvar <- attr(attr(pco, "data"), "pos")[subset]
    Xvar <- .get_pos(attr(pco, "data"))
    Yvar <- pco[["scores"]][, scores, drop = FALSE]
  }
  else if (inherits(pco, "regions_sim")) {
    if (missing(scores)) {
      scores <- seq_len(ncol(pco[["Yvar"]]))
    }
    else {
      chk::chk_whole_numeric(scores)
      chk::chk_range(scores, c(1, ncol(pco[["Yvar"]])))
    }

    Xvar <- pco[["Xvar"]]
    Yvar <- pco[["Yvar"]][, scores, drop = FALSE]
  }
  else {
    chk::err("`pco` must be a `regions_pco` or `regions_sim` object")
  }

  chk::chk_not_missing(noregions, "`noregions`")
  chk::chk_count(noregions)
  chk::chk_gte(noregions, 1)

  chk::chk_count(minvert)
  chk::chk_gte(minvert, 2)
  chk::chk_flag(cont)
  chk::chk_flag(exhaus)
  chk::chk_flag(verbose)

  if (!is.null(omitbp))
    chk::chk_numeric(omitbp)

  # Re-order data
  o <- order(Xvar)
  Xvar <- Xvar[o]
  Yvar <- Yvar[o, , drop = FALSE]

  ncombs <- vapply(seq_len(noregions), function(i) {
    .ncombos(Xvar, minvert, i - 1)
  }, numeric(1L))

  if (any(ncombs == 0)) {
    max.regions <- max(which(ncombs > 0))
    chk::wrn(sprintf("models for %s vertebrae cannot be fit with more than %s region%%s and `minvert` of %s",
                     length(Xvar), max.regions, minvert), n = max.regions)
    noregions <- max.regions
  }

  opb <- pbapply::pboptions(type = if (verbose) "timer" else "none")
  on.exit(pbapply::pboptions(opb))

  names(Xvar) <- seq_along(Xvar) #pos

  noPC <- ncol(Yvar)

  colhead <- c("regions", "breakpoint1", "sumRSS", paste("RSS", 1:noPC, sep = "."))

  # 1 Region (0 breakpoint - bp)

  if (verbose) {
    cat("Fitting model with 1 region (0 breakpoints)... ")
  }

  nbp <- 0		# Number of breakpoints for 1 region
  lines <- .lm.fit(x = cbind(1, Xvar), y = Yvar)

  RSS <- sum(lines$residuals^2)
  if (noPC > 1) {
    rsq <- colSums(lines$residuals^2)
  }
  else {
    rsq <- RSS
  }

  stats <- data.frame(Nregions = nbp + 1L,
                      Nmodel_possible = 1L,
                      Nmodel_tested = 1L,
                      Nmodel_saved = 1L,
                      Comp_method = "Exhaustive",
                      Saving_method = "All")

  res <- as.data.frame(matrix(c(nbp + 1, NA_real_, RSS, rsq), nrow = 1,
                              dimnames = list(NULL, colhead)))
  best_bps <- matrix(NA_character_, nrow = 1, ncol = 1,
                     dimnames = list(NULL, "Best_BPs1"))


  if (verbose) {
    cat("Done.\n")
  }

  rownames(res) <- NULL

  stats <- cbind(stats, as.data.frame(best_bps))

  out <- list(results = res, stats = stats)

  attr(out, "scores") <- Yvar
  attr(out, "cont") <- cont
  attr(out, "minvert") <- minvert
  attr(out, "omitbp") <- omitbp
  attr(out, "pos") <- Xvar

  class(out) <- "regions_results"

  for (noregion in seq_len(noregions)[-1]) {
    out <- .addregions_internal(out, noregion, exhaus = exhaus, cl = cl, verbose = verbose)
  }

  out
}

#' @export
#' @rdname calcregions
addregions <- function(regions_results, noregions, exhaus = FALSE, cl = NULL, verbose = TRUE) {
  # Argument checks
  chk::chk_is(regions_results, "regions_results")

  chk::chk_not_missing(noregions, "`noregions`")
  chk::chk_whole_numeric(noregions)
  chk::chk_gte(noregions, 1)
  noregions <- sort(noregions)

  chk::chk_flag(exhaus)
  chk::chk_flag(verbose)

  opb <- pbapply::pboptions(type = if (verbose) "timer" else "none")
  on.exit(pbapply::pboptions(opb))

  Xvar <- attr(regions_results, "pos")
  minvert <- attr(regions_results, "minvert")
  ncombs <- vapply(seq_len(max(noregions)), function(i) {
    .ncombos(Xvar, minvert, i - 1)
  }, numeric(1L))

  if (any(ncombs == 0)) {
    max.regions <- max(which(ncombs > 0))
    chk::wrn(sprintf("models for %s vertebrae cannot be fit with more than %s region%%s and `minvert` of %s",
                     length(Xvar), max.regions, minvert), n = max.regions)
    noregions <- noregions[noregions <= max.regions]
  }

  for (n in noregions) {
    regions_results <- .addregions_internal(regions_results, n, exhaus, cl, verbose)
  }

  regions_results
}

#' @exportS3Method print regions_results
print.regions_results <- function(x, ...) {
  cat("A `regions_results` object\n")
  cat(" - number of PCOs used:", ncol(attr(x, "scores")), "\n")
  cat(" - number of regions:", paste(sort(unique(x$stats[["Nregions"]])), collapse = ", "), "\n")
  cat(" - model type:", if (attr(x, "cont")) "continuous" else "discontinuous", "\n")
  cat(" - min vertebrae per region:", attr(x, "minvert"), "\n")
  if (!is.null(attr(x, "omitbp"))) {
    cat(" - omitted breakpoints:", paste(attr(x, "omitbp"), collapse = ", "), "\n")
  }
  cat(" - total models saved:", nrow(x$results), "\n")
  cat("Use `summary()` to examine summaries of the fitting process.\n")

  invisible(x)
}

#' @exportS3Method summary regions_results
#' @rdname calcregions
summary.regions_results <- function(object, ...) {
  stats <- object$stats[1:6]
  names(stats) <- c("Regions", "Possible", "Tested", "Saved",
                    "Comp. method", "Saving method")

  class(stats) <- c("summary.regions_results", class(stats))

  stats
}

#' @exportS3Method print summary.regions_results
print.summary.regions_results <- function(x, ...) {
  print.data.frame(x, row.names = FALSE)

  invisible(x)
}

.addregions_internal <- function(regions_results, noregions, exhaus = FALSE, cl = NULL, verbose = TRUE) {

  minvert <- attr(regions_results, "minvert")

  cont <- attr(regions_results, "cont")

  Xvar <- attr(regions_results, "pos")
  Yvar <- attr(regions_results, "scores")

  omitbp <- attr(regions_results, "omitbp")

  nbp <- noregions - 1

  noPC <- ncol(Yvar)

  colhead <- c("regions", paste0("breakpoint", seq_len(max(1, nbp))),
               "sumRSS", paste("RSS", 1:noPC, sep = "."))

  fregions <- function(BPs, Xvar, Yvar) {
    x <- .design_matrix(Xvar, BPs, cont)

    lines <- .lm.fit(x = x, y = Yvar)
    RSS <- sum(lines$residuals^2)

    rsq <- {
      if (noPC > 1) colSums(lines$residuals^2)
      else RSS
    }

    c((length(BPs) + 1), BPs, RSS, rsq)
  }

  if (nbp == 0) {
    # 1 Regions (0 bp)
    if (verbose) {
      cat("Fitting model with 1 region (0 breakpoints)... ")
    }

    lines <- .lm.fit(x = cbind(1, Xvar), y = Yvar)

    RSS <- sum(lines$residuals^2)

    rsq <- {
      if (noPC > 1) colSums(lines$residuals^2)
      else RSS
    }

    res <- as.data.frame(rbind(c(nbp + 1, NA_real_, RSS, rsq)))

    names(res) <- colhead

    bpkeep <- NA_character_

    stats <- data.frame(Nregions = nbp + 1L,
                        Nmodel_possible = 1L,
                        Nmodel_tested = 1L,
                        Nmodel_saved = 1L,
                        Comp_method = "Exhaustive",
                        Saving_method = "All")


    best_bps <- matrix(bpkeep, nrow = 1, dimnames = list(NULL, "Best_BPs1"))

    if (verbose) {
      cat("Done. 1 model tested, 1 model saved.\n")
    }
  }
  else if (nbp == 1) {
    # 2 Regions (1 bp)
    if (verbose) {
      cat("Fitting models with 2 regions (1 breakpoint)...\n")
    }

    # Get all possible combinations for 1 breakpoint:
    goodcomb <- .combosR(Xvar, minvert, nbp, omitbp = omitbp)

    # Run fitting on good combinations:
    res <- as.data.frame(do.call("rbind", pbapply::pbapply(goodcomb, 1, fregions, Xvar = Xvar,
                                                           Yvar = Yvar, simplify = FALSE, cl = cl)))
    names(res) <- colhead

    bpkeep <- {
      if (!exhaus) paste(goodcomb[,1], collapse = "|")
      else NA_character_
    }

    stats <- data.frame(Nregions = nbp + 1L,
                        Nmodel_possible = nrow(res),
                        Nmodel_tested = nrow(res),
                        Nmodel_saved = nrow(res),
                        Comp_method = "Exhaustive",
                        Saving_method = "All")

    best_bps <- matrix(bpkeep, nrow = 1, dimnames = list(NULL, "Best_BPs1"))

    if (verbose) {
      cat(sprintf("Done. %s models tested, %s models saved.\n",
                  nrow(res), nrow(res)))
    }
  }
  else {
    # 3 and more Regions (2+ bp)
    if (verbose) {
      cat(sprintf("Fitting models with %s regions (%s breakpoints)...\n",
                  nbp + 1, nbp))
    }

    # Define BPs to include if non-exhaustive search:
    if (!exhaus) {
      best_bps <- unlist(regions_results$stats[which.max(regions_results$stats$Nregions),
                                               startsWith(names(regions_results$stats), "Best_BPs")])
      bpkeep <- .drop_na(best_bps)
    }

    # Number of combinations
    ncombi <- .ncombos(Xvar, minvert, nbp)

    # Max length for an R object to not exceed 2Gb in size (Memory safety net)
    lim <- 5e8

    # Find nrow max for given number of cols (bps) to avoid exceeding 2Gb objects
    max_rows <- ceiling(lim/nbp)

    if (ncombi <= max_rows) {

      goodcomb <- .combosR(Xvar, minvert, nbp, omitbp = omitbp)

      nmodel_possible <- nrow(goodcomb)

      if (!exhaus) {
        # Keep only probable combinations (non-exhaustive search)
        for (m in strsplit(bpkeep, "|", fixed = TRUE)) {
          mn <- as.numeric(m)
          goodcomb <- goodcomb[apply(goodcomb, 1, function(b) any(b %in% mn)),, drop = FALSE]
        }
      }

      if (nrow(goodcomb) == 0) return(regions_results)

      res <- as.data.frame(do.call("rbind", pbapply::pbapply(goodcomb, 1, fregions, Xvar = Xvar,
                                                             Yvar = Yvar, simplify = FALSE, cl = cl)))

      names(res) <- colhead
      nmodel_tested <- nrow(res)

      if (!exhaus) {
        # If non-exhaustive search, keep only models with sumRSS <= min(sumRSS)+(0.5*sd(sumRSS))
        cutoff <- min(res$sumRSS) + sd(res$sumRSS) / 2
        if (!is.na(cutoff)) {
          # Subsample res only if more than 1 option
          res <- res[res$sumRSS <= cutoff,, drop = FALSE]
        }
      }
    }
    else {
      # Required number of splits of combos
      aa <- ceiling(ncombi / max_rows)

      # File names to store intermediate res objects
      files <- vapply(1:aa, function(i) tempfile("regions_", fileext = ".Rds"),
                      character(1L))

      # Initialize statistics
      nmodel_possible <- nmodel_tested <- 0
      last_combo <- NULL

      if (!exhaus) {
        min_RSS <- Inf
        n_total <- Xbar_total <- M2_total <- 0
      }

      # Initialize progress bar (pb). Won't activate if `verbose = FALSE`
      pb <- pbapply::startpb(0, ncombi)

      for (i in 1:aa) {
        # Only returns next max_rows rows of combos, using last_combo as
        # starting place
        goodcomb <- .combosR(Xvar, minvert, nbp, ncombos = max_rows,
                             last_combo = last_combo, omitbp = omitbp)

        last_combo <- goodcomb[nrow(goodcomb),]
        nmodel_possible <- nmodel_possible + nrow(goodcomb)

        if (!exhaus) {
          # Keep only probable combinations (non-exhaustive search)
          for (m in strsplit(bpkeep, "|", fixed = TRUE)) {
            mn <- as.numeric(m)
            goodcomb <- goodcomb[apply(goodcomb, 1, function(b) any(b %in% mn)),, drop = FALSE]
          }
        }

        if (nrow(goodcomb) == 0) next

        res <- as.data.frame(do.call("rbind", .lapply_selector(seq_len(nrow(goodcomb)), function(g) {
          on.exit(pbapply::setpb(pb, pbapply::getpb(pb) + 1))
          fregions(goodcomb[g,], Xvar = Xvar, Yvar = Yvar)
        }, cl = cl)))

        names(res) <- colhead
        nmodel_tested <- nmodel_tested + nrow(res)

        # Save intermediate object to disk
        saveRDS(res, files[i])

        if (!exhaus) {
          #Accumulate values for cutoff
          n_old <- n_total #old n
          n_i <- nrow(res) #n of new subset
          n_total <- n_old + n_i #total n so far

          Xbar_old <- Xbar_total #old mean
          Xbar_i <- sum(res$sumRSS)/n_i #mean of new subset
          Xbar_total <- (n_old * Xbar_old + n_i * Xbar_i)/(n_old + n_i) #new mean

          delta <- Xbar_i - Xbar_old

          M2_old <- M2_total #old M2
          M2_i <- sum((res$sumRSS - Xbar_i)^2) #M2 of new subset

          M2_total <- M2_old + M2_i + delta^2 * n_old * n_i / (n_old + n_i) #new M2

          min_RSS <- min(min_RSS, res$sumRSS) #new min
        }

      }

      if (!exhaus) {
        sd_RSS <- sqrt(M2_total/(nmodel_tested - 1))
        cutoff <- min_RSS + sd_RSS / 2
      }

      res <- NULL
      for (file in files) {
        res_tmp <- readRDS(file)
        if (!exhaus) {
          # If non-exhaustive search, keep only models with sumRSS <= min(sumRSS)+(0.5*sd(sumRSS))
          res_tmp <- res_tmp[res_tmp$sumRSS <= cutoff,, drop = FALSE]
        }
        res <- rbind(res, res_tmp)
        unlink(file)
      }

      pbapply::closepb(pb)
    }

    nmodel_saved <- nrow(res)

    # Create string with info on best BPs to use for next region:	# BPs kept are
    # all  > cutoff (min+0.5*sd) and BPs of best models +/- 3 vertebrae
    if (!exhaus) {
      bestBPs <- drop(res[which.min(res$sumRSS), startsWith(names(res), "breakpoint")])

      #Convert bp value to position, add range
      bbp <- lapply(bestBPs, function(x) match(x, Xvar) + seq(-3, 3))

      # Convert position of bp to keep to actual value of bp
      bestBPs <- lapply(bbp, function(x) Xvar[x])

      bpkeep <- vapply(seq_along(bestBPs), function(i) {
        paste(
          sort(.drop_na(unique(c(res[, i + 1], bestBPs[[i]]), nmax = length(Xvar)))),
          collapse = "|"
        )
      }, character(1L))

      stats <- data.frame(Nregions = noregions,
                          Nmodel_possible = nmodel_possible,
                          Nmodel_tested = nmodel_tested,
                          Nmodel_saved = nmodel_saved,
                          Comp_method = "Non-exhaus",
                          Saving_method = "SD/2")
    }
    else {
      bpkeep <- rep(NA_character_, nbp)

      stats <- data.frame(Nregions = noregions,
                          Nmodel_possible = nmodel_possible,
                          Nmodel_tested = nmodel_tested,
                          Nmodel_saved = nmodel_saved,
                          Comp_method = "Exhaustive",
                          Saving_method = "All")
    }

    best_bps <- matrix(bpkeep, nrow = 1, dimnames = list(NULL, paste0("Best_BPs", 1:nbp)))

    if (verbose) {
      cat(sprintf("Done. %s models tested, %s models saved.\n",
                  nmodel_tested, nmodel_saved))
    }
  }

  rownames(res) <- NULL

  stats <- cbind(stats, as.data.frame(best_bps))

  #Remove analysis with same number of regions
  if (any(regions_results$results$regions == noregions)) {
    regions_results$results <- regions_results$results[regions_results$results$regions != noregions,, drop = FALSE]
  }
  if (any(regions_results$stats$Nregions == noregions)) {
    regions_results$stats <- regions_results$stats[regions_results$stats$Nregions != noregions,, drop = FALSE]
  }

  #Combine res and stats with exist components of input object
  regions_results$results <- .rbind_larger(regions_results$results, res)
  regions_results$stats <- .rbind_larger(regions_results$stats, stats)

  regions_results$results <- regions_results$results[order(regions_results$results$regions),, drop = FALSE]
  regions_results$stats <- regions_results$stats[order(regions_results$stats$Nregions),, drop = FALSE]

  rownames(regions_results$results) <- NULL
  rownames(regions_results$stats) <- NULL

 regions_results
}

# Returns number of combos possible
.ncombos <- function(y, m, nbp) {
  if (nbp == 0) return(1)

  # Number of possibilities for each BP
  p <- length(y) - m * (nbp + 1) + 1
  if (p <= 0) return(0)

  # Number of combinations - this is magic, kind of, but it works
  prod(p - 1 + 1:nbp)/prod(1:nbp)
}

# Returns matrix of all combos that don't include omitbp
.combosR <- function(vert, m, nbp, ncombos = Inf, last_combo = NULL, omitbp = NULL) {
  n <- length(vert)

  y <- seq_len(n)

  mins <- m * seq_len(nbp)
  maxes <- rep(n, nbp) - rev(mins)

  maxcombos <- .ncombos(y, m, nbp)

  ncombos <- min(ncombos, maxcombos)
  combos <- matrix(NA_integer_, nrow = ncombos, ncol = max(nbp, 1))

  colnames(combos) <- paste0("breakpoint", seq_len(max(nbp, 1)))
  rownames(combos) <- NULL

  if (ncombos == 0 || nbp == 0) return(combos)

  if (is.null(last_combo)) {
    counters <- mins
  }
  else {
    counters <- match(last_combo, y)

    #Increment combo before inserting into combos matrix
    reset <- rep(FALSE, nbp)
    counters[nbp] <- counters[nbp] + 1

    for (i in rev(seq_len(nbp)[-1])) {
      if (counters[i] > maxes[i]) {
        reset[i] <- TRUE
        counters[i - 1] <- counters[i - 1] + 1
      }
    }
    if (any(reset)) {
      for (i in which(reset)) {
        counters[i] <- counters[i - 1] + m
      }
    }
  }

  reset <- rep(FALSE, nbp)
  k <- 1

  while (k <= nrow(combos) && counters[1] <= maxes[1]) {
    reset[] <- FALSE

    if (is.null(omitbp) || !any(vert[y[counters]] %in% omitbp)) {
      combos[k, ] <- y[counters]
      k <- k + 1
    }

    counters[nbp] <- counters[nbp] + 1

    for (i in rev(seq_len(nbp)[-1])) {
      if (counters[i] > maxes[i]) {
        reset[i] <- TRUE
        counters[i - 1] <- counters[i - 1] + 1
      }
    }
    if (any(reset)) {
      for (i in which(reset)) {
        counters[i] <- counters[i - 1] + m
      }
    }
  }

  combos <- combos[seq_len(k - 1),, drop = FALSE]
  combos[] <- vert[combos]

  combos
}

.design_matrix <- function(Xvar, BPs, cont) {
  ## Shoukd we use BPs + .5?
  if (length(BPs) == 0) {
    cbind(1, Xvar)
  }
  else if (cont) {
    # Continuous fit
    do.call("cbind", c(list(1, Xvar),
                       lapply(BPs, function(b) pmax(0, Xvar - b))))
  }
  else {
    # Discontinuous fit
    do.call("cbind", c(list(1 * (Xvar <= BPs[1]),
                            Xvar * (Xvar <= BPs[1])),
                       lapply(seq_along(BPs)[-1], function(i) {
                         cbind(1 * (Xvar > BPs[i - 1] & Xvar <= BPs[i]),
                               Xvar * (Xvar > BPs[i - 1] & Xvar <= BPs[i]))
                       }),
                       list(1 * (Xvar > BPs[length(BPs)]),
                            Xvar * (Xvar > BPs[length(BPs)]))))
  }
}
