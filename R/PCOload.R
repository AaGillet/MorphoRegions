#' Calculate PCO loadings
#'
#' `PCOload()` computes the loadings for each principal coordinates (PCOs) analysis score, which are the correlations between the features used to compute the PCOs and the PCOs.
#'
#' @param x for `PCOload()`, a `regions_pco` object; the output of a call to [svdPCO()]. For `plot()`, a `regions_pco_load` object.
#' @param scores a numeric vector containing the indices of the desired scores.
#' @param \dots ignored.
#'
#' @return
#' `PCOload()` returns a `regions_pco_load` object, which is a matrix with a column for each PCO score requested and a row for each variable in the original dataset; values indicate the correlation between each variable and each PCO score. `plot()` returns a `ggplot` object, which can be manipulated using \pkg{ggplot2} syntax, that displays the loadings visually.
#'
#' @details
#' The loadings for a constructed variable, `vert.size`, are also computed and displayed. This is computed as the mean of the features for each vertebra. This function is not suitable for geometric morphometric data.
#'
#' @seealso
#' [svdPCO()] for computing the PCOs; [plot.regions_pco()] for visualizing the correlations between PCO scores.
#'
#' @example man/examples/example-PCOload.R
#'

#' @export
PCOload <- function(x, scores) {

  arg::arg_is(x, "regions_pco")
  pco_scores <- x[["scores"]]

  # New code: Check if original data is GM:
  if (inherits(attr(x, "data"), "regions_dataGM")) {
    arg::err("loadings cannot be computed for geometric morphometric data")
  }

  # New code: Check if PC scores have been computed outside of MorphoRegions:
  if (identical(attr(x, "metric"), "custom")) {
    arg::wrn("loadings computed as correlations between original data and PC scores; check if this method is relevant for your data and ordination method")
  }

  arg::when_supplied(
    scores,
    arg::arg_whole_numeric,
    arg::arg_between(c(1, ncol(pco_scores)))
  )

  if (missing(scores)) {
    scores <- seq_len(ncol(pco_scores))
  }

  data <- .get_data_without_pos(x)
  vert.size <- rowMeans(data)
  data <- cbind(data, vert.size)

  load.pco <- lapply(scores, function(i) {
    cor(data, pco_scores[, i], use = "pairwise.complete.obs")
  })

  load.pco <- do.call("cbind", load.pco)
  colnames(load.pco) <- paste("PCO", scores, sep = ".")

  class(load.pco) <- c("regions_pco_load", class(load.pco))

  load.pco
}

#' @exportS3Method print regions_pco_load
print.regions_pco_load <- function(x, digits = 3L, ...) {
  d <- as.data.frame.matrix(x)
  cat("- PCO loadings:\n\n")
  print(d[-nrow(x), , drop = FALSE], digits = digits, ...)
  cat("\n - Corr w/ vertebra size:\n\n")
  print(d[nrow(x), , drop = FALSE], digits = digits, ...)

  invisible(x)
}

#' @exportS3Method plot regions_pco_load
#' @rdname PCOload
plot.regions_pco_load <- function(x, ...) {
  d <- as.data.frame.matrix(x)
  rownames(d)[rownames(d) == "vert.size"] <- "Size"

  d$feature <- factor(rownames(d), levels = rev(rownames(d)))
  d$featuren <- as.numeric(d$feature)
  d$featuren[nrow(d)] <- .7

  d_long <- reshape(d, direction = "long", varying = startsWith(names(d), "PCO"),
                    idvar = "feature", v.names = "value",
                    timevar = "PCO")

  ggplot(d_long) +
    geom_tile(aes(x = .data$PCO, y = .data$featuren, fill = .data$value)) +
    labs(y = NULL, x = "PCO", fill = "Loading") +
    scale_x_continuous(position = "top", breaks = seq_len(max(d_long$PCO))) +
    scale_y_continuous(labels = rev(rownames(d)),
                       breaks = sort(unique(d$featuren)),
                       expand = c(0, 0)) +
    scale_fill_gradient2(limits = c(-1, 1), low = scales::muted("blue"),
                         high = scales::muted("red")) +
    coord_fixed() +
    theme_minimal() +
    theme(panel.grid = element_blank())
}
