#' Process vertebra measurements
#'
#' `process_measurements()` initializes the analysis workflow by processing a dataset of vertebra measurements into an object usable by `regions`. Such processing includes identifying the vertebra indices and the measurements and filling in missing values.
#'
#' @param data a dataset containing a column of vertebra indices and measurements for each vertebra.
#' @param pos the name or index of the variable in `data` containing the vertebra indices. Default is to use the first column.
#' @param measurements the names or indices of the variables in `data` containing the relevant vertebra measurements. If unspecified, will use all variables other than that specified in `pos`.
#'
#' @returns A `regions_data` object, which is a data.frame with attributes containing metadata. The vertebra index variable is removed from the data and stored as an attribute.
#'
#' @details
#' When missing values are present, `process_measurements()` fills them in using the mean of the surrounding non-missing values (i.e., a linear interpolation) if the sequence of missing values is no greater than 2 in length. Otherwise, missing values are left as they are.
#'
#' @seealso [svdPCO()] for computing principal coordinate axes from processed vertebra data.
#'
#' @example man/examples/example-process_measurements.R

#' @export
process_measurements <- function(data, pos = 1L, measurements = NULL) {
  if (is.matrix(data)) data <- as.data.frame(data)

  chk::chk_data(data)
  chk::chk_scalar(pos)

  if (chk::vld_whole_number(pos) && chk::vld_gte(pos, 1) && chk::vld_lte(pos, ncol(data))) {
    pos <- as.integer(pos)
  }
  else if (chk::vld_string(pos) && chk::vld_subset(pos, names(data))) {
    pos <- match(pos, names(data))
  }
  else {
    chk::err("`pos` must be a single value indicating the column in `data` containing the vertebra positions")
  }

  pos_var <- data[[pos]]

  if (!is.numeric(pos_var)) {
    chk::err("`pos` must refer to a numeric variable identifying vertebra positions")
  }

  if (!is.null(measurements)) {
    if (chk::vld_whole_numeric(measurements) && chk::vld_subset(measurements, seq_len(ncol(data)))) {
      measurements <- as.integer(measurements)
    }
    else if (chk::vld_character(measurements) && chk::vld_subset(measurements, names(data))) {
      measurements <- match(measurements, names(data))
    }
    else {
      chk::err("if supplied, `measurements` must indicate the columns in `data` containing the measurement values")
    }
  }
  else {
    measurements <- seq_len(ncol(data))[-pos]
  }

  if (pos %in% measurements) {
    chk::err("`pos` and `measurements` cannot overlap")
  }

  data <- data[measurements]

  if (any(!vapply(data, is.numeric, logical(1L)))) {
    chk::err("all measurements must be numeric")
  }

  ord <- order(pos_var)
  data <- data[ord,]
  pos_var <- pos_var[ord]

  data <- .missingval(data)

  if (anyNA(data)) {
    chk::wrn("missing values remain in the dataset because there were sequences of missing values greater than 2 in length")
  }

  attr(data, "pos") <- pos_var
  attr(data, "subset") <- seq_len(nrow(data))
  class(data) <- c("regions_data", class(data))

  data
}

.missingval <- function(data) {
  if (anyNA(data)) {

    ###find strings of NAs with 2 or less missing
    for (i in which(vapply(data, anyNA, logical(1L)))) { #for each variable with NA
      dat <- data[[i]]

      miss.par <- which(is.na(dat)) #find which ones are missing
      seqs <- split(miss.par, cumsum(c(1, diff(miss.par) != 1))) #split them into sequences

      l.seq <- which(lengths(seqs) <= 2) #which strings are two or less
      if (length(l.seq) == 0) next #if no short strings skip
      seqs <- seqs[l.seq]

      for (fill in seqs) { #Fill each string
        before <- min(fill) - 1
        after <- max(fill) + 1
        if (before < 1) before <- after #if at the beginning, use the end points
        if (after > length(dat)) after <- before #if at the end, use beginning points
        val <- mean(c(dat[before], dat[after])) #calculate missing as mean of adjacent
        dat[fill] <- val #fill in the missing
      }

      data[[i]] <- dat
    }
  }
  return(data)
}
