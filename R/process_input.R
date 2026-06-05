#' Process vertebral measurements
#'
#' `process_measurements()` initializes the analysis workflow by processing a dataset of vertebra measurements into an object usable by \pkg{MorphoRegions}. Such processing includes identifying the vertebra indices and the measurements and filling in missing values.
#'
#' @param data a data frame containing a column of vertebra indices and measurements for each vertebra, or a list thereof for multiple specimens.
#' @param pos the name or index of the variable in `data` containing the vertebra indices. Default is to use the first column.
#' @param measurements the names or indices of the variables in `data` containing the relevant vertebra measurements. If unspecified, will use all variables other than that specified in `pos`.
#' @param fillNA `logical`; whether to fill in missing values using a simple linear imputation. Default is `TRUE`. See Details.
#'
#' @return
#' A `regions_data` object, which is a list of data frames (one for each specimen) with attributes containing metadata.
#'
#' @details
#' Any rows with missing values for all measurements will be removed. When missing values in non-removed rows are present and `fillNA` is set to `TRUE`, `process_measurements()` fills them in if the sequence of missing values is no greater than 2 in length. For numeric variables, it uses a linear interpolation, and for categorical variables, it fills in the missing values with the surrounding non-missing values if they are identical and leaves them missing otherwise. Otherwise, missing values are left as they are.
#'
#' When a list of data frames is supplied to `data`, only the variables named in `measurements` that are common across datasets will be stored as measurement variables.
#'
#' @seealso
#' [svdPCO()] for computing principal coordinate axes from processed vertebra data.
#'
#' @example man/examples/example-process_measurements.R

#' @export
process_measurements <- function(data, pos = 1L, measurements, fillNA = TRUE) {
  arg::arg_supplied(data)

  if (is.matrix(data) || is.data.frame(data)) {
    data <- list(as.data.frame(data))
  }
  else if (is.list(data) && all(vapply(data, function(i) is.data.frame(i) || is.matrix(i), logical(1L)))) {
    for (i in seq_along(data)) {
      data[[i]] <- as.data.frame(data[[i]])
    }
  }
  else {
    arg::err("{.arg data} must be a matrix, data frame, or list thereof")
  }

  arg::arg_flag(fillNA)

  pos_names <- character(length(data))
  pos_var_list <- vector("list", length(data))

  for (i in seq_along(data)) {
    arg::arg_index(pos, data[[i]],
                   .msg = "{.arg pos} must be a single value indicating the column in {.arg data} containing the vertebra positions")

    if (is.numeric(pos)) {
      pos_names[i] <- names(data[[i]])[as.integer(pos)]
    }
    else {
      pos_names[i] <- pos
    }

    pos_var_list[[i]] <- data[[i]][[pos_names[i]]]
  }

  if (length(unique(pos_names)) > 1L) {
    arg::err("the variable identified by {.arg pos} must have the same name in all supplied datasets")
  }

  pos <- pos_names[1]

  for (i in seq_along(pos_var_list)) {
    arg::arg_whole_numeric(pos_var_list[[i]],
                           .msg = "{.arg pos} must refer to a variable of whole numbers identifying vertebra positions")
  }

  measurements_names_list <- vector("list", length(data))

  for (i in seq_along(data)) {
    arg::when_supplied(
      measurements,
      arg::arg_or(
        arg::arg_indices(data[[i]], .arg_data = "data"),
        arg::arg_null
      )
    )

    if (missing(measurements)) {
      measurements_names_list[[i]] <- setdiff(names(data[[i]]), pos)
    }
    else if (length(measurements) == 0) {
      measurements_names_list[[i]] <- character()
    }
    else if (is.numeric(measurements)) {
      measurements_names_list[[i]] <- names(data[[i]])[as.integer(measurements)]
    }
    else {
      measurements_names_list[[i]] <- measurements
    }

    if (pos %in% measurements_names_list[[i]]) {
      arg::err("{.arg pos} and {.arg measurements} cannot overlap")
    }
  }

  # Get only measurements that are common across datasets
  measurements <- Reduce(union, measurements_names_list)

  for (i in seq_along(data)) {

    #Subset to specified variables
    data[[i]] <- data[[i]][names(data[[i]]) %in% c(pos, measurements)]

    #Reorder columns to be consistent
    if (i > 1L) {
      data[[i]] <- data[[i]][names(data[[1L]])]
    }

    pos_ind <- match(pos, names(data[[i]]))

    #Rows that are all NA
    all_NA_rows <- which(apply(data[[i]], 1L, function(x) all(is.na(x[-pos_ind]))))

    if (length(all_NA_rows) > 0) {
      data[[i]] <- data[[i]][-all_NA_rows,]
    }

    if (fillNA) {
      #Fill in missing values
      data[[i]] <- .missingval(data[[i]], pos_ind)

      if (anyNA(data[[i]])) {
        if (length(data) == 1L) {
          arg::wrn("missing values remain in the dataset because there were sequences of missing values greater than 2 in length")
        }
        else {
          arg::wrn("missing values remain in dataset {i} because there were sequences of missing values greater than 2 in length")
        }
      }
    }
  }

  attr(data, "pos_ind") <- match(pos, names(data[[1L]]))
  attr(data, "eligible_vertebrae") <- sort(unique(unlist(lapply(data, `[[`, pos))))

  class(data) <- "regions_data"

  data
}


#' Process vertebral measurements and PC scores of traditional morphometric data
#'
#' `process_PC()` initializes the analysis workflow by processing a dataset of vertebral measurements and already computed PC scores into an object usable by \pkg{MorphoRegions}. Such processing includes identifying the vertebra indices, measurements, and PC scores.
#'
#' @inheritParams process_measurements
#' @param data a data frame containing a column of vertebra indices and measurements for each vertebra, or a named list thereof for multiple specimens with the names of the list corresponding to the unique identifiers of each specimen.
#' @param pcscores a matrix or data frame containing PC scores of each vertebra. Note that for multiple specimens, an ordination method should have been performed on the concatenated data across all specimens so that `pcscores` is a single data frame or matrix containing scores for all specimens.
#' @param eigenvals a numeric vector containing the eigenvalues of each PC axis.
#' @param posPC a vector containing the positional information of vertebrae in `pcscores`, or a named list thereof for multiple specimens with the names of the list corresponding to the unique identifiers of each specimen. If not provided, will assume the row order of `pcscores` matches the vertebra order in `pos`.
#'
#' @return
#' A `regions_pco` object, which contains user-provided eigenvectors in the `scores` component and eigenvalues in the `eigen.val` component. The original dataset, including positional information, is stored in the `data` attribute.
#'
#' @details
#' Unlike `process_measurements()`, `process_PC()` does not fill in missing values and these should have been removed or replaced by numeric values before running `process_PC()`.
#'
#' @seealso
#' * [process_gmPC()] for processing PC scores from 2D or 3D geometric morphometric datasets
#' * [svdPCO()] for computing principal coordinate axes from processed vertebra data.
#' * [plot.regions_pco()] for plotting PCO axes

#' @example man/examples/example-process_PC.R

#' @export
process_PC <- function(data, pos = 1L, pcscores, eigenvals, posPC){ #, specimen) {
  arg::arg_supplied(data)
  arg::arg_supplied(pcscores)
  arg::arg_supplied(eigenvals)

  # Check if multiple specimens:
  multiSpec <- is.list(data) && !is.data.frame(data)

  # Check format of posPC and extract specimen names if multi-specimens:
  if (multiSpec && !missing(posPC)) {
    arg::arg_list(posPC)
    arg::arg_named(posPC)
    arg::arg_named(data)

    vertposPC <- setNames(unlist(posPC), NULL)

    arg::arg_whole_numeric(vertposPC,
                           .msg = "{.arg posPC} must only contain whole numbers identifying vertebra positions")

    spec_names_PC <- names(posPC)
    spec_names_data <- names(data)

    if (!setequal(spec_names_PC, spec_names_data)) {
      arg::err("the specimen names in {.arg pcscores} must be the same as those in {.arg data}")
    }
  }

  # Check format data and convert to list if needed:
  if (is.matrix(data) || is.data.frame(data)) {
    data <- list(as.data.frame(data))
  }
  else if (is.list(data) && all(vapply(data, function(i) is.data.frame(i) || is.matrix(i), logical(1L)))) {
    for (i in seq_along(data)) {
      data[[i]] <- as.data.frame(data[[i]])
    }
  }
  else {
    arg::err("{.arg data} must be a matrix, data frame, or list thereof")
  }

  pos_names <- character(length(data))
  pos_var_list <- vector("list", length(data))

  # Get positional info for each specimen:
  for (i in seq_along(data)) {
    arg::arg_index(pos, data[[i]],
                   .msg = "{.arg pos} must be a single value indicating the column in {.arg data} containing the vertebra positions")

    if (is.numeric(pos)) {
      pos_names[i] <- names(data[[i]])[as.integer(pos)]
    }
    else {
      pos_names[i] <- pos
    }

    pos_var_list[[i]] <- data[[i]][[pos_names[i]]]
  }

  if (length(unique(pos_names)) > 1L) {
    arg::err("the variable identified by {.arg pos} must have the same name in all supplied datasets")
  }

  pos <- pos_names[1L]

  for (i in seq_along(pos_var_list)) {
    arg::arg_whole_numeric(pos_var_list[[i]],
                           .msg = "{.arg pos} must refer to a variable of whole numbers identifying vertebra positions")
  }

  # Check pcscores
  if (is.data.frame(pcscores) && all(vapply(pcscores, is.numeric, logical(1L)))) {
    pcscores <- as.matrix(pcscores)
  }
  else if (!is.matrix(pcscores) || !is.numeric(pcscores)) {
    arg::err("{.arg pcscores} must be a numeric matrix or a data frame containing only numeric columns")
  }

  # Check eigenvals
  arg::arg_vector(eigenvals)
  arg::arg_numeric(eigenvals)
  arg::arg_no_NA(eigenvals)

  # Check format posPC and extract positional info from PCs:
  allvertData <- setNames(unlist(lapply(data, `[[`, pos)), NULL)

  nvData <- length(allvertData)

  if (missing(posPC)) {
    arg::wrn("{.arg posPC} not supplied: assuming order of vertebrae in {.arg pcscores} matches the order of {.arg data}")

    if (nrow(pcscores) != nvData) {
      arg::err("the number of rows in {.arg pcscores} must match the number of vertebrae in {.arg data}")
    }

    vertposPC <- allvertData
    spec_names_PC <- names(data)
  }
  else if (!multiSpec) {
    arg::arg_vector(posPC)
    arg::arg_whole_numeric(posPC)

    if (length(posPC) != nvData) {
      arg::err("the number of elements in {.arg posPC} must match the number of vertebrae in {.arg data}")
    }

    vertposPC <- posPC
  }

  # Check number of PC pcscores and eigenvals match:
  if (length(eigenvals) != ncol(pcscores)) {
    arg::err(c("{.arg eigenvals} must have one value per column in {.arg pcscores}.",
               "i" = "Length of {.arg eigenvals}: {length(eigenvals)}",
               "i" = "Number of columns in {.arg pcscores}: {ncol(pcscores)}"))
  }


  # Process data:
  # Note other transformations made to data in process_measurements have been removed here to prevent heavy modifications to data since it won't be used for downstream analyses
  for (i in seq_along(data)) {
    #Reorder columns to be consistent
    if (i > 1L) {
      data[[i]] <- data[[i]][names(data[[1L]])]
    }

    pos_ind <- match(pos, names(data[[i]]))
  }

  # Format 'data' part of output:
  attr(data, "pos_ind") <- match(pos, names(data[[1L]]))
  attr(data, "eligible_vertebrae") <- sort(unique(unlist(lapply(data, `[[`, pos))))
  class(data) <- "regions_data"

  # Process pcscores:
  # Order vertebrae in pcscores according to order in data:
  if (!setequal(vertposPC, allvertData)) {
    arg::err("the vertebrae in {.arg pcscores} must match those included in {.arg data}")
  }

  if (!missing(posPC)) {                                                # Reorder only if posPC is supplied
    if (multiSpec) {
      vert_ord_data <- paste(rep(spec_names_data, unlist(lapply(data, nrow))),
                             allvertData, sep = "_")                     # Create vector of vertebra order in data
      vert_ord_pc <- paste(rep(spec_names_PC, lengths(posPC)),
                           vertposPC, sep = "_")                         # Create vector of vertebra order in pcscores
      pcscores <- pcscores[match(vert_ord_data, vert_ord_pc), , drop = FALSE]         # Re-order pcscores to match order in data
    }
    else {
      pcscores <- pcscores[match(allvertData, vertposPC), , drop = FALSE]         # Re-order pcscores to match order in data
    }
  }

  # Format PCA output:
  out <- list(scores = pcscores,
              eigen.val = eigenvals)

  attr(out,"data") <- data
  attr(out, "metric") <- "custom"
  attr(out, "scale") <- "custom"
  attr(out, "specimen") <- factor(rep(seq_along(data), unlist(lapply(data, nrow))),
                                  labels = paste("Specimen", seq_along(data)),
                                  levels = seq_along(data))
  class(out) <- "regions_pco"

  out
}

#' Process vertebral measurements and PC scores of 2D or 3D geometric morphometric data
#'
#' `process_gmPC()` initializes the analysis workflow by processing vertebral measurements and already computed PC scores from 2D or 3D geometric morphometric datasets into an object usable by \pkg{MorphoRegions}. Such processing includes identifying the vertebra indices, measurements, and PC scores.
#'
#' @inheritParams process_PC
#' @param data a 3D array (\eqn{p \times k \times n}) of Procrustes aligned landmarks where \eqn{p} is the number of of landmark points per vertebra, \eqn{k} is the number of landmark dimensions (2D or 3D), and \eqn{n} is the number of vertebrae sampled for the specimen.
#' @param pos either character string `"names"` in which case the function will extract vertebral position info form the names of the array (using `dimnames(data)[[3]]`), or a vector containing numeric values of vertebral position. Default is `"names"`.
#' @param pcscores a matrix or data frame containing PC scores of each vertebra.
#' @param specimens if `data` contains multiple specimens, a vector of specimen names of length equal to the number of vertebrae in `data`.
#'
#' @return
#' A `regions_pco` object, which contains user-provided eigenvectors in the `scores` component and eigenvalues in the `eigen.val` component. The original dataset, including landmark coordinates and positional information, is stored in the `data` attribute.
#'
#' @details
#' Unlike `process_measurements()`, `process_gmPC()` does not fill in missing values and these should have been removed or replaced by numeric values before running `process_gmPC()`.
#'
#' @seealso
#' * [process_PC()] for processing PC scores from traditional morphometric datasets
#' * [svdPCO()] for computing principal coordinate axes from processed vertebra data.
#' * [plot.regions_pco()] for plotting PCO axes

#' @example man/examples/example-process_gmPC.R

#' @export
process_gmPC <- function(data, pos = "names", pcscores, eigenvals, specimens = NULL) {

  # Check format of landmark data:
  if (!is.array(data) || !is.numeric(data) || length(dim(data)) != 3L) {
    arg::err("{.arg data} must be a 3-dimensional numeric array")
  }

  # Check format of vertebral position info:
  arg::arg_or(
    pos,
    arg::arg_element("names"),
    arg::arg_numeric
  )

  if (identical(pos, "names")) {
    raw_names <- dimnames(data)[[3L]]

    # Check no NULL values in names of array
    if (length(raw_names) == 0L) {
      arg::err("when {.arg pos} is {.val names}, {.code dimnames(data)[[3]]} cannot be empty")
    }

    # If vertebra position info taken from names of array, convert them to numeric
    pos <- suppressWarnings(as.numeric(raw_names))

    # Explicitly test for coercion-induced NAs
    if (anyNA(pos)) {
      arg::err(c("Cannot derive numeric vertebral positions from data.",
                 "i" = "The following names could not be converted to numbers: {.val {raw_names[is.na(pos)]}}"))
    }
  }
  else if (length(pos) == dim(data)[3L]) {
    arg::arg_no_NA(pos)
    # Warn that function cannot verify that order of pos and data match each other
    arg::wrn("{.arg pos} supplied as numeric vector; assuming its order matches the order of {.arg data}")
    pos <- as.vector(pos)
  }
  else {
    # If user-supplied vertebral position, check if has same length as number of vertebrae in array
    arg::err("when supplied as a numeric vector, {.arg pos} must have length equal to that of {.code dim(data)[3]}")
  }

  # Check if specimen provided:
  if (!missing(specimens)) {
    if (length(specimens) != dim(data)[3L]) {
      arg::err("when supplied, {.arg specimens} must have length equal to {.code dim(data)[3]}")
    }

    specimens <- as.factor(specimens)
  }
  else if (anyDuplicated(pos) == 0L) {
    specimens <- factor(rep.int("Specimen1", length(pos)), nmax = 1L)
  }
  else {
    arg::err("{.arg specimens} must be provided when the same vertebral position is sampled more than once in {.arg data} or {.arg pos}")
  }

  # Check pcscores
  if (is.data.frame(pcscores) && all(vapply(pcscores, is.numeric, logical(1L)))) {
    pcscores <- as.matrix(pcscores)
  }
  else if (!is.matrix(pcscores) || !is.numeric(pcscores)) {
    arg::err("{.arg pcscores} must be a numeric matrix or a data frame containing only numeric columns")
  }

  # Check eigenvals
  arg::arg_vector(eigenvals)
  arg::arg_numeric(eigenvals)
  arg::arg_no_NA(eigenvals)

  # Check number of PC pcscores and eigenvals match:
  if (length(eigenvals) != ncol(pcscores)) {
    arg::err(
      c("{.arg eigenvals} must have one value per column in {.arg pcscores}.",
      "i" = "Length of {.arg eigenvals}: {.val {length(eigenvals)}}",
      "i" = "Columns in {.arg pcscores}: {.val {ncol(pcscores)}}")
    )
  }

  # Check number of vertebrae in pos and pcscores match:
  if (length(pos) != nrow(pcscores)) {
    arg::err("the number of vertebrae in {.arg pcscores} must match number of vertebrae in {.arg pos}")
  }

  # Check rownames pcscores:
  if (length(rownames(pcscores)) == 0L) {
    # If no rownames, add pos as rownames
    arg::wrn("{.arg pcscores} has no row names; assuming row order matches {.arg pos}")
    rownames(pcscores) <- pos
  }
  else if (!setequal(rownames(pcscores), as.character(pos))) {
    # If rownames different from pos, warning but do not replace original row names
    arg::wrn("the row names of {.arg pcscores} do not match {.arg pos}; assuming row order matches {.arg pos}")
  }

  # # Order pos, data, and pcscores from most anterior to most posterior vertebra:
  # ord <- order(pos)
  # pos <- pos[ord]
  # data <- data[,,ord, drop=F]
  # pcscores <- pcscores[ord,, drop=F]


  # Format 'data' part of output:
  attr(data, "Xvar") <- pos
  attr(data, "eligible_vertebrae") <- sort(unique(pos))
  class(data) <- "regions_dataGM"

  # Format PCA output:
  out <- list(scores = pcscores,
              eigen.val = eigenvals)

  attr(out,"data") <- data
  attr(out, "metric") <- "custom"
  attr(out, "scale") <- "custom"
  attr(out, "specimen") <- specimens

  class(out) <- "regions_pco"

  out
}

# Takes in a dataframe and fills in missing values
# `max_NA_seq_len` controls max number of consecutive missing values that can
# be imputed
.missingval <- function(data, pos_ind, max_NA_seq_len = 2L) {

  if (!anyNA(data[-pos_ind])) {
    return(data)
  }

  arg::arg_count(max_NA_seq_len)

  # Add missing vertebrae
  missing_vertebrae <- setdiff(seq(min(data[[pos_ind]]),
                                   max(data[[pos_ind]])),
                               data[[pos_ind]])

  if (length(missing_vertebrae) > 0) {
    missing_rows <- as.data.frame(matrix(NA, nrow = length(missing_vertebrae),
                                         ncol = ncol(data)))
    names(missing_rows) <- names(data)
    missing_rows[[pos_ind]] <- missing_vertebrae

    data <- rbind(data, missing_rows)
  }

  # Order by vertebra
  o <- order(data[[pos_ind]])

  pos <- data[[pos_ind]]

  # Find strings of 1-2 NAs
  variables_with_missing <- which(vapply(data[!(pos %in% missing_vertebrae), -pos_ind],
                                         anyNA, logical(1L)))

  for (i in variables_with_missing) {
    # Extract variable ordered by vertebra
    dat <- data[o, -pos_ind][[i]]

    miss.par <- which(is.na(dat)) #find which ones are missing
    seqs <- split(miss.par, cumsum(c(1, diff(miss.par) != 1))) #split them into sequences

    l.seq <- which(lengths(seqs) <= max_NA_seq_len) #which strings are two or less

    if (length(l.seq) == 0) {
      #if no short strings skip
      next
    }

    seqs <- seqs[l.seq]

    for (fill in seqs) { #Fill each string
      before <- min(fill) - 1
      after <- max(fill) + 1
      if (before < 1) before <- after #if at the beginning, use the end points
      if (after > length(dat)) after <- before #if at the end, use beginning points

      if (!is.numeric(dat) && dat[before] != dat[after]) next

      if (is.numeric(dat)) {
        val <- seq(dat[before], dat[after], length.out = length(fill) + 2) #calculate missing as mean of adjacent
        dat[fill] <- val[-c(1, length(val))] #fill in the missing
      }
      else {
        val <- rep(dat[before], length(fill))
        dat[fill] <- val #fill in the missing
      }
    }

    data[o, -pos_ind][[i]] <- dat
  }

  # Remove extra fully missing rows
  subset(data, !pos %in% missing_vertebrae)
}
