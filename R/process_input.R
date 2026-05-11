#' Process vertebral measurements
#'
#' `process_measurements()` initializes the analysis workflow by processing a dataset of vertebra measurements into an object usable by \pkg{MorphoRegions}. Such processing includes identifying the vertebra indices and the measurements and filling in missing values.
#'
#' @param data a data.frame containing a column of vertebra indices and measurements for each vertebra, or a list thereof for multiple specimens.
#' @param pos the name or index of the variable in `data` containing the vertebra indices. Default is to use the first column.
#' @param measurements the names or indices of the variables in `data` containing the relevant vertebra measurements. If unspecified, will use all variables other than that specified in `pos`.
#' @param fillNA `logical`; whether to fill in missing values using a simple linear imputation. Default is `TRUE`. See Details.
#'
#' @returns A `regions_data` object, which is a list of data.frames (one for each specimen) with attributes containing metadata.
#'
#' @details
#' Any rows with missing values for all measurements will be removed. When missing values in non-removed rows are present and `fillNA` is set to `TRUE`, `process_measurements()` fills them in if the sequence of missing values is no greater than 2 in length. For numeric variables, it uses a linear interpolation, and for categorical variables, it fills in the missing values with the surrounding non-missing values if they are identical and leaves them missing otherwise. Otherwise, missing values are left as they are.
#'
#' When a list of data frames is supplied to `data`, only the variables named in `measurements` that are common across datasets will be stored as measurement variables.
#'
#' @seealso [svdPCO()] for computing principal coordinate axes from processed vertebra data.
#'
#' @example man/examples/example-process_measurements.R

#' @export
process_measurements <- function(data, pos = 1L, measurements, fillNA = TRUE) {

  if (is.matrix(data) || is.data.frame(data)) {
    data <- list(as.data.frame(data))
  }
  else if (is.list(data) && all(vapply(data, function(i) is.data.frame(i) || is.matrix(i), logical(1L)))) {
    for (i in seq_along(data)) {
      data[[i]] <- as.data.frame(data[[i]])
    }
  }
  else {
    chk::err("`data` must be a matrix, dataframe, or list thereof")
  }

  chk::chk_scalar(pos)
  chk::chk_flag(fillNA)

  pos_names <- character(length(data))
  pos_var_list <- vector("list", length(data))

  for (i in seq_along(data)) {
    if (chk::vld_whole_number(pos) && chk::vld_gte(pos, 1) && chk::vld_lte(pos, ncol(data[[i]]))) {
      pos_names[i] <- names(data[[i]])[as.integer(pos)]
    }
    else if (chk::vld_string(pos) && chk::vld_subset(pos, names(data[[i]]))) {
      pos_names[i] <- pos
    }
    else {
      chk::err("`pos` must be a single value indicating the column in `data` containing the vertebra positions")
    }

    pos_var_list[[i]] <- data[[i]][[pos_names[i]]]
  }

  if (length(unique(pos_names)) > 1) {
    chk::err("the variable identified by `pos` must have the same name in all supplied datasets")
  }

  pos <- pos_names[1]

  for (i in seq_along(pos_var_list)) {
    if (!chk::vld_whole_numeric(pos_var_list[[i]])) {
      chk::err("`pos` must refer to a variable of whole numbers identifying vertebra positions")
    }
  }

  measurements_names_list <- vector("list", length(data))

  for (i in seq_along(data)) {
    if (missing(measurements)) {
      measurements_names_list[[i]] <- setdiff(names(data[[i]]), pos)
    }
    else {
      if (length(measurements) == 0) {
        measurements_names_list[[i]] <- character()
      }
      else if (chk::vld_whole_numeric(measurements) && chk::vld_subset(measurements, seq_len(ncol(data[[i]])))) {
        measurements_names_list[[i]] <- names(data[[i]])[as.integer(measurements)]
      }
      else if (chk::vld_character(measurements) && chk::vld_subset(measurements, names(data[[i]]))) {
        measurements_names_list[[i]] <- measurements
      }
      else {
        chk::err("if supplied, `measurements` must indicate the columns in `data` containing the measurement values")
      }

      if (pos %in% measurements_names_list[[i]]) {
        chk::err("`pos` and `measurements` cannot overlap")
      }
    }
  }

  # Get only measurements that are common across datasets
  measurements <- Reduce(union, measurements_names_list)

  for (i in seq_along(data)) {

    #Subset to specified variables
    data[[i]] <- data[[i]][names(data[[i]]) %in% c(pos, measurements)]

    #Reorder columns to be consistent
    if (i > 1)
      data[[i]] <- data[[i]][names(data[[1]])]

    pos_ind <- match(pos, names(data[[i]]))

    #Rows that are all NA
    all_NA_rows <- which(apply(data[[i]], 1, function(x) all(is.na(x[-pos_ind]))))

    if (length(all_NA_rows) > 0)
      data[[i]] <- data[[i]][-all_NA_rows,]

    if (fillNA) {
      #Fill in missing values
      data[[i]] <- .missingval(data[[i]], pos_ind)

      if (anyNA(data[[i]])) {
        chk::wrn(sprintf("missing values remain in %s because there were sequences of missing values greater than %s in length", if (length(data) == 1) "the dataset" else paste("dataset", i), 2))
      }
    }
  }

  attr(data, "pos_ind") <- match(pos, names(data[[1]]))
  attr(data, "eligible_vertebrae") <- sort(unique(unlist(lapply(data, `[[`, pos))))

  class(data) <- "regions_data"

  data
}




#' Process vertebral measurements and PC scores of traditional morphometric data
#'
#' `process_PC()` initializes the analysis workflow by processing a dataset of vertebral measurements and already computed PC scores into an object usable by \pkg{MorphoRegions}. Such processing includes identifying the vertebra indices, measurements, and PC scores.
#'
#' @param data a data.frame containing a column of vertebra indices and measurements for each vertebra, or a named list thereof for multiple specimens with names of list corresponding to the unique identifier of each specimen.
#' @param pos the name or index of the variable in `data` containing the vertebra indices. Default is to use the first column.
#' @param pcscores a matrix or data.frame containing PC scores of each vertebra. Note that for multiple specimens, ordination method should have been performed on concatenated data across all specimens so that `pcscores` is a single data.frame or matrix containing scores of all specimens.
#' @param eigenvals a numeric vector containing eigenvalues of each PC axis.
#' @param posPC a vector containing the positional information of vertebrae in `pcscores`, or a named list thereof for multiple specimens with names of list corresponding to the unique identifier of each specimen. If not provided, function will assume row order of `pcscores` matches vertebra order in `pos`.
#'
#' @returns A `regions_pco` object, which contains user-provided eigenvectors in the `scores` component and eigenvalues in the `eigen.val` component. The original dataset, including positional information, is stored in the `data` attribute.
#'
#' @details
#' Unlike `process_measurements`, `process_PC` does not fill in missing values and these should have been removed or replaced by numeric values before running `process_PC`.
#'
#' @seealso
#' [process_gmPC()] for processing PC scores from 2D or 3D geometric morphometric datasets
#'
#' [svdPCO()] for computing principal coordinate axes from processed vertebra data.
#'
#' [plot.regions_pco()] for plotting PCO axes

#' @example man/examples/example-process_PC.R

#' @export



process_PC <- function(data, pos = 1L, pcscores, eigenvals, posPC){ #, specimen) {

  # Check if multiple specimens:
  multiSpec <- is.list(data) && !is.data.frame(data)

  # Check format of posPC and extract specimen names if multi-specimens:
  if(multiSpec && !missing(posPC)){
    chk::chk_list(posPC)
    chk::chk_named(posPC)
    chk::chk_named(data)
    vertposPC <- setNames(unlist(posPC),NULL)
    if(!chk::vld_whole_numeric(vertposPC)){
      chk::err("`posPC` must only contain whole numbers identifying vertebra positions")
    }
    spec_names_PC <- names(posPC)
    spec_names_data <- names(data)
    chk::chk_setequal(spec_names_PC, spec_names_data, x_name="specimen names in `pcscores`")
  }


  # Check format data and convert to list if needed:
  if (is.matrix(data) || is.data.frame(data)) {
    data <- list(as.data.frame(data))
  } else if (is.list(data) && all(vapply(data, function(i) is.data.frame(i) || is.matrix(i), logical(1L)))) {
    for (i in seq_along(data)) {
      data[[i]] <- as.data.frame(data[[i]])
    }
  } else {
    chk::err("`data` must be a matrix, dataframe, or list thereof")
  }

  # Check format pos:
  chk::chk_scalar(pos)

  pos_names <- character(length(data))
  pos_var_list <- vector("list", length(data))

  # Get positional info for each specimen:
  for (i in seq_along(data)) {
    if (chk::vld_whole_number(pos) && chk::vld_gte(pos, 1) && chk::vld_lte(pos, ncol(data[[i]]))) {
      pos_names[i] <- names(data[[i]])[as.integer(pos)]
    } else if (chk::vld_string(pos) && chk::vld_subset(pos, names(data[[i]]))) {
      pos_names[i] <- pos
    } else {
      chk::err("`pos` must be a single value indicating the column in `data` containing the vertebra positions")
    }
    pos_var_list[[i]] <- data[[i]][[pos_names[i]]]
  }

  # Check names and format of positional info:
  if (length(unique(pos_names)) > 1) {
    chk::err("the variable identified by `pos` must have the same name in all supplied datasets")
  }

  pos <- pos_names[1]

  for (i in seq_along(pos_var_list)) {
    if (!chk::vld_whole_numeric(pos_var_list[[i]])) {
      chk::err("`pos` must refer to a variable of whole numbers identifying vertebra positions")
    }
  }

  # Check format pcscores:
  if (is.matrix(pcscores)) {
    chk::chk_numeric(pcscores)
  } else if (is.data.frame(pcscores)) {
    chk::chk_data(pcscores)
    chk::chk_all(pcscores, chk::chk_numeric)
    pcscores <- as.matrix(pcscores)
  } else {
    chk::err("`pcscores` must be a numeric matrix or data.frame")
  }

  # Check eigenvals is numeric vector:
  chk::chk_vector(eigenvals)
  chk::chk_numeric(eigenvals)
  chk::chk_not_any_na(eigenvals)


  # Check format posPC and extract positional info from PCs:
  allvertData <- setNames(unlist(lapply(data,'[[', pos)), NULL)
  nvData <- length(allvertData)
  if (missing(posPC)) {
    chk::wrn("`posPC` not supplied: assuming order of vertebrae in `pcscores` matches the order of `data`")
    chk::chk_equal(nrow(pcscores),nvData, x_name="number of rows in `pcscores` must match number of vertebrae in `data`;")
    vertposPC <- allvertData
    spec_names_PC <- names(data)
  }

  if(!multiSpec && !missing(posPC)){
    chk::chk_vector(posPC)
    chk::chk_whole_numeric(posPC)
    chk::chk_equal(length(posPC),nvData, x_name="number of elements in `posPC` must match number of vertebrae in `data`;")
    vertposPC <- posPC
  }


  # Check number of PC pcscores and eigenvals match:
  if (length(eigenvals) != ncol(pcscores)) {
    chk::err(
      paste0("`eigenvals` must have one value per column in `pcscores`.\nLength of `eigenvals`: ", length(eigenvals),
             "; number of columns in `pcscores`: ", ncol(pcscores), ".")
    )
  }


  # Process data:
  # Note other transformations made to data in process_measurements have been removed here to prevent heavy modifications to data since it won't be used for downstream analyses
  for (i in seq_along(data)) {
    #Reorder columns to be consistent
    if (i > 1) {
      data[[i]] <- data[[i]][names(data[[1]])]
    }
    pos_ind <- match(pos, names(data[[i]]))
  }

  # Format 'data' part of output:
  attr(data, "pos_ind") <- match(pos, names(data[[1]]))
  attr(data, "eligible_vertebrae") <- sort(unique(unlist(lapply(data, `[[`, pos))))
  class(data) <- "regions_data"


  # Process pcscores:
  # Order vertebrae in pcscores according to order in data:
  if(setequal(vertposPC, allvertData)){
    if(!missing(posPC)){                                                # Reorder only if posPC is supplied
      if(multiSpec){
        vert_ord_data <- paste(rep(spec_names_data,lapply(data, nrow)),
                               allvertData, sep="_")                     # Create vector of vertebra order in data
        vert_ord_pc <- paste(rep(spec_names_PC, lapply(posPC,length)),
                             vertposPC, sep="_")                         # Create vector of vertebra order in pcscores
        pcscores <- pcscores[match(vert_ord_data, vert_ord_pc),]         # Re-order pcscores to match order in data
      } else {
        pcscores <- pcscores[match(allvertData, vertposPC),]         # Re-order pcscores to match order in data
      }
    }

  } else {
    chk::err("Vertebrae in `pcscores` must match vertebrae in included in `data`")
  }


  # Format PCA output:
  out <- list(scores=pcscores,
              eigen.val=eigenvals)
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
#' @param data a 3D array (p x k x n) of Procrustes aligned landmarks where p is the number of of landmark points per vertebra, k is the number of landmark dimensions (2D or 3D), and n is the number of vertebrae sampled for the specimen
#' @param pos either character string 'names' in which case the function will extract vertebral position info form the names of the array (using `dimnames[[3]]`), or a vector containing numeric values of vertebral position. Default to 'names'
#' @param pcscores a matrix or data.frame containing PC scores of each vertebra.
#' @param eigenvals a numeric vector containing eigenvalues of each PC axis.
#' @param specimens if data contains multiple specimens, a vector of specimen names of length equal to the number of vertebrae in `data`
#'
#' @returns A `regions_pco` object, which contains user-provided eigenvectors in the `scores` component and eigenvalues in the `eigen.val` component. The original dataset, including landmark coordinates and positional information, is stored in the `data` attribute.
#'
#' @details
#' Unlike `process_measurements`, `process_gmPC` does not fill in missing values and these should have been removed or replaced by numeric values before running `process_gmPC`.
#'
#' @seealso
#' [process_PC()] for processing PC scores from traditional morphometric datasets
#'
#' [svdPCO()] for computing principal coordinate axes from processed vertebra data.
#'
#' [plot.regions_pco()] for plotting PCO axes

#' @example man/examples/example-process_gmPC.R

#' @export



process_gmPC <- function(data, pos = 'names', pcscores, eigenvals, specimens) {

  # Check format of landmark data:
  chk::chk_array(data)

  # Check format of vertebral position info:
  if(is.vector(pos) && is.numeric(pos)){
    chk::chk_length(pos,dim(data)[3])       # If user-supplied vertebral position, check if has same length as number of vertebrae in array
    chk::chk_not_any_na(pos)                # If user-suppplied, check no NA value supplied
    chk::wrn("`pos` supplied as numeric vector: assuming its order matches the order of `data`")  # put warning that function cannot verify that order of pos and data match each other

  } else if(pos=='names'){
    raw_names <- dimnames(data)[[3]]
    chk::chk_not_null(raw_names)                     # Check no NULL values in names of array
    pos <- suppressWarnings(as.numeric(raw_names))       # If vertebra position info taken from names of array, convert them to numeric

    if(anyNA(pos)) {                           # Explicitly test for coercion-induced NAs:
      bad <- raw_names[is.na(pos)]
      chk::err(glue::glue("Cannot derive numeric vertebral positions from data.\n",
                          "The following names could not be converted to numbers:\n",
                          paste(bad, collapse = ", ")))
    }
  } else {
    chk::err("`pos` must be numeric or the character string `names`")
  }

  # Check if specimen provided:
  if(!missing(specimens)){
    chk::chk_length(specimens, dim(data)[3])  # Check if its length matches length of data
    specimens <- as.factor(specimens)
  } else {
    if(any(duplicated(pos))){
      chk::err("`specimens` must be provided when the same vertebral position is sampled more than once in `data` or `pos`")
    }
    specimens <- factor(rep('Specimen1', length(pos)))
  }

  # Check pcscores is a matrix or dataframe:
  if (is.matrix(pcscores)) {
    chk::chk_numeric(pcscores)
  } else if (is.data.frame(pcscores)) {
    chk::chk_data(pcscores)
    chk::chk_all(pcscores, chk::chk_numeric)
    pcscores <- as.matrix(pcscores)
  } else {
    chk::err("`pcscores` must be a numeric matrix or data.frame")
  }

  # Check eigenvals is numeric vector:
  chk::chk_vector(eigenvals)
  chk::chk_numeric(eigenvals)
  chk::chk_not_any_na(eigenvals)

  # Check number of PC pcscores and eigenvals match:
  if (length(eigenvals) != ncol(pcscores)) {
    chk::err(
      paste0("`eigenvals` must have one value per column in `pcscores`.\nLength of `eigenvals`: ", length(eigenvals),
             "; number of columns in `pcscores`: ", ncol(pcscores), ".")
    )
  }

  # Check number of vertebrae in pos and pcscores match:
  if (length(pos) != nrow(pcscores)) {
    chk::err("Number of vertebrae in `pcscores` must match number of vertebrae in `pos`.")
  }

  # Check rownames pcscores:
  if(is.null(rownames(pcscores))){
    chk::wrn("`pcscores` has no row names; assuming row order matches `pos`.")  # If no rownames, add pos as rownames
    rownames(pcscores) <- pos
  } else if (!setequal(rownames(pcscores), as.character(pos))){
    chk::wrn("Row names of `pcscores` do not match `pos`; assuming row order matches `pos`")  # If rownames different from pos, warning but do not replace original row names
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
  out <- list(scores=pcscores,
              eigen.val=eigenvals)
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
.missingval <- function(data, pos_ind, max_NA_seq_len = 2) {

  if (!anyNA(data[-pos_ind])) return(data)

  chk::chk_count(max_NA_seq_len)

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
    if (length(l.seq) == 0) next #if no short strings skip
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
