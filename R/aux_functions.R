# Clean version of match.arg() for a single choice
.match_arg <- function(arg, choices) {
  arg_name <- deparse1(substitute(arg))

  if (is.null(arg)) return(choices[1L])

  chk::chk_string(arg, arg_name)

  i <- pmatch(arg, choices, nomatch = 0L, duplicates.ok = TRUE)
  if (all(i == 0L)) {
    if (length(choices) == 2L) {
      chk::err("`", arg_name, "` should be one of ", paste(dQuote(choices, FALSE), collapse = " or "))
    }
    else {
      chk::err("`", arg_name, "` should be one of ", paste(dQuote(choices[-length(choices)], FALSE), collapse = ", "),
               ", or ", dQuote(choices[length(choices)], FALSE))
    }
  }

  choices[i]
}

# Similar to na.omit() without appending attributes
.drop_na <- function(x) {
  x[!is.na(x)]
}

# Similar to pbapply::pblapply() for handling cluster argument, but
# with no progress bar.
.lapply_selector <- function(..., cl = NULL) {
  if (is.null(cl)) {
    lapply(...)
  }
  else if (inherits(cl, "cluster")) {
    if (isTRUE(getOption("pboptions")$use_lb))
      parallel::parLapplyLB(cl, ...)
    else parallel::parLapply(cl, ...)
  }
  else {
    parallel::mclapply(..., mc.cores = as.integer(cl))
  }
}

# Similar to rbind() or dplyr::bind_rows(), but specifically to append a df
# to another with fewer columns
.rbind_larger <- function(x, y) {
  if (is.data.frame(x) || is.data.frame(y)) {
    x <- as.data.frame(x)
    y <- as.data.frame(y)

    nam <- {
      if (ncol(x) > ncol(y)) names(x)
      else names(y)
    }

    .expand <- function(d, nam) {
      nam_not_in_d <- setdiff(nam, names(d))
      to_add <- as.data.frame(matrix(nrow = nrow(d), ncol = length(nam_not_in_d),
                                     dimnames = list(rownames(d), nam_not_in_d)))
      return(cbind(d, to_add)[nam])
    }
  }
  else if (is.matrix(x) && is.matrix(y)) {
    nam <- {
      if (ncol(x) > ncol(y)) names(x)
      else names(y)
    }

    .expand <- function(d, nam) {
      nam_not_in_d <- setdiff(nam, colnames(d))
      to_add <- matrix(nrow = nrow(d), ncol = length(nam_not_in_d),
                       dimnames = list(rownames(d), nam_not_in_d))
      return(cbind(d, to_add)[, nam, drop = FALSE])
    }
  }
  else {
    stop("`x` and `y` must be data frames or matrices")
  }

  if (ncol(x) > ncol(y)) {
    rbind(x, .expand(y, nam))
  }
  else if (ncol(x) < ncol(y)) {
    rbind(.expand(x, nam), y)
  }
  else {
    rbind(x, y)
  }
}

.word_list <- function(word.list = NULL, and.or = "and") {
  #When given a vector of strings, creates a string of the form "a and b"
  #or "a, b, and c"

  word.list <- word.list[!word.list %in% c(NA_character_, "")]
  L <- length(word.list)

  if (L == 0) return("")

  if (L == 1) return(word.list)


  and.or <- .match_arg(and.or, c("and", "or"))
  if (L == 2) {
    out <- paste(word.list, collapse = paste0(" ", and.or, " "))
  }
  else {
    out <- paste(paste(word.list[seq_len(L - 1)], collapse = ", "),
                 word.list[L], sep = paste0(", ", and.or, " "))

  }

  return(out)
}

# Checks if supplied argument is a color; vectorized
.is_color <- function(x) {
  vapply(x, function(z) {
    tryCatch(is.matrix(grDevices::col2rgb(z)),
             error = function(e) FALSE)
    }, logical(1L))
}

#To pass CRAN checks:
utils::globalVariables(c("L.beg", "L.end", "L.pct.beg", "L.pct.end", "PCO",
                         "Region", "beg", "bp", "featuren", "ind", "ind.beg",
                         "ind.end", "pct.beg", "pct.end", "value", "vname"))