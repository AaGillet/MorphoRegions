#Functions from regions V1

#'Calculate PCO (principal co-ordinates analysis) based on svd
#'
#'Calculates distance matrix from raw data, then conducts a PCO ordination using a
#'single value decomposition.
#'
#'This differs from other PCO functions which use \code{CMDscale} that depends on a
#'spectral decomposition.
#'
#' @param x Raw data matrix.
#' @param metric Distance matrix calculation metric. "euclidean", "manhattan", or "gower".
#'
#' @return PCO scores (scores) and eigenvalues (eigen.val)
#' @export

svdPCO <- function(x, metric) {

  ### PCO analysis based on SVD not eigen as in CMDscale

  dis.gower <- cluster::daisy(x, metric = c(metric))  #set distance metric
  step2 <- (dis.gower^2)

  # Double centering code
  x <- as.matrix(step2)
  n <- dim(x)[1]
  k <- dim(x)[2]
  rowMeans <- matrix(apply(x, 1, mean, na.rm = TRUE), n, k, byrow = TRUE)
  colMeans <- matrix(apply(x, 2, mean, na.rm = TRUE), n, k, byrow = FALSE)
  matrixMean <- matrix(mean(x, na.rm = TRUE), n, k)
  step3 <- (x - rowMeans - colMeans + matrixMean)/-2

  # Single value decomp
  step4 <- svd(step3)
  step5 <- step4$v

  #step6<-scale(x%*%step5, center=T, scale=F)

  # Scale by root eigenval
  for (i in 1:ncol(step5)) {
    step5[, i] <- step5[, i] * sqrt(step4$d[i])
  }
  eigen.val <- step4$d
  eigen.vect <- step5[, 1:(ncol(step5) - 1)]  #remove the last PCO which has no variance

  return(list(scores = eigen.vect, eigen.val = eigen.val))
}


#' Find signficant PCOs
#'
#' Compares eigenvalue distribution of PCO to that with randomized data in order to extract
#' axes with significant signal.
#'
#' Significant PCO's are defined as those with greater eigenvalues than from random data. Therefore
#' the PCO cutoff is axes whose eigenvalues fall below the mean eigenvalue for that axis from the
#' randomized data. Data are randomly sampled by row.
#'
#' Warning: Bootstrapping is sensitive to unequale variances of columns. You may want to use
#' \code{scale} to normalize before using this approach
#'
#' @param data Raw data.
#' @param nreps Number of iterations.
#' @param metric Distance matrix calculation metric. "euclidean", "manhattan", or "gower".
#'
#'
#' @return Eigenvalues of original PCO (eigen.true), mean eigenvalues from the randomized data (eigen.mean),
#' standard deviation of randomized eigenvalues (eigen.sd), PCOs with eigenvalues exceeding the
#' randomized data (sigpco), matrix of all eigenvalues produced by randomization (eigen.boot).
#' @export

.PCOcutoff <- function(data, nreps, metric) {

  pcotrue <- svdPCO(data, metric)
  eigen.true <- pcotrue$eigen.val  #calculate 'true' eigenvalues
  eigen.true <- prop.table(eigen.true) #calculate as percentage variance

  eigen.boot <- NULL
  for (i in 1:nreps) {

    randdata <- data
    for (a in 1:nrow(data)) {
      randdata[a, ] <- sample(data[a, ])  #randomize the order of each row
    }

    if (length(which(is.na(cluster::daisy(randdata, metric = c("gower")))) > 0)) {
      # discard incomplete distance matrices
      i - 1
    } else {
      pco <- svdPCO(randdata, metric)  #calculate bootstrapped eigenvalues
      eigen.booti<-prop.table(pco$eigen.val) #calculate as percentage variance
      eigen.boot <- cbind(eigen.boot, eigen.booti)
    }
  }
  eigen.mean <- rowMeans(eigen.boot)  #calculate mean and SD of bootstrapped values
  eigen.sd <- apply(eigen.boot, 1, stats::sd)
  diff <- eigen.true - eigen.mean  #figure out which PCOs have greater eigenvalues for the 'true' dataset
  diff[diff < 0] <- 0
  sigpco <- length(split(diff, cumsum(diff == 0))$"0")  #split the dataset at the zeros, and calculate number of pcos in the first string


  return(list(eigen.true = eigen.true, eigen.mean = eigen.mean, eigen.sd = eigen.sd, sigpco = sigpco,
              eigen.boot = eigen.boot))

}

#' Plot eigenvalues
#'
#' Produces a scree plot of eigenvalues from a PCO and compares them with eigenvalues of
#' randomized data, obtained from \code{PCOcutoff}.
#'
#' @param eigen.true Vector of actual eigenvectors.
#' @param eigen.boot Matrix of bootstrapped eigenvectors.
#'
#' @importFrom graphics plot lines boxplot
#' @export

eigenplot <- function(eigen.true, eigen.boot) {

  plot(c(1:length(eigen.true)), eigen.true, main = "Eigenvalue cutoff",
       xlab = "PCO axis", ylab = "eigenvalue")
  lines(c(1:length(eigen.true)), eigen.true)
  boxplot(t(eigen.boot), xaxt = "n", add = T, outline = F)

}

#' Find the number of PCOs which yeilds the Maximum regionscore
#'
#' Searches top 10 PCOs to find the peak of regionscore
#'
#' @param regiondata Object returned from \code{compileregions}
#' @param noregions Maximum number of regions
#' @param nvert Number of vertebrae
#' @param nvar No of variables/pcos to evaluate. Default is all PCOs.
#'
#' @return pco.max No of PCOs producing the largest regionscore
#' @return pco.dist Region score for each cumulative addition of PCOs
#' @export
#'

#' Plot PCO axes
#'
#' Compare any PCO axes along the series.
#'
#' @param data PCO scores or any multivariate data matrix.
#' @param x Variable for x axis.
#' @param y Variable for y axis.
#' @param Xvar Variable for labelling points, usually positional variable such as vertebral count.
#' @export

axesplot <- function(data, x, y, Xvar) {

  plot(data[, x], data[, y], col = "white", xlab = "", ylab = "")
  graphics::text(data[, x], data[, y], labels = Xvar)
  graphics::title(main = "Principal coordinates analysis", xlab = paste("PCO", x), ylab = paste("PCO",
                                                                                                y))

}

#' Estimate and fill missing values
#'
#' Estimates missing values by averaging values from nearest adjacent vertebrae.
#'
#' In the case of a first or last vertebra, the nearest value is taken. Will only fill up to two missing data points.
#'
#' @param data data with missing elements
#'
#' @return data Data with missing values filled
#' @export
#'
#
Missingval<-function(data){

  if(any(is.na(data)==TRUE)){

    ###find strings of NAs with 2 or less missing
    for(i in 1:ncol(data)){ #for each variable
      dat<-data[,i]
      if(!any(is.na(dat)==TRUE)){ #Is there any missing data?
        next} else{
          miss.par<-which(is.na(dat))#find which ones are missing
          seqs<-split(miss.par, cumsum(c(1, diff(miss.par) != 1)))#split them into sequences
          if(length(seqs)==1){ #if theres one string
            if(length(unlist(seqs))>2){ #if the string is longer than 2 skip
              next}else{
                miss<-seqs#otherwise add to miss to fill
              }
          }
          if(length(seqs)>1){ #if theres more than one string
            l.seq<-which(sapply(seqs, length)<3)#which strings are two or less
            if(length(l.seq)==0){ next} #if no short strings skip
            seqs<-seqs[l.seq]
            miss<-seqs
          }

          for(a in 1:length(miss)){ #Fill each string
            fill<-unlist(miss[a])
            before<-min(fill)-1
            after<-max(fill)+1
            if(before<1){before<-after}#if at the beginning, use the end points
            if(after>length(dat)){after<-before}#if at the end, use beginning points
            val<-mean(c(dat[before],dat[after])) #calculate missing as mean of adjacent
            dat[fill]<-val #fill in the missing
          }

          data[,i]<-dat

        }
    }
    return(data)

  } else{return(data)}
}
