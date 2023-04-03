
###############
### CODE AG ###
###############



#############################
#############################
## Computational functions ##
#############################
#############################



###########
## Function: calculate regions ##
###########
	# Using tibble instead of dataframe



# ERROR to correct for non-parallel computing: 
	# issues with environments in 3+ regions calculations for fcomb and fregions, fregions do not find Yvar...

# Note: for this function and addregionsAG, code changed so that function doesn't return R² value of each PCO anymore but RSS instead (needed for the plotpcoregAG function)
	# but object containing RSS of each PCO still named "rsq" for simplicity


calcregionsAG <- function(Xvar, Yvar, noregions, minvert=3, cont=T, exhaus=F, par=T, numCores=2, verbose=T){
		# Xvar: vector with positional vertebra info (ie vertebral number)
		# Yvard: 
		# noregions: maximal number of regions wants to try (eg, if noregions=9, will test all options from 1 to 9 regions)
		# minvert: minimal number of vertebrae / region, if not provided, default is 3
		# cont: TRUE/FALSE - choose if want to use continuous (TRUE) or discontinuous segmentation fitting, default to TRUE
		# exhaus: TRUE/FALSE - choose if want to use an exhaustive search of all possible breakpoints (TRUE) or non-exhaustive (FALSE), default to F
		#		will automatically use the exhaustive search for 1 and 2 region analyses (ie, 0 or 1 breakpoint)
		# par: TRUE/FALSE - choose if want to use parallel computing, default to FALSE
  # Default should be changed to FALSE but noont working for now!!
		# numCores: if par=TRUE, starting from 3 regions, code is parallelized, 
		#	chose the number of cores wants to use to parallelisation (good could be half of total numCores of computer)
		#	if no number of cores provided, default is 2
		# verbose: TRUE/FALSE - choose if print progress of analysis, default to TRUE

  library(RcppAlgos)
  library(data.table)
  library(parallel)
  library(dplyr)
  if(verbose==T){library(tictoc)}

  #print(paste0("Parameters used: Min number of vert./region:", minvert, " - Continuous fitting:", cont, " - Exhaustive search:", exhaus))

  Yvar <<- Yvar; Xvar <<- Xvar	# !! TO FIX !!  Temporary fix to be able to run without parallel (assigns Xvar and Yvar as global variables) - TO BE FIXED!!!!

  noBP <- noregions-1		# Get max number of breakpoints for given number of regions
  Yvar <- as.matrix(Yvar)

  rownames(Yvar) <- Xvar
  names(Xvar) <- 1:length(Xvar)


  noverts <- nrow(Yvar)
  noPC <- ncol(Yvar)
  if(noregions==1){colhead <- c("regions", "breakpoint1", "sumRSS", paste("RSS", 1:noPC, sep="."))
  } else { colhead <- c("regions", paste("breakpoint", 1:noBP, sep=""), "sumRSS", paste("RSS", 1:noPC, sep="."))
  }


  a <- Xvar[minvert:(length(Xvar)-minvert)]



# 1 Region (0 breakpoint - bp)

  if(verbose==T){tic()}

  nbp <- 0		# Number of breakpoints for 1 region
  lines <- lm(Yvar ~ Xvar)
  RSS <- sum(lines$residuals^2)
  if (noPC > 1){
	#rsq <- c()
	#for (i in 1:noPC){
	#  rsq[i] <- summary(lines)[[i]]$adj.r.squared
	#}
	rsq <- colSums(lines$residuals^2)
  } else { 
	#rsq <- summary(lines)$adj.r.squared
	rsq <- RSS
  }

  if(noregions==1){regionsAG <- data.table(t(c((nbp+1), 0, RSS, rsq)))	# Create data frame with results
  } else {regionsAG <- data.table(t(c((nbp+1), rep(0, (noBP-nbp)), RSS, rsq)))	# Create data frame with results
  }
  rownames(regionsAG) <- NULL; colnames(regionsAG) <- colhead
  if(noregions==1){
    stats <- data.frame(t(c(nbp+1, rep(1, 3), "Exhaustive", "All",NA)))
    colnames(stats) <- colstats <- c("Nregions", "Nmodel_possible", "Nmodel_tested", "Nmodel_saved", "Comp_method", "Saving_method", "Best_BPs1")
  } else {
    stats <- data.frame(t(c(nbp+1, rep(1, 3), "Exhaustive", "All",rep(NA,noBP))))
    colnames(stats) <- colstats <- c("Nregions", "Nmodel_possible", "Nmodel_tested", "Nmodel_saved", "Comp_method", "Saving_method",paste0("Best_BPs", 1:noBP))
  }

  if(verbose==T){print("1 region fitted: 1 model tested"); toc()}

  if (noregions == 1) {
     return(list(results=regionsAG, stats=stats))
  }



# 2 Regions (1 bp)

  if(verbose==T){tic()}

 # Get all possible combinations for 1 breakpoint:
  nbp <- 1
  goodcomb <- a

 # Run fitting on good combinations:
  res <- lapply(goodcomb, function(x){
	bp1 <- x		# First breakpoint

	if(cont==F){
	   lines<-lm(Yvar ~ Xvar*(Xvar<=bp1)
			  +Xvar*(Xvar>bp1))
	} else {
	   lines<-lm(Yvar ~ Xvar + pmax(0,Xvar-bp1))
	}

	RSS <- sum(lines$residuals^2)
	if (noPC > 1){
	   #rsq <- c()
	   #for (i in 1:noPC){
	   #  rsq[i] <- summary(lines)[[i]]$adj.r.squared
	   #}
	   rsq <- colSums(lines$residuals^2)
	} else { 
	   #rsq <- summary(lines)$adj.r.squared
	   rsq <- RSS
	}	
      c((nbp+1), x, rep(0, (noBP-nbp)), RSS, rsq)
  })
  res <- data.table(t(as.data.frame(res)))
  colnames(res) <- colhead; rownames(res) <- NULL

  bpkeep <- goodcomb
  bpkeep <- paste(bpkeep, collapse="|")

  if(exhaus==F){
    stat <- c(nbp+1, rep(nrow(res),3), "Exhaustive", "All",bpkeep,rep(NA,noBP-nbp))
  } else {
    stat <- c(nbp+1, rep(nrow(res),3), "Exhaustive", "All",rep(NA,noBP))
  }
  stats <- rbind(stats,stat)

  regionsAG <- rbind(regionsAG, res)

  if(verbose==T){print(paste0("2 regions fitted: ",nrow(res)," models tested")); toc()}

  if (noregions == 2) {
     return(list(results=regionsAG, stats=stats))
  }


# 3 and more Regions (2+ bp)

  for (j in 2:noBP){

	if(verbose==T){tic()}
	nbp <- j

	# Create function to filter combinations (to respect minimum number of vertebrae / region):
	f <- c()
	for(k in 1:(nbp-1)){
		f[k] <- paste0("x[",k+1,"]>=(x[",k,"]+minvert)")	# repeat the condition for as many bps as necessary
	}
	form <- paste(f, collapse=" && ")		# collapse everything as a string
	form <-parse(text=form)				# convert string to expression
	args <- alist(x=)					# define arguments of the expression
	fcomb <- as.function(c(args, form[[1]]))		# convert expression to function

	fregions <- function(x, Xvar, Yvar){
	  BPs <- x
	    if(cont==F){								# Discontinuous fit
		beg <- paste("Yvar ~ Xvar*(Xvar<=",BPs[1],")", sep="")
		end <- paste("+Xvar*(Xvar>",BPs[nbp],")", sep="")
		f <- c()
		for (i in 1:(nbp-1)){
		  f[i] <- paste("+ Xvar*(Xvar>",BPs[i]," & Xvar<=", BPs[i+1], ")", sep="")
		}
		form <- paste(beg, paste(f, collapse=""), end, sep="")

	    } else {								# Continuous fit
		f <- c()
		for (i in 1:nbp){
		  f[i] <- paste("+ pmax(0,Xvar-",BPs[i],")", sep="")
		}
		form <- paste("Yvar ~ Xvar", paste(f, collapse=""), sep="")
	    }
	  lines <- lm(form)
	  RSS <- sum(lines$residuals^2)
	  if (noPC > 1){
	    rsq <- colSums(lines$residuals^2)
	  } else { 
	    rsq <- RSS
	  }
	  if(nbp >= 7){		# erase fit and clean memory if more than 7 regions
	    rm(lines)
	    gc(verbose=F, full=F)
	  }
	  return(c((nbp+1), x, rep(0, (noBP-nbp)), RSS, rsq))
	}

	# Define BPs to include if non-exhaustive search:
	if(exhaus==F){
	  bpkeep <- unlist(stats[stats$Nregions==max(stats$Nregions),grep("Best_BPs", colnames(stats))])
	  bpkeep <- bpkeep[!is.na(bpkeep)]
	}

	ncombi <- comboCount(a, nbp)
	lim <- 0.5e+09			# max length for an R object to not exceed 2Gb in size (Memory safety net)
	div <- ceiling(lim/nbp)		# find nrow max for given number of cols (bps) to avoid exceeding 2Gb objects
	aa <- ceiling(ncombi/div)

	if(aa==1){
	  if(par==T){
	    cl <- makeCluster(numCores, type="PSOCK")
	    clusterEvalQ(cl,{library(regions)})
	    clusterExport(cl=cl, varlist=c('Xvar', 'Yvar', 'nbp', 'cont', 'noPC', 'noBP', 'minvert'), envir=environment())
	  }

	  goodcomb <- comboGeneral(as.numeric(names(a)), nbp)					# Generates all combinations by their position
	  goodcomb <- goodcomb[which((apply(goodcomb, 1, fcomb))==TRUE),]		# Keep only combs respecting the minimum number of vertebrae / region (minvert) using their name to account for unsampled vertebrae
	  if(is.null(nrow(goodcomb)) & !is.null(length(goodcomb))){goodcomb <- t(data.frame(goodcomb)); rownames(goodcomb)<-NULL} 	# Convert goodcomb to data.frame if only 1 goodcomb
	  goodcomb <- t(apply(goodcomb,1,function(x){a[names(a) %in% x]}))			# Convert position of bp to keep to actual value of bp
	  gc()

	  nmodel_possible <- nrow(goodcomb)

	  if(exhaus==F){ for(m in 1:length(bpkeep)){						# Keep only probable combinations (non-exhaustive search)
	    goodcomb <- as.matrix(as.data.frame(goodcomb) %>% filter_all(any_vars(grepl(bpkeep[m], .))))
	  }}
	  if(nrow(goodcomb)==0) next
	  if(par==T){
	    res <- data.table(t(parApply(cl, goodcomb, 1, fregions, Xvar=Xvar, Yvar=Yvar)))	# Fit all regressions with parallel computing
	  } else {
	    res <- data.table(t(apply(goodcomb, 1, fregions, Xvar=Xvar, Yvar=Yvar)))	# Fit all regressions without parallel
	  }
	  gc()
	  colnames(res) <- colnames(regionsAG)
	  nmodel_tested <- nrow(res)

	  if(exhaus==F){		# If non-exhaustive search, keep only models with sumRSS >= min(sumRSS)+(0.5*sd(sumRSS))
	    cutoff <- min(res$sumRSS)+(sd(res$sumRSS)/2)
	    if(!is.na(cutoff)){ res <- res[which(res$sumRSS <= cutoff),]}	# Subsample res only if more than 1 option
	  }

	  rm(goodcomb)


	  if (par==T){stopCluster(cl)}
	  gc()


	} else {		# Subset combinations to test in several smaller vectors (limit memory usage), results stored in tmp folder before merging all of them 
	  dir.create("tmp")
	  if(par==T){
	    cl <- makeCluster(numCores, type="PSOCK")
	    clusterEvalQ(cl,{library(regions)})
	    clusterExport(cl=cl, varlist=c('Xvar', 'Yvar', 'nbp', 'cont', 'noPC', 'noBP', 'minvert'), envir=environment())
	  }

	  nmodel_possible_sub <- c()

	  for(k in 1:aa){
	    low_lim <- (k-1)*div+1
	    up_lim <- k*div
	    if(up_lim > ncombi){up_lim <- ncombi}

	    goodcomb_sub <- comboGeneral(as.numeric(names(a)), nbp, lower=low_lim, upper=up_lim)		# Generates all combinations
	      gc()
	    goodcomb_sub <- goodcomb_sub[which((apply(goodcomb_sub, 1, fcomb))==TRUE),]	# Keep only combs respecting the minimum number of vertebrae / region (minvert)
	      gc()
	    if(is.null(nrow(goodcomb_sub)) & !is.null(length(goodcomb_sub))){goodcomb_sub <- t(data.frame(goodcomb_sub)); rownames(goodcomb_sub)<-NULL} 	# Convert goodcomb to data.frame if only 1 goodcomb
	    goodcomb_sub <- t(apply(goodcomb_sub,1,function(x){a[names(a) %in% x]}))			# Convert position of bp to keep to actual value of bp

	    nmodel_possible_sub <- c(nmodel_possible_sub,nrow(goodcomb_sub))


	    if(exhaus==F){ for(m in 1:length(bpkeep)){							# Keep only probable combinations (non-exhaustive search)
	      goodcomb_sub <- as.matrix(as.data.frame(goodcomb_sub) %>% filter_all(any_vars(grepl(bpkeep[m], .))))
	    }}
	    if(nrow(goodcomb_sub)==0) next
	    if(par==T){
	      res <- data.table(t(parApply(cl, goodcomb_sub, 1, fregions, Xvar=Xvar, Yvar=Yvar)))	# Fit all regressions with parallel computing
	    } else {
	      res <- data.table(t(apply(goodcomb_sub, 1, fregions, Xvar=Xvar, Yvar=Yvar)))	# Fit all regressions without parallel
	    }
	    colnames(res) <- colnames(regionsAG)

	    file <- paste0("tmp/res_part",k,".csv")
	    fwrite(res, file)
	    rm(res)
	    rm(goodcomb_sub)
	    gc()	
	  }

	  nmodel_possible <- sum(nmodel_possible_sub)

	  if(par==T){stopCluster(cl)}

	  # Import results of all smaller vectors tested:
	  if(exhaus==F){					# If non-exhaustive search, keep only models with sumRSS >= min(sumRSS)+(0.5*sd(sumRSS))
	    res_filenames <- list.files(path="tmp/", pattern=".csv")
	    res_RSS <- lapply(res_filenames, function(x){fread(paste0("tmp/", x), select="sumRSS")})
	    res_RSS <- unlist(res_RSS)
	    nmodel_tested <- length(res_RSS)

	    if(length(res_RSS) < 10^8){		# function sd doesn't work on very long vectors
	       cutoff <- min(res_RSS)+(sd(res_RSS)/2)
	    } else {
	       minim <- min(res_RSS)
	       var_all <- function(x) {       # Create function for population variance
	           (sum((x - mean(x))^2))/(length(x)-1)
	        }
	        sdd <- var_all(res_RSS)
	        sd2 <- sqrt(sdd)/2
	        cutoff <- minim+sd2
	    }

	    rm(res_RSS)
	    gc()

	    res <- bind_rows(lapply(res_filenames, function(x){
			y <- fread(paste0("tmp/", x))
			y <- y[which(y$sumRSS <= cutoff),]
			return(y)
	    }))

	  } else {
	    res <- bind_rows(lapply(res_filenames, function(x){fread(paste0("tmp/", x))}))	# note: not tested on extremely large datasets, might not be supported?
	  }

	  unlink("tmp", recursive=TRUE)
	  gc()

	}

 	nmodel_saved <- nrow(res)

	# Create string with info on best BPs to use for next region:	# BPs kept are all  > cutoff (min+0.5*sd) and BPs of best models +/- 3 vertebrae 
	if(exhaus==F){
	  bestBPs <- res[which(res$sumRSS==min(res$sumRSS)),2:(nbp+1)]

	  bbp <- as.numeric(names(a[a %in% bestBPs]))
	  bbp <- as.list(data.frame(rbind(bbp, bbp-3,bbp-2,bbp-1,bbp+1,bbp+2,bbp+3)))
	  bestBPs <- lapply(bbp,function(x){a[names(a) %in% x]})			# Convert position of bp to keep to actual value of bp
	  bpkeep <- as.list(res[,2:(nbp+1)])
	  for(z in 1:length(bpkeep)){bpkeep[[z]]<-c(bestBPs[[z]],bpkeep[[z]])}
	  bpkeep <- lapply(bpkeep, function(x)unique(x))

	  if(class(bpkeep)[1]=="list"){
	    bpkeep <- lapply(bpkeep, function(x){x[x>0]})			# Remove potential bp <= 0
	    bpkeep <- unlist(lapply(bpkeep, function(x){x<-x[order(x)];paste(x, collapse="|")}))
	  } else {
	    bpkeep <- apply(bpkeep, 2, function(x){x[x>0]})			# Remove potential bp <= 0
	    bpkeep <- unlist(apply(bpkeep, 2, function(x){x<-x[order(x)];paste(x, collapse="|")}))
	  }
	  stat <- c(nbp+1, nmodel_possible, nmodel_tested, nmodel_saved, "Non-exhaus", "SD/2", bpkeep, rep(NA,noBP-nbp))
	} else {
	  stat <- c(nbp+1, nmodel_possible, nmodel_tested, nmodel_saved, "Exhaustive", "All", rep(NA,noBP))
	}

	stats <- rbind(stats,stat)

	regionsAG <- rbind(regionsAG,res)
	rm(res)

	todrop <- which(regionsAG$regions==0)		# drop first line with 0 if not prevreg file provided
	if(length(todrop)!=0){regionsAG <- regionsAG[-todrop,]}
	
	gc(full=T)

	if(verbose==T){print(paste0((nbp+1)," regions fitted: ",nmodel_tested," models tested")); toc()}

  }
  return(list(results=regionsAG, stats=stats))
}







###########
## Function: calculate additional regions ##
###########
	# With data.table instead of data.frames (should require less RAM)
	# With tmp file


## NOTE: add option to input stats from previous run and merge it with new results
# Check if function works if wants to add more than 1 additional region at a time


addregionsAG <- function(Xvar, Yvar, prevreg=NULL, prevnoregions=NULL, bpkeep=NULL, noregions, minvert=3, cont=T, exhaus=T, par=T, numCores=2, verbose=T){
	# prevreg (optional): table with results from previous analysis with less regions
	# prevnoregions: number of regions previously fitted
	# bpkeep: a vector of character from output of previous region search with breakpoints to test if doesn't want exhaustive search, similar structure as output of calcregions but only with breakpoint columns (cols= bp1, bp2, ...; rows= each bp)

  library(RcppAlgos)
  library(data.table)
  library(parallel)
  library(dplyr)
  library(tibble)
  if(verbose==T){library(tictoc)}

  Yvar <<- Yvar; Xvar <<- Xvar	# !! TO FIX !!  Temporary fix to be able to run without parallel (assigns Xvar and Yvar as global variables) - TO BE FIXED!!!!

  Yvar <- as.matrix(Yvar)
  noverts <- nrow(Yvar)
  noPC <- ncol(Yvar)
  noBP <- noregions-1		# Get max number of breakpoints

  rownames(Yvar) <- Xvar
  names(Xvar) <- 1:length(Xvar)



  if(!is.null(prevreg)){
    prevBP <- grep('breakpoint', colnames(prevreg), value=TRUE)
    prevnoBP <- length(prevBP)		# Get max number of breakpoints from previous analysis
    BPstoadd <- (prevnoBP+1):noBP		# BPs to add
    newcols <- as.data.frame(matrix(0,ncol=length(BPstoadd), nrow=nrow(prevreg)))	# Create matrix of 0 for new BPs
    colnames(newcols) <- paste0("breakpoint", BPstoadd)		# Set colnames for new BPs
    regionsAG <- add_column(prevreg, newcols, .after = prevBP[length(prevBP)])	# Insert new cols into table of previous analyses
  } else if(!is.null(prevnoregions)){
    prevnoBP <- prevnoregions-1
    BPstoadd <- (prevnoBP+1):noBP		# BPs to add
    regionsAG <- data.table(as.data.frame(matrix(0, ncol=(1+noBP+1+noPC), nrow=1)))
    colhead <- c("regions", paste("breakpoint", 1:noBP, sep=""), "sumRSS", paste("RSS", 1:noPC, sep="."))
    colnames(regionsAG) <- colhead
  } else {stop("Previous number of regions not provided")}


  # Ensure that bpkeep is provided is non-exhaustive search:
  if(exhaus==F & is.null(bpkeep)){stop("Breakpoints need to be provided for non-exhaustive search")}

  a <- Xvar[minvert:(length(Xvar)-minvert)]


  stats <- c()

  for (j in (prevnoBP+1):noBP){

	if(verbose==T){tic()}
	nbp <- j

	# Create function to filter combinations:
	f <- c()
	for(k in 1:(nbp-1)){
		f[k] <- paste0("x[",k+1,"]>=(x[",k,"]+minvert)")	# repeat the condition for as many bps as necessary
	}
	form <- paste(f, collapse=" && ")		# collapse everything as a string
	form <-parse(text=form)				# convert string to expression
	args <- alist(x=)					# define arguments of the expression
	#fcomb <- as.function(c(args, form[[1]]),env=parent.frame())		# convert expression to function
	fcomb <- as.function(c(args, form[[1]]))		# convert expression to function


	fregions <- function(x, Xvar, Yvar){
	  BPs <- x
	    if(cont==F){								# Discontinuous fit
		beg <- paste("Yvar ~ Xvar*(Xvar<=",BPs[1],")", sep="")
		end <- paste("+Xvar*(Xvar>",BPs[nbp],")", sep="")
		f <- c()
		for (i in 1:(nbp-1)){
		  f[i] <- paste("+ Xvar*(Xvar>",BPs[i]," & Xvar<=", BPs[i+1], ")", sep="")
		}
		form <- paste(beg, paste(f, collapse=""), end, sep="")

	    } else {								# Continuous fit
		f <- c()
		for (i in 1:nbp){
		  f[i] <- paste("+ pmax(0,Xvar-",BPs[i],")", sep="")
		}
		form <- paste("Yvar ~ Xvar", paste(f, collapse=""), sep="")
	    }
	  lines <- lm(form)

	  RSS <- sum(lines$residuals^2)
	  if (noPC > 1){
	    rsq <- colSums(lines$residuals^2)
	  } else { 
	    rsq <- RSS
	  }
	  if(nbp >= 7){		# erase fit and clean memory if 8 regions or more
	    rm(lines)
	    gc(verbose=F, full=F)
	  }
	  return(c((nbp+1), x, rep(0, (noBP-nbp)), RSS, rsq))
	}


	ncombi <- comboCount(a, nbp)
	lim <- 0.5e+09			# max length for an R object to not exceed 2Gb in size (Memory safety net)
	div <- ceiling(lim/nbp)		# find nrow max for given number of cols (bps) to avoid exceeding 2Gb objects
	aa <- ceiling(ncombi/div)

	if(aa==1){
	  if(par==T){
	    cl <- makeCluster(numCores, type="PSOCK", setup_strategy = "sequential")
	    clusterEvalQ(cl,{library(regions)})
	    clusterExport(cl=cl, varlist=c('Xvar', 'Yvar', 'nbp', 'cont', 'noPC', 'noBP'), envir=environment())
	  }

	  goodcomb <- comboGeneral(as.numeric(names(a)), nbp)			# Generates all combinations by their position
	  goodcomb <- goodcomb[which((apply(goodcomb, 1, fcomb))==TRUE),]	# Keep only combs respecting the minimum number of vertebrae / region (minvert) using their name to account for unsampled vertebrae
	  if(is.null(nrow(goodcomb)) & !is.null(length(goodcomb))){goodcomb <- t(data.frame(goodcomb)); rownames(goodcomb)<-NULL} 	# Convert goodcomb to data.frame if only 1 goodcomb
	  goodcomb <- t(apply(goodcomb,1,function(x){a[names(a) %in% x]}))			# Convert position of bp to keep to actual value of bp
	  gc()

	  nmodel_possible <- nrow(goodcomb)


	  if(exhaus==F){ for(m in 1:length(bpkeep)){						# Keep only probable combinations (non-exhaustive search)
	    goodcomb <- as.matrix(as.data.frame(goodcomb) %>% filter_all(any_vars(grepl(bpkeep[m], .))))
	  }}
	  if(nrow(goodcomb)==0) next
	  if(par==T){
	    res <- data.table(t(parApply(cl, goodcomb, 1, fregions, Xvar=Xvar, Yvar=Yvar)))
	  } else {
	    res <- data.table(t(apply(goodcomb, 1, fregions, Xvar=Xvar, Yvar=Yvar)))	# Fit all regressions without parallel
	  }
	  gc()
	  colnames(res) <- colnames(regionsAG)
	  nmodel_tested <- nrow(res)

	  if(exhaus==F){		# If non-exhaustive search, keep only models with sumRSS >= min(sumRSS)+(0.5*sd(sumRSS))
	    cutoff <- min(res$sumRSS)+(sd(res$sumRSS)/2)
	    res <- res[which(res$sumRSS <= cutoff),]
	  }

	  rm(goodcomb)

	  if (par==T){stopCluster(cl)}
	  gc()


	} else {
	  dir.create("tmp")

	  if(par==T){
	    cl <- makeCluster(numCores, type="PSOCK", setup_strategy = "sequential")
	    clusterEvalQ(cl,{library(regions)})
	    clusterExport(cl=cl, varlist=c('Xvar', 'Yvar', 'nbp', 'cont', 'noPC', 'noBP'), envir=environment())
	  }


	  nmodel_possible_sub <- c()
	  for(k in 1:aa){
	    low_lim <- (k-1)*div+1
	    up_lim <- k*div
	    if(up_lim > ncombi){up_lim <- ncombi}

	    goodcomb_sub <- comboGeneral(as.numeric(names(a)), nbp, lower=low_lim, upper=up_lim)		# Generates all combinations by their position
	      gc()
	    goodcomb_sub <- goodcomb_sub[which((apply(goodcomb_sub, 1, fcomb))==TRUE),]	# Keep only combs respecting the minimum number of vertebrae / region (minvert)
	      gc()
	    if(is.null(nrow(goodcomb_sub)) & !is.null(length(goodcomb_sub))){goodcomb_sub <- t(data.frame(goodcomb_sub)); rownames(goodcomb_sub)<-NULL} 	# Convert goodcomb to data.frame if only 1 goodcomb
	    goodcomb_sub <- t(apply(goodcomb_sub,1,function(x){a[names(a) %in% x]}))			# Convert position of bp to keep to actual value of bp

	    nmodel_possible_sub <- c(nmodel_possible_sub,nrow(goodcomb_sub))


	    if(exhaus==F){ for(m in 1:length(bpkeep)){
	      goodcomb_sub <- as.matrix(as.data.frame(goodcomb_sub) %>% filter_all(any_vars(grepl(bpkeep[m], .))))
	    }}
	    if(nrow(goodcomb_sub)==0) next

	    if(par==T){
	      res <- data.table(t(parApply(cl, goodcomb_sub, 1, fregions, Xvar=Xvar, Yvar=Yvar)))
	    } else {
	      res <- data.table(t(apply(goodcomb_sub, 1, fregions, Xvar=Xvar, Yvar=Yvar)))	# Fit all regressions without parallel
	    }
	    colnames(res) <- colnames(regionsAG)

	    file <- paste0("tmp/res_part",k,".csv")
	    fwrite(res, file)
	    rm(res)
	    rm(goodcomb_sub)
	    gc()	
	  }

	  nmodel_possible <- sum(nmodel_possible_sub)


	  if(par==T){stopCluster(cl)}

	  # Import results of all smaller vectors tested:
	  if(exhaus==F){					# If non-exhaustive search, keep only models with sumRSS >= min(sumRSS)+(0.5*sd(sumRSS))
	    res_filenames <- list.files(path="tmp/", pattern=".csv")
	    res_RSS <- lapply(res_filenames, function(x){fread(paste0("tmp/", x), select="sumRSS")})
	    res_RSS <- unlist(res_RSS)
	    nmodel_tested <- length(res_RSS)

	    if(length(res_RSS) < 10^8){		# function sd doesn't work on very long vectors
	      cutoff <- min(res_RSS)+(sd(res_RSS)/2)
	    } else {
	      minim <- min(res_RSS)
	      var_all <- function(x) {       # Create function for population variance
	           (sum((x - mean(x))^2))/(length(x)-1)
	      }
	      sdd <- var_all(res_RSS)
	      sd2 <- sqrt(sdd)/2
	      cutoff <- minim+sd2
	    }

	    rm(res_RSS)
	    gc()

	    res <- bind_rows(lapply(res_filenames, function(x){
			y <- fread(paste0("tmp/", x))
			y <- y[which(y$sumRSS <= cutoff),]
			return(y)
	    }))

	  } else {
	    res <- bind_rows(lapply(res_filenames, function(x){fread(paste0("tmp/", x))}))	# note: not tested on extremely large datasets, might not be supported?
	  }

	  unlink("tmp", recursive=TRUE)
	  gc()

	}

	nmodel_saved <- nrow(res)

	# Create string with info on best BPs to use for next region:	# BPs kept are all  > cutoff (min+0.5*sd) and BPs of best models +/- 3 vertebrae 
	if(exhaus==F){
	  bestBPs <- res[which(res$sumRSS==min(res$sumRSS)),2:(nbp+1)]
	  bbp <- as.numeric(names(a[a %in% bestBPs]))
	  bbp <- as.list(data.frame(rbind(bbp, bbp-3,bbp-2,bbp-1,bbp+1,bbp+2,bbp+3)))
	  bestBPs <- lapply(bbp,function(x){a[names(a) %in% x]})			# Convert position of bp to keep to actual value of bp
	  bpkeep <- as.list(res[,2:(nbp+1)])
	  for(z in 1:length(bpkeep)){bpkeep[[z]]<-c(bestBPs[[z]],bpkeep[[z]])}
	  bpkeep <- lapply(bpkeep, function(x)unique(x))

	  if(class(bpkeep)[1]=="list"){
	    bpkeep <- lapply(bpkeep, function(x){x[x>0]})			# Remove potential bp <= 0
	    bpkeep <- unlist(lapply(bpkeep, function(x){x<-x[order(x)];paste(x, collapse="|")}))
	  } else {
	    bpkeep <- apply(bpkeep, 2, function(x){x[x>0]})			# Remove potential bp <= 0
	    bpkeep <- unlist(apply(bpkeep, 2, function(x){x<-x[order(x)];paste(x, collapse="|")}))
	  }
	  stat <- data.frame(t(c(nbp+1, nmodel_possible, nmodel_tested, nmodel_saved, "Non-exhaus", "SD/2", bpkeep, rep(NA,noBP-nbp))))
	} else {
	  stat <- data.frame(t(c(nbp+1, nmodel_possible, nmodel_tested, nmodel_saved, "Exhaustive", "All", rep(NA,noBP))))
	}
	rownames(stat) <- NULL
	colnames(stat) <- c("Nregions", "Nmodel_possible", "Nmodel_tested", "Nmodel_saved", "Comp_method", "Saving_method", paste0("Best_BPs", 1:noBP))
	stats <- rbind(stats,stat)

	regionsAG <- rbind(regionsAG,res)
	rm(res)

	todrop <- which(regionsAG$regions==0)		# drop first line with 0 if not prevreg file provided
	if(length(todrop)!=0){regionsAG <- regionsAG[-todrop,]}

	gc(full=T)

	if(verbose==T){print(paste0((nbp+1)," regions fitted: ",nmodel_tested," models tested")); toc()}
  }
  return(list(results=regionsAG, stats=stats))
}





###########
### Function to calculate a single predefined model ###
###########
	# can be used for model testing/comparison


calcmodel <- function(Xvar, data, BPs, cont=T){
	# Xvar = vector with vertebral number included in the analysis
	# data = dataframe with pco scores (columns = pco axes, rows = vertebrae)
	# BPs vector containing breakpoint(s) position
	# cont: TRUE/FALSE - choose if want to use continuous (TRUE) or discontinuous segmentation fitting, default to TRUE

  Yvar <<- as.matrix(data); Xvar <<- Xvar
  rownames(Yvar) <- Xvar
  names(Xvar) <- 1:length(Xvar)

  nBPs <- length(BPs)
  noregions <- nBPs+1
  noPC <- ncol(Yvar)


 # Fit model and calculate multivariate and univariate RSS:

  if(noregions==1){colhead <- c("regions", "breakpoint1", "sumRSS", paste("RSS", 1:noPC, sep="."))
  } else { colhead <- c("regions", paste("breakpoint", 1:nBPs, sep=""), "sumRSS", paste("RSS", 1:noPC, sep="."))
  }

  if(nBPs ==0){
    form <- "Yvar ~ Xvar"

  } else {
    if(cont==F){	# Discontinuous fit
	beg <- paste("Yvar ~ Xvar*(Xvar<=",BPs[1],")", sep="")
	end <- paste("+Xvar*(Xvar>",BPs[nBPs],")", sep="")

	if(nBPs==1){
	   form <- paste(beg, end, sep="")
	} else if(nBPs >1){
	   f <- c()
	   for (i in 1:(nBPs-1)){
	      f[i] <- paste("+ Xvar*(Xvar>",BPs[i]," & Xvar<=", BPs[i+1], ")", sep="")
	   }
	   form <- paste(beg, paste(f, collapse=""), end, sep="")
	}
    } else {		# Continuous fit
	f <- c()
	for (i in 1:nBPs){
	   f[i] <- paste("+ pmax(0,Xvar-",BPs[i],")", sep="")
	}
	form <- paste("Yvar ~ Xvar", paste(f, collapse=""), sep="")
    }
  }

  lines <- lm(form)

  RSS <- sum(lines$residuals^2)
  if (noPC > 1){
   rsq <- colSums(lines$residuals^2)
  } else { 
   rsq <- RSS
  }

  regions <- data.frame(t(c((nBPs+1), BPs, RSS, rsq)  ))
  colnames(regions) <- colhead

 # Calculate AICc & BIC:
  stat <- AICcalcAG(regions$sumRSS,nPC=noPC,nvert=length(Xvar),noregions=noregions,cont=cont)


  return(list(results=regions,AICc=stat[1],BIC=stat[2]))
}





###########
### Function model select ###
###########


modelselectAG<-function(regiondata, PCOs.no=NULL){
	# regiondata: $results value from calcregionsAG output
	# PCOs.no (optional): if null, will select models based on analysis on all PCOs included in calregionsAG
		#	if specified (either as a single numeric value, i.e., 1 or as a range, i.e., 1:3) will run the analysis only on PCOs specified in the range (will exclude RSS values from PCOs not included in the range)

  noregions <- max(regiondata$regions)
  regiondata<-as.data.frame(regiondata)
  models<-numeric()
  for (i in 1:noregions){
    allmodels<-subset(regiondata, regiondata$regions==i)	#select only models with correct region no
    if(is.null(PCOs.no)){
      best<-allmodels[which(allmodels$sumRSS==min(allmodels$sumRSS)),]	#select the lowest RSS
    } else {
      rss.keep <- data.frame(allmodels[,grep("RSS.",colnames(allmodels))[PCOs.no]])
      if(ncol(rss.keep)==1){colnames(rss.keep) <- paste0("RSS.",PCOs.no)}
      sumRSS <- rowSums(rss.keep)
      allmodels <- cbind(allmodels[,grep("breakpoint|regions", colnames(allmodels))], sumRSS,rss.keep)
      best<-allmodels[which(allmodels$sumRSS==min(allmodels$sumRSS)),]	#select the lowest RSS
    }
    rm(allmodels)
    gc()
    models<-rbind(models, best)#fill in model table
  }
  return(models)
}



###########
### Function model support ###
###########


model_supportAG<-function(models, nvert, cont){
  nPC <- length(grep("RSS." ,colnames(models)))
  AICc=numeric()
  BIC=numeric()
  for (i in 1:nrow(models)){
    probs<-AICcalcAG(models$sumRSS[i],nPC, nvert, models$regions[i], cont) #Calculate AIC score
		# Only change compared to KJ code: use AICcalcAG function (which accounts for continuous or discont. fitting
    AICc<-rbind(AICc,probs[1])
    BIC<-rbind(BIC,probs[2])
  }

 # AICc:
  AICmin<-min(AICc)	#AIC of best model
  deltaAIC<-sapply(AICc, function(x) x-AICmin)	#Calculate AIC difference
  model_lik<-sapply(deltaAIC, function(x) exp(-0.5*x))	#Likelihood of the model
  tot_lik<-sum(model_lik)
  Ak_weight<-sapply(model_lik, function(x) x/tot_lik)	#Akaike Weights
  AIC_models<-cbind(models, AICc, deltaAIC, model_lik, Ak_weight)
  AIC_models<-AIC_models[order(AIC_models$AICc),]#Sort so best at top

  weight_region<-AIC_models$regions*AIC_models$Ak_weight
  Regions_score<-sum(weight_region)

 # BIC:
  BICmin<-min(BIC)	#AIC of best model
  deltaBIC<-sapply(BIC, function(x) x-BICmin)	#Calculate AIC difference
  model_lik<-sapply(deltaBIC, function(x) exp(-0.5*x))	#Likelihood of the model
  tot_lik<-sum(model_lik)
  BIC_weight<-sapply(model_lik, function(x) x/tot_lik)	#BIC Weights
  BIC_models<-cbind(models, BIC, deltaBIC, model_lik, BIC_weight)
  BIC_models<-BIC_models[order(BIC_models$BIC),]	#Sort so best at top

  weight_regionB<-BIC_models$regions*BIC_models$BIC_weight
  Regions_scoreB<-sum(weight_regionB)

  return(list(Model_support=AIC_models,Region_score=Regions_score,Model_support_BIC=BIC_models,Region_score_BIC=Regions_scoreB))
}



###########
### Function AICc ###
###########



AICcalcAG<-function(RSS, nPC, nvert, noregions, cont){
  n=nPC*nvert 	# No of variables used
  var=RSS/n 	# Variance calculated ML way
  if(cont==F){	# Calculating the number of parameters (k) being fitted in the model
    k=(2*noregions*nPC)+(noregions-1) 	# For discontinuous fit, a slope and an intercept (2 params) are evaluated for region and for each PCO, + a noregions-1 number of breakpoint is estimated
  } else {
    k=(noregions*nPC)+nPC+(noregions-1)	# For a continuous fit, a slope is estimated for each region and each PCO, only the intercept of the fist segment is estimated for each PCO since other intercepts are forced by the previous slope and breakpoint position, + a noregions-1 number of breakpoint is estimated
  }
  if(n<(k+2)) stop('ratio of variables to parameters too small. Reduce number of regions or increase variables')
  AIC=n*log(var)+(2*k)
  corr=(2*k*(k+1))/(n-k-1)	# Correct for number of parameters and small sample
  AICc=AIC+corr 	# Calculate AICc
  BIC=n*log(var)+(log(n)*k)
  return(c(AICc, BIC))
}





################################
## Function simulate regions ##
################################

#' @param nvert Number of verts
#' @param noregions Number of regions
#' @param ersd Amount of error (sd)
#' @param plot Should the data be plot?
#' @param nvar number of variables



sim_region<-function(nvert, noregions, ersd=0.075, plot=T, nvar=1, minvert=3, cont, sl.lims=c(-2,2), sl.dif=NULL){
		# sl.lims = minimal and maximal values of simulated slopes (default to -2 and 2)

  Xvar<-c(1:nvert)
  br<-sort(sample(seq(from=minvert, to=(nvert-minvert), by=minvert),(noregions-1)))
  xs<-split(Xvar, cut(Xvar,c(-Inf,br,Inf)))

  slope <- int <- matrix(nrow=noregions, ncol=nvar)
  if(cont==F){			# if no continuous fitting, just generates random slopes and intercepts
    for(i in 1:nvar){
	good <- FALSE
	while(!good){
	  slope[,i] <- runif(noregions, min=sl.lims[1],max=sl.lims[2])
	  if(is.null(sl.dif)){
	    good <- slope
	  } else {
	    good <- (!FALSE %in% (abs(diff(slope[,i], lag=1)) > sl.dif))
	  }
	}
	int[,i] <- runif(noregions, 0,1)
    }
  } else {			# if continuous fitting needs to match intercepts of regions 2 to n
    for(i in 1:nvar){
	good <- FALSE
	while(!good){
	  slope[,i] <- runif(noregions, min=sl.lims[1],max=sl.lims[2])
	  if(is.null(sl.dif)){
	    good <- slope
	  } else {
	    good <- (!FALSE %in% (abs(diff(slope[,i], lag=1)) > sl.dif))
	  }
	}
	for(j in 1:noregions){
	  if(j==1){
	    int[j,i] <- runif(1, 0,1)		# generates random intercept for fist region
	  } else {
		int[j,i] <- int[j-1,i] + (br[j-1]+0.5) * (slope[j-1,i]-slope[j,i])
	  }
	}
    }
  }

  y<-matrix(NA, nrow=nvert, ncol=nvar)
  for(j in 1:noregions){
    for(a in 1:nvar){
      x<-xs[[j]]
      e<-rnorm(length(x), mean=0, sd=ersd)
      c<-(noregions*j)+a-noregions
      yi<-(slope[j,a]*x)+int[j,a]+e
      y[x,a]<-yi
    }
  }

  if(isTRUE(plot)){
    p<-plot(Xvar,y[,1])
    abline(v=(br+0.5))
    title(main="Simulated segmented regression")
  }

  return(list(Xvar=Xvar, y=y, br=br, slope=slope, int=int))

}




###########
### Function for multivariate R² calculation ###
###########


multivarrsqAG <- function(Xvar, Yvar, cont, bps=NULL, modelsupport=NULL, model=NULL){
	# cont: TRUE/FALSE
	# BPs (optional): vector containing BP position
	# modelsupport (optional): object $Model_support or $Model_support_BIC from model_supportAG output
	# model (optional): number corresponding to the row of modelsupport to calculate (1 = best model, 2 = 2nd best model, etc)
	# note: BPs & modelsupport are optional but at least one of them must be provided! If both are provided, modelsupport & model will be ignored

  if(is.null(bps) & is.null(modelsupport)){stop("bps or modelsupport argument must be provided")}

  if(!is.null(modelsupport)){
    if(is.null(model)){model <- 1}
    bestmodel <- modelsupport[model,]
    keep <- grep('breakpoint', colnames(bestmodel), value=TRUE)
    BPs <- unlist(bestmodel[,keep])
    nBPs <- length(BPs)
  }

  if(!is.null(bps)){
    BPs <- bps
    nBPs <- length(BPs)
  }


  totrsq<-data.frame(matrix(NA, nrow=ncol(Yvar), ncol=6, dimnames=list(c(1:ncol(Yvar)), c("rsq","adj.rsq", "SSres", "df", "SStot", "dfe"))))


 # Calculate R² on each PCO univariately:
  for(j in 1:ncol(Yvar)){
    pco <<- Yvar[,j]	# Not sure why if simply using <- instead of <<- the variable pco is not passed along for lm

    if(nBPs ==0){
      form <- "pco ~ Xvar"

    } else {
      if(cont==F){	# Discontinuous fit
	beg <- paste("pco ~ Xvar*(Xvar<=",BPs[1],")", sep="")
	end <- paste("+Xvar*(Xvar>",BPs[nBPs],")", sep="")

	if(nBPs==1){
	   form <- paste(beg, end, sep="")
	} else if(nBPs >1){
	   f <- c()
	   for (i in 1:(nBPs-1)){
	      f[i] <- paste("+ Xvar*(Xvar>",BPs[i]," & Xvar<=", BPs[i+1], ")", sep="")
	   }
	   form <- paste(beg, paste(f, collapse=""), end, sep="")
	}
     } else {		# Continuous fit
	f <- c()
	for (i in 1:nBPs){
	   f[i] <- paste("+ pmax(0,Xvar-",BPs[i],")", sep="")
	}
	form <- paste("pco ~ Xvar", paste(f, collapse=""), sep="")
      }
    }
    fit <- lm(form)
    aov <- anova(fit)
    totrsq[j,] <- c(summary(fit)$r.squared, summary(fit)$adj.r.squared, aov$`Sum Sq`[length(aov$`Sum Sq`)],
			sum(aov$Df),sum(aov$`Sum Sq`),aov$Df[length(aov$Df)])
  }
  rownames(totrsq) <- paste("PC",1:ncol(Yvar),sep=".")

 # Calculate R² for all PCOs multivariately:
  tot.rsq <- as.data.frame(t(colSums(totrsq)))
  ord.tot.rsq <- 1-(tot.rsq$SSres/tot.rsq$SStot)
  adj.tot.rsq <- 1-((tot.rsq$SSres/tot.rsq$dfe)/(tot.rsq$SStot/tot.rsq$df))

  return(list(univariate=totrsq[,1:2],multivariate=cbind(rsq=ord.tot.rsq, adj.rsq=adj.tot.rsq)))

}



###########
### Function for calculating weighted mean and SD of BP position for a given number of regions ###
###########



calcBPvar <- function(regiondata, nreg, nmodel, pct, nvert, cont){
	# regiondata = $results output from calcregionsAG
	# nreg = number of regions for which want to get weighted mean and SD (must be a number of regions that is included in regiondata object)
	# nmodel = total number of possible models for that number of regions and vertebrae (provided in the "Nmodel_possible" column of the $stats object from calcregionsAG output) 
	# pct = percentage of best models to keep from the original total number of possible models (nmodel)
	# nvert = number of vertebrae that were included in the calcregionsAG analysis
	# cont = TRUE/FALSE, should be the same value as the one provided for calcregionsAG analysis


  nmodel <- as.numeric(nmodel)
  ntop5 <- ceiling(nmodel/100*pct)	# Define number of models to keep to correspond to given pct of models wants to keep

  nPC <- length(grep("RSS.",colnames(regiondata)))	# Get number of PCs on which analysis was run

  dat <- regiondata[regiondata$regions==nreg,]	# Extract models corresponding to the given number of regions
  dat <- dat[order(dat$sumRSS),]	# Order by increasing sumRSS (from best to worst model)
  
  if(nrow(dat)>=ntop5){
    dat <- dat[1:ntop5,]
  } else {warning(paste0("Number of models provided lower than percentage requested.\n Weighted means and SD calculated on: ",
		(round(nrow(dat)/nmodel*100,2)),"% of total number of models."))}

 # Calculate probability (with AICc) and weight of each model:
  AICcs <- sapply(dat$sumRSS,FUN=AICcalcAG,nPC=nPC,nvert=nvert,noregions=nreg,cont=cont)[1,]	# Calculate AICc of each model
  AICmin <- min(AICcs)	# AICc of best model
  deltaAIC <- sapply(AICcs, function(x) x-AICmin)	# Calculate AIC difference
  model_lik <- sapply(deltaAIC, function(x) exp(-0.5*x))
  tot_lik <- sum(model_lik)
  Ak_weight <- sapply(model_lik, function(x) x/tot_lik)	# Akaike weights

 # Calculate weighted mean and weighted SD of each BP position using Akaike weights: (formulae from Symonds & Moussalli, Behav Ecol Sociobiol 2011)
  bps <- data.table(data.frame(dat)[,grep("breakpoint", colnames(dat))])
  wMean <- apply((bps*Ak_weight),2,sum)	# Weighted mean is the sum of bps values multiplied by the corresponding Akaike weight
  wSD <- c()
  for(i in 1:ncol(bps)){	# Working 1 bp at a time
    y <- bps[,..i]
    wmean <- wMean[i]
    sdy <- sapply(y, function(x){(x-wmean)^2})	# residuals calculated as difference from the weighted mean
    wsdy <- sdy*Ak_weight	# residuals are multiplied by AIC weight
    wSD <- c(wSD,sqrt(sum(wsdy)))	# weighted SD = squared root of sum of weighted residuals
  }
  means <- data.frame(rbind(wMean,wSD))

  dat <- cbind(dat,AICweight=Ak_weight,CumWeight=cumsum(Ak_weight))

  return(list(WeightedBp = means, BestModels=dat))
	# WeightedBP = data.frame with row1 = weighted means and row2 = weighted SD
	# BestModels = top pct % models ordered from best to worst with corresponding AIC weight and cumulative weight
}





###########
### Function for calculating PCO loadings ###
###########

pco.loadAG <- function(data, PCOscores){
	# data = data frame used for calculating PCOscore in formula svdPCO
	# PCOscores = data frame with scores from svdPCO

  vert.size<-apply(data,1,mean)
  data<-cbind(data, vert.size)
  load.pco <- list()
  for(i in 1:ncol(PCOscores)){
   #Figure out loadings on each PCO using correlation:
    PCOscore <- PCOscores[,i]
    load.pco[[i]]<-stats::cor(data[,1:ncol(data)],PCOscore, use="pairwise.complete.obs")
		# QUESTION FOR KATRINA: Original code started at 3rd column of data: "data[,3:ncol(data)] ==> WHY??????
  }
  load.pco <- do.call(cbind,load.pco)
  colnames(load.pco) <- paste("PCO",1:ncol(PCOscores),sep=".")
  return(loadings=load.pco)
}



###########
### Function for finding optimal number of PCO to maximizing region score ###
###########

PCOmaxAG <- function(regiondata, nvert, cont){

  noregions <- max(regiondata$regions)
  nvar <- length(grep("RSS." ,colnames(regiondata)))

  pco.no.test<-data.frame(matrix(NA,nrow=nvar,ncol=3, dimnames=list(c(1:nvar),c("PCO", "RS_AICc", "RS_BIC"))))

  for (i in 1:nvar){

   #Run for cumulative PCs
    models.cum <- modelselectAG(regiondata, PCOs.no=1:i)
    support.cum<-model_supportAG(models.cum,nvert,cont)

    pco.no.test[i,1]<-i
    pco.no.test[i,2]<-support.cum$Region_score
    pco.no.test[i,3]<-support.cum$Region_score_BIC
  }

  pco.max.AICc <- min(which(pco.no.test[,"RS_AICc"]==max(pco.no.test[,"RS_AICc"])))
  pco.max.BIC <- min(which(pco.no.test[,"RS_BIC"]==max(pco.no.test[,"RS_BIC"])))
  return(list(pco.max.AICc=pco.max.AICc, pco.max.BIC=pco.max.BIC, pco.dist=pco.no.test))
}








########################
########################
## Plotting functions ##
########################
########################




###########
### Function plot segmented rgeression ###
###########

plotsegregAG <- function(Xvar, pcono, data, modelsupport, model, BPs, cont, xlim, ylim, title, pch, col.pts, col.fit, col.BP, lwd.fit, lwd.BP) {
	# Xvar = vector with vertebral number included in the analysis
	# pcono = number of the pco want to plot (i.e., 1, 2, 3)
	# data = dataframe with pco scores (columns = pco axes, rows = vertebrae)
	# modelsupport (OPTIONAL) = object Model_support from model_support function output   /!\ modelsupport and BPs optional but at least one of them must be provided!

	# Added from KJ code:
		# model (OPTIONAL): number to select which model to plot (1 = best, 2 = 2nd best, etc)
		# cont: TRUE/FALSE depending if calculated regions using continuous fitting or not
		# BPs (OPTIONAL): if modelsupport and model are not provided, can provide directly vector with ordered breakpoints - if no BPs (ie, 1 region) use "NULL"
		# Plotting options: xlim, ylim, title, pch, col.pts, col.fit, col.BP, lwd.fit, lwd.BP
	# Modified from KJ code: 
		# changed modelsupport column subsetting to allow variable number of regions
		# changed code to adapt it to different number of regions
		# changes code to allow inputing BPs directly instead of output from model_support function

  if(!missing(modelsupport)){
    if(missing(model)){model <-1; warning("No model provided, selecting best model by default")}
    if(!missing(BPs)){warning("modelsupport provided, breakpoints from BPs argument ignored")}
    bestmodel <- modelsupport[model,]  #bestmodel is first row from analysis
    nBPs <- bestmodel$regions-1
    keep <- grep('breakpoint', colnames(bestmodel), value=TRUE)
    BPs <- unlist(bestmodel[,keep])
  } else if (!missing(BPs)){
    nBPs <- length(BPs)
  } else {stop("No modelsupport or BPs provided")}
  firstvert <- Xvar[1]
  lastvert <- Xvar[length(Xvar)]
  Yvar <- data[,pcono]
  df <- as.data.frame(cbind(Xvar, Yvar))

  if(nBPs ==0){
    form <- "Yvar ~ Xvar"

  } else {
    if(cont==F){	# Discontinuous fit
	beg <- paste("Yvar ~ Xvar*(Xvar<=",BPs[1],")", sep="")
	end <- paste("+Xvar*(Xvar>",BPs[nBPs],")", sep="")

	if(nBPs==1){
	   form <- paste(beg, end, sep="")
	} else if(nBPs >1){
	   f <- c()
	   for (i in 1:(nBPs-1)){
	      f[i] <- paste("+ Xvar*(Xvar>",BPs[i]," & Xvar<=", BPs[i+1], ")", sep="")
	   }
	   form <- paste(beg, paste(f, collapse=""), end, sep="")
	}
    } else {		# Continuous fit
	f <- c()
	for (i in 1:nBPs){
	   f[i] <- paste("+ pmax(0,Xvar-",BPs[i],")", sep="")
	}
	form <- paste("Yvar ~ Xvar", paste(f, collapse=""), sep="")
    }
  }

  fit <- lm(form, data=df)

  label.pco<-paste0("PCO",pcono)

  if(missing(xlim)){xlim <- c(min(Xvar), max(Xvar))}
  if(missing(ylim)){ylim <- c(min(Yvar), max(Yvar))}
  if(missing(title)){title <- label.pco}
  if(missing(pch)){pch <- 16}
  if(missing(col.pts)){col.pts <- "darkgrey"}
  if(missing(col.fit)){col.fit <- "darkturquoise"}
  if(missing(col.BP)){col.BP <- "coral"}
  if(missing(lwd.fit)){lwd.fit <- 2}
  if(missing(lwd.BP)){lwd.BP <- 2}

  BP.pos <- which(Xvar %in% BPs)	# Added to account for potential missing values in Xvar

  plot(Xvar, Yvar, pch=pch, col=col.pts, xlab="Vertebral position", ylab = label.pco, xlim=xlim, ylim=ylim)
  title(main = title)

  if(nBPs==0){
	lines(Xvar,predict(fit), col=col.fit, lwd=lwd.fit)
  } else if(nBPs==1){
	lines(Xvar[1:BP.pos[1]],predict(fit)[1:BP.pos[1]], col=col.fit, lwd=lwd.fit)
	lines(Xvar[(BP.pos[1]+1):length(Xvar)],predict(fit)[(BP.pos[1]+1):length(Xvar)],col=col.fit, lwd=lwd.fit)
	abline(v = BPs[1] + 0.5, col = col.BP, lty = 3, lwd=lwd.BP)
  } else {
	lines(Xvar[1:BP.pos[1]],predict(fit)[1:BP.pos[1]], col=col.fit, lwd=lwd.fit)
	for(i in 1:(nBPs-1)){
	   lines(Xvar[(BP.pos[i]+1):BP.pos[i+1]],predict(fit)[(BP.pos[i]+1):BP.pos[i+1]], col=col.fit, lwd=lwd.fit)
	}
	lines(Xvar[(BP.pos[nBPs]+1):length(Xvar)],predict(fit)[(BP.pos[nBPs]+1):length(Xvar)], col=col.fit, lwd=2)
	for(i in 1:nBPs){abline(v = BPs[i] + 0.5, col = col.BP, lty = 3, lwd=lwd.BP)}
  }
}






###########
## Function: plot vertebral map of region ##
###########


plotvertmap <- function(name, Xvar, modelsupport, model=NULL, model.nreg=NULL, BPs, drop.na=FALSE, plotType, centraL=NULL, col.by.block=FALSE, cols=NULL, blocklim=NULL, border.col="white", reglimits=NULL, lim.col=NULL, text=FALSE, bp.sd=NULL, sd.col=NULL, sd.jit=0.1){
		# name = character with specimen/species name
		# Xvar = vector with vertebral number included in the analysis (missing vertebrae should be skipped, not with NA)
		# modelsupport (optional) = data.frame with best models ordered according to AICcW (object Model_support from model_support function output) /!\ modelsupport and BPs optional but at least one of them must be provided!
		# model (optional) = model to plot (1=best, 2=2nd best, etc), if not specified, will plot best model
		# model.nreg (optional) = model to plot selected by number of regions (1 = 1 region model, 2 = 2 regions models, etc), if specified, will overwritte "model" argument
		# BPs (OPTIONAL): can provide vector with ordered breakpoints directly if modelsupport not provided - if no BPs (ie, 1 region) use "NULL"
		# drop.na = how to treat missing vertebrae, if FALSE (default) then, missing vertebrae are added to the mapping and grayed, if TRUE, missing vertebrae are ignored and specimen is plotted as given
		# plotType = "count" OR "percent" OR "length" to choose the type of value on X axis (vertebral count or % vertebral count or % total length). If "length" need to provide length of all centra with "centraL"
		# centraL (has to be provided is plotType==length) = vecotr with length of centrum of all vertebrae included in plot
		# col.by.block (optional) = TRUE/FALSE wether to group region colors by predefined blocks. If FALSE (default), cols should be a vector. If TRUE, cols should be a named list in which each element is a vector with colors to use for each region, names of the list will be used as names for blocks. 
		# cols (optional) = vector of colors to use for plotting (length must be equal to number of regions)
		# blocklim = Only needed if "col.by.block" is TRUE. Vector with position (vertebral number) of limit between blocks (limit between 2 blocks must corresponds to the number of the last vertebra of the anterior block) 
		# border.col (optional) = color to use for the border of each vertebra, if set to "NA" no border will be drawn (default to white)
		# reglimits (optional) = vector of traditional regions limits to plot (limit between 2 regions must correspond to the number of the last vertebra of the region)
		# lim.col (optional) = color to be used to plot traditional regions limits. Default is black
		# text (optional) = TRUE/FALSE wether to plot the vertebral number on each square. Default is FALSE
		# bp.sd (optional) = standard deviation (in number of vertebrae no percent) of regions limits
		# sd.col (optional) = if bp.sd provided, color to be used to plot standard deviation of regions limits. Default is black
		# sd.jit (optional) = amount of vertical jitter to add to sd bars if they overlap, value between 0 & 1 (0= no jitter, 1=maximal jitter)

 library(RColorBrewer)

  if(plotType=="length" & is.null(centraL)){
    stop("Provide length of vertebral centra for length plotType option")
  }

  if(!drop.na){  # Fill missing vertebrae:
    Xvar.all <- 1:max(Xvar)
  } else {
    Xvar.all <- setNames(1:length(Xvar), Xvar)
  }

  if(!missing(modelsupport)){
    if(is.null(model) & is.null(model.nreg)){model <-1; warning("No model provided, selecting best model by default")}
    if(!missing(BPs)){warning("modelsupport provided, breakpoints from BPs argument ignored")}
    if(!is.null(model.nreg)){ best <- modelsupport[modelsupport$regions==model.nreg,]	# If model.nreg provided, select model with determined number of regions
    } else { best <- modelsupport[model,] }		# If model.nreg not provided, select model by order of best support
    nreg <- unlist(best$regions)
    BPs <- unlist(best[2:nreg])
  } else if (!missing(BPs)){
    BPs <- round(BPs)
    nreg <- length(BPs)+1
  } else {stop("No modelsupport or BPs provided")}

  if(col.by.block==TRUE & is.null(blocklim)){stop("No block limits provided")}

  if(sd.jit > 1 | sd.jit <0){stop("SD jittering value must be between 0 and 1")}

 # Set beginning position of each vert in count:
  vertmap <- as.data.frame(cbind(Xvar.all, Xvar.all.deb=Xvar.all-1))

 # Calculate % vertebral count:
  vertmap$pct <- vertmap$Xvar.all/max(Xvar.all)*100

  # Set beginning position of each vert in %:
  vertmap$pct.deb <- vertmap$pct
  for(i in 2:nrow(vertmap)){vertmap$pct.deb[i] <- vertmap$pct[i-1]}
  vertmap$pct.deb[1] <- 0	# Set beginning of first vert as 0


 # Set beginning & end position of each vert in % length centra:
  if(plotType=="length"){
    Lc.end <- centraL[1]
    Lc.beg <- 0
    for(i in 2:length(centraL)){
      Lc.end[i] <- centraL[i]+Lc.end[i-1]
      Lc.beg[i] <- Lc.end[i-1]
    }
    TCL <- sum(centraL)
    Lc.end <- Lc.end/TCL*100
    Lc.beg <- Lc.beg/TCL*100
    vertmap <- cbind(vertmap,Lc.beg,Lc.end)
  }

 # Add corresponding region to each vertebra:
  vertmap$reg <- rep(NA, length(Xvar.all))
  regions <- paste0("Region", 1:nreg)
  if(!drop.na){
    vertmap[which(vertmap$Xvar.all>= min(Xvar)), "reg"] <- regions[1]	# Place 1st regions
    for(i in 1:length(BPs)){							  	# Add other regions
	bp <- BPs[i]
	vertmap[which(vertmap$Xvar.all > bp), "reg"] <- regions[i+1]
    }
  } else {
    vertmap$Xvar <- Xvar
    vertmap[which(vertmap$Xvar>= min(Xvar)), "reg"] <- regions[1]	# Place 1st regions
    for(i in 1:length(BPs)){							  	# Add other regions
	bp <- BPs[i]
	vertmap[which(vertmap$Xvar > bp), "reg"] <- regions[i+1]
    }
  }

 # Specify missing vert:
  if(!drop.na){vertmap[!Xvar.all %in% Xvar,"reg"] <- "Missing"}

 # Add indiv identity:
  vertmap <- cbind(indiv=rep(name, nrow(vertmap)), vertmap)


 ## Define position of blocks if requested:
  if(col.by.block){
    reg.names <- names(cols)
    vertmap$trad.reg <- rep(NA, nrow(vertmap))
    # Add traditional regions to table:
    for(i in 1:length(reg.names)){		
      if(i==1){vertmap[vertmap$Xvar<=blocklim[i],"trad.reg"] <- reg.names[i]
      } else if(i==length(reg.names)){vertmap[vertmap$Xvar>blocklim[i-1],"trad.reg"] <- reg.names[i]
      } else {vertmap[vertmap$Xvar<=blocklim[i] & vertmap$Xvar>blocklim[i-1],"trad.reg"] <- reg.names[i]}
    }

    # Match traditional regions boundaries and new regions:
    regs <- levels(as.factor(vertmap$reg))
    vertmap$tradreg.corr <- rep(NA, nrow(vertmap))
    for(i in 1:length(regs)){
      t <- table(vertmap[vertmap$reg==regs[i],"trad.reg"])
      vertmap[vertmap$reg==regs[i],"tradreg.corr"] <- names(which.max(t))
    }

    # Name new regions based on their appartenance to traditional regions:
    reg.corr <- unique(vertmap$tradreg.corr)
    vertmap$Region <- rep(NA,nrow(vertmap))
    for(i in 1:length(reg.corr)){
      reg_sub <- unique(vertmap[vertmap$tradreg.corr==reg.corr[i],"reg"])
      for(j in 1:length(reg_sub)){
  	vertmap[vertmap$reg==reg_sub[j],"Region"] <- paste0(reg.corr[i],j)
      }
    }
  } else {colnames(vertmap)[colnames(vertmap)=="reg"] <- "Region"}


 ## Traditional regions limit position ##

  if(!is.null(reglimits)){
    if (is.null(lim.col)){
	lim.col <- "black"
    }
    reglimits <- as.data.frame(reglimits)
    colnames(reglimits) <- "lims"

    if(drop.na){
    for(m in 1:nrow(reglimits)){
      reglimits[m,] <- vertmap[max(which(vertmap$Xvar<=reglimits[m,])),"Xvar.all"]	# Replace limit position if subsampling
    }}

 
   # Calculate limits in % count:
    reglimits.pct <- reglimits/max(Xvar.all)*100

   # Calculate limits in % centra length:
    if(plotType=="length"){
      reglimits.Lc <- c()
      for(i in 1:nrow(reglimits)){
        bound <- reglimits[i,]
        if(is.na(bound)){
          reglimits.Lc[i] <- NA
        } else {
          reglimits.Lc[i] <- vertmap[vertmap$Xvar.all==bound,"Lc.end"]
        }
        names(reglimits.Lc)[i] <- rownames(reglimits)[i]
      }
      reglimits.Lc <- data.frame(reglimits.Lc); colnames(reglimits.Lc) <- "lims"
    }
  }

 # Get BP corrected position if drop.na vertebrae:
  if(drop.na){
	BPplot <- vertmap[(vertmap$Xvar %in% BPs),"Xvar.all"]
  } else {
	BPplot <- BPs
  }


 # Correct BP sd position according to % count and % length:
  if(!is.null(bp.sd)){
    if(is.null(sd.col)){sd.col <- "black"}
    sd.jit <- sd.jit/2

   # count:
    bp.sd.c <- data.frame(cbind(BPs=BPplot,bp.sd)); 
    bp.sd.c$beg <- bp.sd.c$BPs - bp.sd.c$bp.sd; 
    bp.sd.c$end <- bp.sd.c$BPs + bp.sd.c$bp.sd

   # test for overlap of sd bars & add jitter if needed:
    overlap <- sapply(seq(nrow(bp.sd.c)), function(i) {with(bp.sd.c[seq(i),], length(which(beg < beg[i] & end > beg[i]))) })
    for(i in 1:length(overlap)){if((i%%2)==0){bp.sd.c$y[i] <- 1.5-sd.jit*overlap[i]}else{bp.sd.c$y[i] <- 1.5+sd.jit*overlap[i]}}

   # percent:
    bp.sd.pct <- data.frame(cbind(BPs=(BPplot/max(Xvar.all)*100),bp.sd=(bp.sd/max(Xvar.all)*100)))
    bp.sd.pct$beg <- bp.sd.pct$BPs - bp.sd.pct$bp.sd; bp.sd.pct$end <- bp.sd.pct$BPs + bp.sd.pct$bp.sd
    for(i in 1:length(overlap)){if((i%%2)==0){bp.sd.pct$y[i] <- 1.5-sd.jit*overlap[i]}else{bp.sd.pct$y[i] <- 1.5+sd.jit*overlap[i]}}

    if(plotType=="length"){
     # % length:
      bp.sd.Lc <- data.frame(matrix(ncol=4,nrow=length(bp.sd)))
      colnames(bp.sd.Lc) <- c("BPs","bp.sd","beg","end")
      for(i in 1:length(bp.sd)){
        sd <- bp.sd[i]
        bp <- BPplot[i]
        bp.sd.Lc$BPs[i] <- vertmap[vertmap$Xvar.all %in% bp,"Lc.end"] 
        if(sd <= 1){
          vant <- vertmap[vertmap$Xvar.all %in% (bp),c("Lc.beg","Lc.end")]
          pant <- sd*(vant$Lc.end-vant$Lc.beg)
          bp.sd.Lc$beg[i] <- bp.sd.Lc$BPs[i] - pant
          vpost <- vertmap[vertmap$Xvar.all %in% (bp+1),c("Lc.beg","Lc.end")]
          ppost <- sd*(vpost$Lc.end-vpost$Lc.beg)
          bp.sd.Lc$end[i] <- bp.sd.Lc$BPs[i] + ppost
        } else{
          sdfull <- floor(sd)
          sdres <- sd-sdfull
          vant <- vertmap[vertmap$Xvar.all %in% (bp-(sdfull)),c("Lc.beg","Lc.end")]
          pant <- sdres*(vant$Lc.end-vant$Lc.beg)
          bp.sd.Lc$beg[i] <- vant$Lc.end - pant
          vpost <- vertmap[vertmap$Xvar.all %in% (bp+(sdfull+1)),c("Lc.beg","Lc.end")]
          ppost <- sdres*(vpost$Lc.end-vpost$Lc.beg)
          bp.sd.Lc$end[i] <- vpost$Lc.beg + ppost
        }
      }
      overlap <- sapply(seq(nrow(bp.sd.Lc)), function(i) {with(bp.sd.Lc[seq(i),], length(which(beg < beg[i] & end > beg[i]))) })
      for(i in 1:length(overlap)){if((i%%2)==0){bp.sd.Lc$y[i] <- 1.5-sd.jit*overlap[i]}else{bp.sd.Lc$y[i] <- 1.5+sd.jit*overlap[i]}}
    }
  }

 ## Set color for each region
  if(col.by.block){
    colors <- unlist(cols)
    colors <- c(Missing="grey", colors)
  } else {
   # Set color for each region
    if (!is.null(cols)){
  	pal <- cols
    } else {
  	if (nreg %% 2 != 0){
  	  pal <- brewer.pal((nreg+1), "Paired")
  	  pal <- pal[-(length(pal)/2)]
	} else {
	  pal <- brewer.pal(nreg, "Paired")
	}
    }
    colors <- pal[1:nreg]
    names(colors) <- regions
    colors <- c(Missing="grey", colors)		# Add color for missing vertebrae
  }


 # Vertebrae names for plot:
  if(!drop.na){
    vname <- Xvar.all
  } else {
    vname <- as.numeric(names(Xvar.all))
  }
  vertmap$vname <- vname

 # Create plot
  if(plotType=="count"){
    p <- ggplot(vertmap)+
	geom_rect(vertmap, mapping=aes(xmin=Xvar.all.deb, xmax=Xvar.all, ymin=1, ymax=2, fill=Region),color=border.col)+ 
	scale_fill_manual("Legend", values = colors) +
	labs(y=unique(vertmap$indiv), x="Vertebral count")+
	scale_x_continuous(expand=c(0,0))+
	theme_classic()+
	theme(axis.text.y=element_blank(),
		axis.ticks.y=element_blank(),
		axis.line.y=element_blank(),
		axis.line.x=element_blank(),
		axis.title.y = element_text(angle=0, size=10, vjust=0.5, hjust=1),
		legend.position="none",
		plot.margin=unit(c(0,10,0,0), "pt"))
    if(!is.null(reglimits)){p <- p + geom_vline(data=reglimits, aes(xintercept=lims), color=lim.col, lwd=1)}
    if(text){p <- p+geom_text(aes(x=(Xvar.all.deb + Xvar.all)/2, y = 1.5, label = vname))}   # Add vertebral number on each colored square
    if(!is.null(bp.sd)){													# Add horizontal line for BPs standard deviation
	p <- p+geom_segment(data=bp.sd.c,aes(x=beg,xend=end,y=y,yend=y), lwd=1, col=sd.col)
    }

  } else if(plotType=="percent"){
    p <- ggplot(vertmap)+
	geom_rect(vertmap, mapping=aes(xmin=pct.deb, xmax=pct, ymin=1, ymax=2, fill=Region),color=border.col)+ 
	scale_fill_manual("Legend", values = colors) +
	labs(y=unique(vertmap$indiv), x="% Vertebral count")+
	scale_x_continuous(expand=c(0,0))+
	theme_classic()+
	theme(axis.text.y=element_blank(),
		axis.ticks.y=element_blank(),
		axis.line.y=element_blank(),
		axis.line.x=element_blank(),
		axis.title.y = element_text(angle=0, size=10, vjust=0.5, hjust=1),
		legend.position="none",
		plot.margin=unit(c(0,10,0,0), "pt"))
    if(!is.null(reglimits)){p <- p + geom_vline(data=reglimits, aes(xintercept=lims), color=lim.col, lwd=1)}
    if(text){p <- p+geom_text(aes(x=(pct.deb + pct)/2, y = 1.5, label = vname))}   
    if(!is.null(bp.sd)){
	p <- p+geom_segment(data=bp.sd.pct,aes(x=beg,xend=end,y=y,yend=y), lwd=1, col=sd.col)
    }
  } else if(plotType=="length") {
    p <- ggplot(vertmap)+
	geom_rect(vertmap, mapping=aes(xmin=Lc.beg, xmax=Lc.end, ymin=1, ymax=2, fill=Region),color=border.col)+ 
	scale_fill_manual("Legend", values = colors) +
	labs(y=unique(vertmap$indiv), x="% Total centrum length")+
	scale_x_continuous(expand=c(0,0))+
	theme_classic()+
	theme(axis.text.y=element_blank(),
		axis.ticks.y=element_blank(),
		axis.line.y=element_blank(),
		axis.line.x=element_blank(),
		axis.title.y = element_text(angle=0, size=10, vjust=0.5, hjust=1),
		legend.position="none",
		plot.margin=unit(c(0,10,0,0), "pt"))
    if(!is.null(reglimits)){p <- p + geom_vline(data=reglimits, aes(xintercept=lims), color=lim.col, lwd=1)}
    if(text){p <- p+geom_text(aes(x=(Lc.beg + Lc.end)/2, y = 1.5, label = vname))}   
    if(!is.null(bp.sd)){
	p <- p+geom_segment(data=bp.sd.Lc,aes(x=beg,xend=end,y=y,yend=y), lwd=1, col=sd.col)
    }
  } else {stop("invalid plotType: choose count, percent or length")}
  return(p)
}





###########
## Function: plot region score of eahc PCO ##
###########


plotpcoregAG <-function(eigenvals, nvert, regiondata, cont, title=NULL){
	# eigenvals: eigenvalues from the PCO (or other) ($eigen.val value from svdPCO function output)
	# nvert: number of vertebrae included in the analysis (must be the same as calcregionsAG)
	# regiondata: $results value from calcregionsAG output
	# nPCO: number of PCOs used in the calcregionsAG analysis
	# cont (TRUE/FALSE): if continuous (TRUE) or discontinuous (FALSE) fit was use to calculate regiondata
	# title (optional): title to give to the graph

  require(ggplot2)
  require(tidyr)

  noregions <- max(regiondata$regions)
  nPCO <- length(grep("RSS." ,colnames(regiondata)))

  # Calculate variance of each PCO:
  var.exp <- cumsum((eigenvals/sum(eigenvals)*100))[1:nPCO]

  pco.no.test<-data.frame(matrix(NA,nrow=nPCO,ncol=6, dimnames=list(c(1:nPCO),c("PCO", "RSind.AICc", "RScum.AICc", "RSind.BIC", "RScum.BIC", "CumulVar"))))


  for (a in 1:nPCO){

    #Run for individual PCs
    models.ind<-modelselectAG(regiondata, PCOs.no=a)
    support.ind<-model_supportAG(models.ind,nvert,cont)

    #Run for cumulative PCs
    models.cum<-modelselectAG(regiondata, PCOs.no=1:a)
    support.cum<-model_supportAG(models.cum,nvert, cont)

    pco.no.test[a,1]<-a
    pco.no.test[a,2]<-support.ind$Region_score
    pco.no.test[a,3]<-support.cum$Region_score
    pco.no.test[a,4]<-support.ind$Region_score_BIC
    pco.no.test[a,5]<-support.cum$Region_score_BIC
  }

  pco.no.test[,6]<-var.exp
  pco.no.test.long <- pivot_longer(pco.no.test, cols=2:5, names_to=c("PCOtype","Testtype"),names_sep="\\.")

  p <- ggplot(pco.no.test.long,aes(x=PCO,y=value))+
    geom_point(aes(col=PCOtype), pch=19)+
    scale_color_manual(values=c("#fc8d62","#8da0cb"), labels=c("Cumulated PCOs","Single PCO"))+
    geom_smooth(formula=y~x,method="loess",aes(x=PCO, y=(CumulVar/(100/noregions))), se=F, col="darkgrey")+
    scale_y_continuous(name="Region score", sec.axis=sec_axis(~.*(100/noregions), name="Cumultaed variance explained (%)"))+
    facet_wrap(~Testtype)+
    theme_bw()+
    theme(panel.grid.minor=element_blank(), legend.position="bottom")+
    labs(color="Region score")

  if(!is.null(title)){p <- p+labs(title=title)}

  plot(p)
  return(table=pco.no.test)
}










