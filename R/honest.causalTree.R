#
#  The honest re-estimation function.
#

honest.causalTree <- function(formula, data, weights, treatment, subset, 
							  est_data, est_weights, est_treatment, est_subset,
							  na.action = na.causalTree, split.Rule, split.Honest,
							  HonestSampleSize, split.Bucket, bucketNum = 10,
							  bucketMax = 40, cv.option, cv.Honest, minsize = 2L, model = FALSE,
							  x = FALSE, y = TRUE, propensity, control, split.alpha = 0.5, 
							  cv.alpha = 0.5,cv.gamma=0.5,split.gamma=0.5, cost, ...)  { 

	Call <- match.call()

	indx <- match(c("formula", "data", "weights", "subset"),
				  names(Call), nomatch = 0L)

	if (indx[1] == 0L) stop("a 'formula' argument is required")
	temp <- Call[c(1L, indx)]      
	temp$na.action <- na.action  
	temp[[1L]] <- quote(stats::model.frame) 
	m <- eval.parent(temp)


	Terms <- attr(m, "terms")
	if (any(attr(Terms, "order") > 1L))
		stop("Trees cannot handle interaction terms")

	Y <- model.response(m)
	wt <- model.weights(m)
	if (any(wt < 0)) stop("negative weights not allowed")
	if (!length(wt)) wt <- rep(1, nrow(m))
	offset <- model.offset(m)
	X <- causalTree.matrix(m)

	nobs <- nrow(X)
	nvar <- ncol(X)
	#treatment <- m$`(treatment)`

	# requirement for treatment status
	if (missing(treatment)) {
		stop("You should input the treatment status vector for data:
			 1 represent treated and 0 represent controlled.")   
	}
	if (sum(treatment %in% c(0,1)) != nobs) {
		stop("The treatment status should be 1 or 0 only: 1 represent treated and 0 represent controlled.")
	}
	if (sum(treatment) == 0 || sum(treatment) == nobs) {
		stop("The data only contains treated cases or controlled cases, please check 'treatment' again.") 
	}

	# ---------------------------------------------------------------------------------------
	# check the honest re-estimation data set:
	if (missing(est_data)) {
		stop("Note give the honest estimation data set!\n")
	}
	
	indx2 <- match(c("formula", "est_data", "est_weights", "est_subset"),
				   names(Call), nomatch = 0L)

	temp2 <- Call[c(1L, indx2)]
	temp2$na.action <- na.action
	names(temp2) <- gsub("est_", "", names(temp2))
	temp2[[1L]] <- quote(stats::model.frame)
	m2 <- eval.parent(temp2)
	# honest data set used for later:
	est_Y <- model.response(m2)
	est_wts <- model.weights(m2)
	if (any(est_wts < 0)) stop("negative weights not allowed")
	if (!length(est_wts)) est_wts <- rep(1, nrow(m2))
	est_offset <- model.offset(m2)
	est_X <- causalTree.matrix(m2)
	est_nobs <- nrow(est_X)
	est_nvar <- ncol(est_X)
	
	if (missing(est_treatment)) {
	    stop("Not give the treatment status of honest estimation data set!\n ")
	}
	if (sum(est_treatment %in% c(0,1)) != est_nobs) {
	    stop("The treatment status should be 1 or 0 only: 1 represent treated and 0 represent controlled.")
	}
	if (sum(est_treatment) == 0 || sum(est_treatment) == est_nobs) {
	    stop("The data only contains treated cases or controlled cases, please check 'est_treatment' again.") 
	}
	

	if (est_nvar != nvar) {
		stop("Honest estimation data set should have same variables as training data!\n")
	}
	# ---------------------------------------------------------------------------------------
	if (missing(propensity)) {
		propensity <- sum(treatment) / nobs
	}

	## check the Split.Rule:
	if (missing(split.Rule)) {
		split.Rule <- "TOT"
		warning("The default split rule is 'TOT'.")
	}

	# check split.Bucket:
	if (missing(split.Bucket)) {
		split.Bucket <- FALSE
		bucketNum <- 0
		bucketMax <- 0
	} 
	split.Bucket.num <- pmatch(split.Bucket, c(T, F))
	if (is.na(split.Bucket.num))
		stop("Invalid split.Honest input, split.Honest can be only TRUE or FALSE.")

	if (!split.Bucket) {
		# split.Bucket = F
		bucketNum <- 0
		bucketMax <- 0
	} else {
		# split.Bucket = T
		if (missing(bucketMax)) {
			# maximum number of buckets
			bucketMax <- 100
		} 
		if (missing(bucketNum)) {
			# number of treat cases or control cases in one bucket
			# Numbuckets = max(minsize, min(round(numtreated/bucketNum),round(numcontrol/bucketNum),bucketMax))
			bucketNum <- 5
			# at least 5 obs in one bucket
		}
		split.Rule <- paste(split.Rule, 'D', sep = '') 
	}

	split.Rule.int <- pmatch(split.Rule, c("TOT", "CT", "fit", "tstats", "TOTD", "CTD", "fitD", "tstatsD", "user", "userD","policy","policyD"))
	if (is.na(split.Rule.int)) stop("Invalid splitting rule.")
	split.Rule <- c("TOT", "CT", "fit", "tstats", "TOTD", "CTD", "fitD", "tstatsD", "user", "userD","policy","policyD")[split.Rule.int]
   print(split.Rule.int)
   print(split.Rule)
	## check the Split.Honest, for convenience
	if (split.Rule.int %in% c(1, 5)) {
		if (!missing(split.Honest)) {
			warning("split.Honest is not used in your chosen splitting rule.")
		}
		if (!missing(split.alpha)) {
			warning("split.alpha is not used in your chosen splitting rule. split.Honest set to FALSE")
		}
		split.Honest <- FALSE
	} else {
		if (missing(split.Honest)) {
			split.Honest <- TRUE
			warning("The default split.Honest = TRUE for your chosen splitting rule.")
		}
	}

	## check the Split.Honest == T/F
	split.Honest.num <- pmatch(split.Honest, c(T, F))
	if(is.na(split.Honest.num)) 
		stop("Invalid split.Honest input, split.Honest can be only TRUE or FALSE.")

	if (split.Honest == TRUE && split.Rule.int %in% c(2, 3, 4, 6, 7, 8, 9, 10,11,12)) {
		# ct, fit, tstats, ctd, fitd, tstatsd, user, userd,policy,policyD:
		if(missing(split.alpha)) {
			# set default honest splitting alpha to 0.5
			split.alpha <- 0.5
		} else {
			# check split.alpha in [0, 1]
			if (split.alpha > 1 || split.alpha < 0) {
				stop("Invalid input for split.alpha. split.alpha should between 0 and 1.")
			}
		}
	  #check for gamma for policy
	  if(missing(split.gamma)) {
	    # set default honest splitting alpha to 0.5
	    split.gamma <- 0.5
	  } else {
	    # check split.alpha in [0, 1]
	    if (split.gamma > 1 || split.gamma < 0) {
	      stop("Invalid input for split.gamma. split.gamma should between 0 and 1.")
	    }
	  }
	} else if (split.Rule.int %in% c(2, 3, 4, 6, 7, 8, 9, 10,11,12)){
		# split.Honest = False
		if (split.alpha != 1) 
			warning("For dishonest(adaptive) splitting, split.alpha =  1.");
		split.alpha <- 1
		#added gamma check for policy
		if(missing(split.gamma)) {
		  # set default honest splitting alpha to 0.5
		  split.gamma <- 0.5
		} else {
		  # check split.alpha in [0, 1]
		  if (split.gamma > 1 || split.gamma < 0) {
		    stop("Invalid input for split.gamma. split.gamma should between 0 and 1.")
		  }
		}
		
	}


	# check propensity score:
	if (split.Rule %in% c("TOT", "TOTD", "fit", "fitD")) {
		if (propensity > 1 || propensity < 0) {
			stop("Propensity score should be between 0 and 1.")
		}
	}

	xvar <- apply(X, 2, var)
	method <- "anova"
	method.int <- 1

	# ------------------------------------- cv begins -------------------------------------- #
	if (missing(cv.option)) {
		# temporarily, no crossvalidation 
		warning("Miss 'cv.option', choose not to do cross validations.")
		cv.option <-"none"
		xval <- 0
	} 

	if(missing(cv.Honest)) {
		cv.Honest <- TRUE
	}
	cv.Honest.num <- pmatch(cv.Honest, c(T, F))
	if (is.na(cv.Honest.num)) 
		stop ("Invalid cv.Honest. cv.Honest should be TRUE or FALSE.")

	if (cv.option == 'CT' || cv.option == 'fit') {
		if (cv.Honest) {
			# cv.Honest = T
			cv.option <- paste(cv.option, 'H', sep = '')
		} else {
			cv.option <- paste(cv.option, 'A', sep = '')
		}
	}

	cv.option.num <- pmatch(cv.option, c("TOT", "matching", "fitH", "fitA", "CTH", "CTA", "userH", "userA", "none"))
	if(is.na(cv.option.num)) stop("Invalid cv option.") 

	# check cv.alpha
	if (cv.option.num %in% c(1, 2, 4, 6)) {
		if (!missing(cv.alpha))
			warning("cv.alpha is not used in your chosen cross validation method.")
	} 

	if (missing(cv.alpha)) {
		cv.alpha <- 0.5
	}
	#for policy, set gamma (set for all presently)
	if (missing(cv.gamma)) {
	  cv.gamma <- 0.5
	}
	if (missing(HonestSampleSize)) {
		# default should be the # of samples in training 
		HonestSampleSize <- est_nobs
	}

	if(HonestSampleSize != est_nobs) {
		warning("HonestSampleSize shoud be the number of observations in estimation sample.")
		HonestSampleSize <- est_nobs
	}

	HonestSampleSize <- as.integer(HonestSampleSize)
	if (is.na(HonestSampleSize))
		stop("HonestSampleSize should be an integer.")
	# -------------------------------- cv checking ends -------------------------------- #

	init <- get(paste("causalTree", method, sep = "."), envir = environment())(Y, offset, wt) 

	ns <- asNamespace("causalTree")
	if (!is.null(init$print)) environment(init$print) <- ns
	if (!is.null(init$summary)) environment(init$summary) <- ns
	if (!is.null(init$text)) environment(init$text) <- ns

	Y <- init$y

	xlevels <- .getXlevels(Terms, m)
	cats <- rep(0L, ncol(X))
	if (!is.null(xlevels))
		cats[match(names(xlevels), colnames(X))] <-
			unlist(lapply(xlevels, length))

		extraArgs <- list(...)
		if (length(extraArgs)) {
			controlargs <- names(formals(rpart.control)) # legal arg names
			indx <- match(names(extraArgs), controlargs, nomatch = 0L)
			if (any(indx == 0L))
				stop(gettextf("Argument %s not matched",
							  names(extraArgs)[indx == 0L]),
					 domain = NA)
		}

		controls <- causalTree.control(...)
		if (!missing(control)) controls[names(control)] <- control

		xval <- controls$xval
		if (is.null(xval) || (length(xval) == 1L && xval == 0L) || method=="user") {
			xgroups <- 0L
			xval <- 0L
		} else if (length(xval) == 1L) {
			# make random groups
			control_idx <- which(treatment == 0)
			treat_idx <- which(treatment == 1)
			xgroups <- rep(0, nobs)
			xgroups[control_idx] <- sample(rep(1L:xval, length = length(control_idx)), length(control_idx), replace = F)
			xgroups[treat_idx] <- sample(rep(1L:xval, length = length(treat_idx)), length(treat_idx), replace = F)  
		} else if (length(xval) == nobs) {
			## pass xgroups by xval 
			xgroups <- xval
			xval <- length(unique(xgroups))
		} else {
			## Check to see if observations were removed due to missing
			if (!is.null(attr(m, "na.action"))) {
				## if na.causalTree was used, then na.action will be a vector
				temp <- as.integer(attr(m, "na.action"))
				xval <- xval[-temp]
				if (length(xval) == nobs) {
					xgroups <- xval
					xval <- length(unique(xgroups))
				} else stop("Wrong length for 'xval'")
			} else stop("Wrong length for 'xval'")
		}


		##
		## Incorprate costs
		##
		if (missing(cost)) cost <- rep(1, nvar)
		else {
			if (length(cost) != nvar)
				stop("Cost vector is the wrong length")
		if (any(cost <= 0)) stop("Cost vector must be positive")
		}

		##
		## Have C code consider ordered categories as continuous
		##  A right-hand side variable that is a matrix forms a special case
		## for the code.
		##
		tfun <- function(x)
			if (is.matrix(x)) rep(is.ordered(x), ncol(x)) else is.ordered(x)

		labs <- sub("^`(.*)`$", "\\1", attr(Terms, "term.labels")) # beware backticks
		isord <- unlist(lapply(m[labs], tfun))

		storage.mode(X) <- "double"
		storage.mode(wt) <- "double"
		storage.mode(treatment) <- "double"
		minsize <- as.integer(minsize) # minimum number of obs for treated and control cases in one leaf node

		ctfit <- .Call(C_causalTree,
					   ncat = as.integer(cats * !isord),
					   split_Rule = as.integer(split.Rule.int), # tot, ct, fit, tstats, totD, ctD, fitD, tstatsD
					   bucketNum = as.integer(bucketNum), # if == 0, no discrete; else do discrete
					   bucketMax = as.integer(bucketMax),
					   method = as.integer(method.int),  # "anova" not changed yet
					   crossmeth = as.integer(cv.option.num), # tot, ct, fit, tstats
					   crossHonest = as.integer(cv.Honest.num),
					   as.double(unlist(controls)), # control list in rpart
					   minsize, # minsize = min_node_size
					   as.double(propensity),
					   as.integer(xval),
					   as.integer(xgroups),
					   as.double(t(init$y)),
					   X, # X features for model data
					   wt, # for model data
					   treatment, # for model data
					   as.integer(init$numy),
					   as.double(cost),
					   as.double(xvar), # for model daa
					   as.double(split.alpha),
					   as.double(cv.alpha),
					   as.integer(HonestSampleSize),
					   as.double(cv.gamma)
					   )

		nsplit <- nrow(ctfit$isplit) # total number of splits, primary and surrogate
		## total number of categorical splits
		ncat <- if (!is.null(ctfit$csplit)) nrow(ctfit$csplit) else 0L

		if (nsplit == 0L) xval <- 0L # No xvals were done if no splits were found

		numcp <- ncol(ctfit$cptable)
		temp <- if (nrow(ctfit$cptable) == 3L) c("CP", "nsplit", "rel error")
			else c("CP", "nsplit", "rel error", "xerror", "xstd")
		dimnames(ctfit$cptable) <- list(temp, 1L:numcp)

		tname <- c("<leaf>", colnames(X))
		splits <- matrix(c(ctfit$isplit[, 2:3], ctfit$dsplit), ncol = 5L,
						 dimnames = list(tname[ctfit$isplit[, 1L] + 1L],
										 c("count", "ncat", "improve", "index", "adj")))
		index <- ctfit$inode[, 2L]  # points to the first split for each node

		## Now, make ordered factors look like factors again (a printout choice)
		nadd <- sum(isord[ctfit$isplit[, 1L]])
		if (nadd > 0L) { # number of splits at an ordered factor.
			newc <- matrix(0L, nadd, max(cats))
			cvar <- ctfit$isplit[, 1L]
			indx <- isord[cvar]             # vector of TRUE/FALSE
			cdir <- splits[indx, 2L]        # which direction splits went
			ccut <- floor(splits[indx, 4L]) # cut point
			splits[indx, 2L] <- cats[cvar[indx]] # Now, # of categories instead
			splits[indx, 4L] <- ncat + 1L:nadd # rows to contain the splits

			## Next 4 lines can be done without a loop, but become indecipherable
			for (i in 1L:nadd) {
				newc[i, 1L:(cats[(cvar[indx])[i]])] <- -as.integer(cdir[i])
				newc[i, 1L:ccut[i]] <- as.integer(cdir[i])
			}
			catmat <- if (ncat == 0L) newc
				else {
					## newc may have more cols than existing categorical splits
					## the documentation says that levels which do no exist are '2'
					## and we add 2 later, so record as 0 here.
					cs <- ctfit$csplit
					ncs <- ncol(cs); ncc <- ncol(newc)
					if (ncs < ncc) cs <- cbind(cs, matrix(0L, nrow(cs), ncc - ncs))
					rbind(cs, newc)
				}
			ncat <- ncat + nadd
		} else catmat <- ctfit$csplit

		if (nsplit == 0L) {   
			frame <- data.frame(row.names = 1L,
								var = "<leaf>",
								n = ctfit$inode[, 5L],
								wt = ctfit$dnode[, 3L],
								dev = ctfit$dnode[, 1L],
								yval = ctfit$dnode[, 4L],
								complexity = ctfit$dnode[, 2L],
								ncompete = 0L,
								nsurrogate = 0L)
		} else {
			temp <- ifelse(index == 0L, 1L, index)
			svar <- ifelse(index == 0L, 0L, ctfit$isplit[temp, 1L]) # var number
			frame <- data.frame(row.names = ctfit$inode[, 1L],
								var = tname[svar + 1L],
								n = ctfit$inode[, 5L],
								wt = ctfit$dnode[, 3L],
								dev = ctfit$dnode[, 1L],
								yval = ctfit$dnode[, 4L],
								complexity = ctfit$dnode[, 2L],
								ncompete = pmax(0L, ctfit$inode[, 3L] - 1L),
								nsurrogate = ctfit$inode[, 4L])
		}


		if (is.null(init$summary))
			stop("Initialization routine is missing the 'summary' function")
		functions <- if (is.null(init$print)) list(summary = init$summary)
			else list(summary = init$summary, print = init$print)
		if (!is.null(init$text)) functions <- c(functions, list(text = init$text))
		if (method == "user") functions <- c(functions, mlist)

		where <- ctfit$which
		names(where) <- row.names(X)

		ans <- list(frame = frame,
					where = where,
					call = Call, terms = Terms,
					cptable = t(ctfit$cptable),
					method = method,
					control = controls,
					functions = functions,
					numresp = init$numresp)
		if (nsplit) ans$splits = splits
		if (ncat > 0L) ans$csplit <- catmat + 2L
		if (nsplit) ans$variable.importance <- importance(ans)
		if (y) ans$y <- Y
		if (x) {
			ans$x <- X
			ans$wt <- wt
		}
		ans$ordered <- isord
		if (!is.null(attr(m, "na.action"))) ans$na.action <- attr(m, "na.action")
		if (!is.null(xlevels)) attr(ans, "xlevels") <- xlevels
		if (method == "class") attr(ans, "ylevels") <- init$ylevels
		class(ans) <- "rpart"
		if(ncol(ans$cptable) >= 4) {
			ans$cptable[,4]  <- ans$cptable[,4] / ans$cptable[1, 4]
		}
		ans <- honest.est.causalTree(ans, est_X, est_wts, est_treatment, est_Y)
		#estimate honest causaltree with train X and compare with est.causaltree after pruning
		ans
}
