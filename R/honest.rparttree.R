#
#  The recursive partitioning function, for R
#
honest.rparttree <-
    function(formula, data, weights, subset, est_data,est_weights, na.action = na.rpart, method,
             model = FALSE, x = FALSE, y = TRUE, parms, control, cost, ...)
{
    Call <- match.call()
    if (is.data.frame(model)) {
        m <- model
        model <- FALSE
    } else {
        indx <- match(c("formula", "data", "weights", "subset"),
                      names(Call), nomatch = 0L)
        if (indx[1] == 0L) stop("a 'formula' argument is required")
        temp <- Call[c(1L, indx)]      # only keep the arguments we wanted
        temp$na.action <- na.action    # This one has a default
        temp[[1L]] <- quote(stats::model.frame) # change the function called
        m <- eval.parent(temp)
    }

    Terms <- attr(m, "terms")
    if (any(attr(Terms, "order") > 1L))
	stop("Trees cannot handle interaction terms")

    Y <- model.response(m)
    wt <- model.weights(m)
    if (any(wt < 0)) stop("negative weights not allowed")
    if (!length(wt)) wt <- rep(1, nrow(m))
    offset <- model.offset(m)
    # X <- rpart.matrix(m)
    X <- causalTree.matrix(m)
    nobs <- nrow(X)
    nvar <- ncol(X)
###

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
    
    # if (missing(est_treatment)) {
    #   stop("Not give the treatment status of honest estimation data set!\n ")
    # }
    # if (sum(est_treatment %in% c(0,1)) != est_nobs) {
    #   stop("The treatment status should be 1 or 0 only: 1 represent treated and 0 represent controlled.")
    # }
    # if (sum(est_treatment) == 0 || sum(est_treatment) == est_nobs) {
    #   stop("The data only contains treated cases or controlled cases, please check 'est_treatment' again.") 
    # }
    # 
    
    if (est_nvar != nvar) {
      stop("Honest estimation data set should have same variables as training data!\n")
    }
    
    print("done setting estimation data")
###
    if (missing(method)) {
	method <- if (is.factor(Y) || is.character(Y)) "class"
        else if (inherits(Y, "Surv")) "exp"
	else if (is.matrix(Y)) "poisson"
	else "anova"
    }

    if (is.list(method)) {
        ## User-written split methods
	mlist <- method
	method <- "user"

        ## Set up C callback.  Assign the result to a variable to avoid
        ## garbage collection
	init <- if (missing(parms)) mlist$init(Y, offset, wt = wt)
	        else mlist$init(Y, offset, parms, wt)
        keep <- rpartcallback(mlist, nobs, init)

	method.int <- 4L             # the fourth entry in func_table.h
#	numresp <- init$numresp
#	numy <- init$numy
	parms <- init$parms
    } else {
	method.int <- pmatch(method, c("anova", "poisson", "class", "exp"))
	if (is.na(method.int)) stop("Invalid method")
	method <- c("anova", "poisson", "class", "exp")[method.int]
	if (method.int == 4L) method.int <- 2L

        ## If this function is being retrieved from the rpart package, then
        ##   preferentially "get" the init function from there.  But don't
        ##   lock in the rpart package otherwise, so that we can still do
        ##   standalone debugging.
	init <- if (missing(parms))
            get(paste("causalTree", method, sep = "."),
	              envir = environment())(Y, offset)
	            # envir = environment())(Y, offset, , wt)
        else
            get(paste("causalTree", method, sep = "."),
                envir = environment())(Y, offset, parms)
                # envir = environment())(Y, offset, parms, wt)
	            
        ## avoid saving environment on fitted objects
        ns <- asNamespace("causalTree")
        if (!is.null(init$print)) environment(init$print) <- ns
        if (!is.null(init$summary)) environment(init$summary) <- ns
        if (!is.null(init$text)) environment(init$text) <- ns
    }

    Y <- init$y

    xlevels <- .getXlevels(Terms, m)
    cats <- rep(0L, ncol(X))
    if (!is.null(xlevels))
	cats[match(names(xlevels), colnames(X))] <-
            unlist(lapply(xlevels, length))

    ## We want to pass any ... args to rpart.control, but not pass things
    ##  like "dats = mydata" where someone just made a typo.  The use of ...
    ##  is simply to allow things like "cp = 0.05" with easier typing
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
        ## make random groups
        xgroups <- sample(rep(1L:xval, length = nobs), nobs, replace = FALSE)
    } else if (length(xval) == nobs) {
	xgroups <- xval
	xval <- length(unique(xgroups))
    } else {
        ## Check to see if observations were removed due to missing
	if (!is.null(attr(m, "na.action"))) {
            ## if na.rpart was used, then na.action will be a vector
	    temp <- as.integer(attr(m, "na.action"))
	    xval <- xval[-temp]
	    if (length(xval) == nobs) {
		xgroups <- xval
		xval <- length(unique(xgroups))
            } else stop("Wrong length for 'xval'")
        } else stop("Wrong length for 'xval'")
    }

    ##
    ## Incorporate costs
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
    temp <- as.double(unlist(init$parms))
    if (!length(temp)) temp <- 0    # if parms is NULL pass a dummy
    rpfit <- .Call(C_rpart,
                   ncat = as.integer(cats * !isord),
                   method = as.integer(method.int),
                   as.double(unlist(controls)),
                   temp,
                   as.integer(xval),
                   as.integer(xgroups),
                   as.double(t(init$y)),
                   X,
                   wt,
                   as.integer(init$numy),
                   as.double(cost))

    nsplit <- nrow(rpfit$isplit) # total number of splits, primary and surrogate
    ## total number of categorical splits
    ncat <- if (!is.null(rpfit$csplit)) nrow(rpfit$csplit) else 0L
#    nodes <- nrow(rpfit$inode)
    if (nsplit == 0L) xval <- 0L # No xvals were done if no splits were found

    numcp <- ncol(rpfit$cptable)
    temp <- if (nrow(rpfit$cptable) == 3L) c("CP", "nsplit", "rel error")
            else c("CP", "nsplit", "rel error", "xerror", "xstd")
    dimnames(rpfit$cptable) <- list(temp, 1L:numcp)

    tname <- c("<leaf>", colnames(X))
    splits <- matrix(c(rpfit$isplit[, 2:3], rpfit$dsplit), ncol = 5L,
                     dimnames = list(tname[rpfit$isplit[, 1L] + 1L],
                     c("count", "ncat", "improve", "index", "adj")))
    index <- rpfit$inode[, 2L]  # points to the first split for each node

    ## Now, make ordered factors look like factors again (a printout choice)
    nadd <- sum(isord[rpfit$isplit[, 1L]])
    if (nadd > 0L) { # number of splits at an ordered factor.
	newc <- matrix(0L, nadd, max(cats))
	cvar <- rpfit$isplit[, 1L]
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
            cs <- rpfit$csplit
            ncs <- ncol(cs); ncc <- ncol(newc)
            if (ncs < ncc) cs <- cbind(cs, matrix(0L, nrow(cs), ncc - ncs))
            rbind(cs, newc)
        }
	ncat <- ncat + nadd
    } else catmat <- rpfit$csplit

    ## NB: package adabag depends on 'var' being a factor.
    if (nsplit == 0L) {                    # tree with no splits
	frame <- data.frame(row.names = 1L,
			    var = "<leaf>",
			    n = rpfit$inode[, 5L],
			    wt = rpfit$dnode[, 3L],
			    dev = rpfit$dnode[, 1L],
			    yval = rpfit$dnode[, 4L],
			    complexity = rpfit$dnode[, 2L],
			    ncompete = 0L,
			    nsurrogate = 0L)
    } else {
	temp <- ifelse(index == 0L, 1L, index)
	svar <- ifelse(index == 0L, 0L, rpfit$isplit[temp, 1L]) # var number
	frame <- data.frame(row.names = rpfit$inode[, 1L],
                            ## maybe better to specify tname as the level?
			    var = tname[svar + 1L],
			    n = rpfit$inode[, 5L],
			    wt = rpfit$dnode[, 3L],
			    dev = rpfit$dnode[, 1L],
			    yval = rpfit$dnode[, 4L],
			    complexity = rpfit$dnode[, 2L],
			    ncompete = pmax(0L, rpfit$inode[, 3L] - 1L),
			    nsurrogate = rpfit$inode[, 4L])
    }
    if (method.int == 3L) {
        ## Create the class probability vector from the class counts, and
        ##   add it to the results
        ## Also scale the P(T) result
        ## The "pmax" 3 lines down is for the case of a factor y which has
        ##   no one at all in one of its classes.  Both the prior and the
        ##   count will be zero, which led to a 0/0.
        numclass <- init$numresp - 2L
        nodeprob <- rpfit$dnode[, numclass + 5L] / sum(wt) # see ginidev.c
        temp <- pmax(1L, init$counts)   # overall class freq in data
        temp <- rpfit$dnode[, 4L + (1L:numclass)] %*% diag(init$parms$prior/temp)
        yprob <- temp /rowSums(temp)    # necessary with altered priors
        yval2 <- matrix(rpfit$dnode[, 4L + (0L:numclass)], ncol = numclass + 1L)
	frame$yval2 <- cbind(yval2, yprob, nodeprob)
    } else if (init$numresp > 1L)
        frame$yval2 <- rpfit$dnode[, -(1L:3L), drop = FALSE]

    if (is.null(init$summary))
        stop("Initialization routine is missing the 'summary' function")
    functions <- if (is.null(init$print)) list(summary = init$summary)
                 else list(summary = init$summary, print = init$print)
    if (!is.null(init$text)) functions <- c(functions, list(text = init$text))
    if (method == "user") functions <- c(functions, mlist)

    where <- rpfit$which
    names(where) <- row.names(m)
    print("end of honest.rpart tree growing,before ans")
    ans <- list(frame = frame,
             where = where,
             call = Call, terms = Terms,
             cptable = t(rpfit$cptable),
             method = method,
             parms = init$parms,
             control = controls,
             functions = functions,
             numresp = init$numresp)
    if (nsplit) ans$splits = splits
    if (ncat > 0L) ans$csplit <- catmat + 2L
    if (nsplit) ans$variable.importance <- importance(ans)
    if (model) {
	ans$model <- m
	if (missing(y)) y <- FALSE
    }
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
    #save 
    save(list="ans","est_X","est_wts","est_Y",file="test1_honest_rpart_before_est_call.rData")
    ans <- honest.est.rparttree(ans, est_X, est_wts, est_Y)
    
    ans
    
}
