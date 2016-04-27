## rescaled exponential splitting
##  The survival object 'y' is rescaled so that
##    a. overall death rate is 1.0
##    b. death rate within small intervals of time is 1
##  The first makes printouts easier to read, since the rates in subnodes are
##  numbers like "1.23" (23% higher event rate than the root node) instead
##  of "0.0014" (which requires looking back at the root node rate and
##  dividing).  The second makes the data appear to be Poisson, and causes
##  the early splits at least to be equivalent to the local full likelihood
##  method of LeBlanc and Crowley
##
causalTree.exp <- function(y, offset, parms, wt)
{
    if (!inherits(y, "Surv"))
        stop("Response must be a 'survival' object - use the 'Surv()' function")

    ny <- ncol(y)
    n <- nrow(y)

    status <- y[, ny]
    if (any(y[, 1L] <= 0)) stop("Observation time must be > 0")
    if (all(status == 0)) stop("No deaths in data set")
    time <- y[ , ny - 1L]

    ## Make my table of time intervals.  The first goes from 0 to the first
    ##   death, the next from death 2 to death 3, ..., and the last from
    ##   "next to last death" to "max time in dataset".
    ## We also need to avoid a pathological condition in some data sets, where
    ##   two death times differ by a trivial amount, e.g., 10^-16, perhaps due
    ##   to roundoff error in creating the input data.  Ammalgamate such
    ##   intervals.  This turns out to be hard to do in S, but easy in C
    dtimes <- sort(unique(time[status == 1])) # unique death times
    temp <- .Call(C_causalTreeexp2, as.double(dtimes), as.double(.Machine$double.eps))
    dtimes <- dtimes[temp == 1]

    ## For the sake of speed, restrict the number of intervals to be < 1000.
    ##   (Actually, anything > 100 is probably overkill for the
    ##   actual task at hand, which is to approximately scale to exponential).
    if (length(dtimes) > 1000) dtimes <- quantile(dtimes, 0:1000/1000)

    ## The last interval goes to the max time in the data set
    itable <- c(0, dtimes[-length(dtimes)], max(time)) # set of intervals


    drate2 <- function(n, ny, y, wt, itable) {
	## An alternative to the drate1 function
	## Why?  The pyears2 routine changed in 6/2001, with the inclusion
	##  of case weights.  We need the newer version.  If you have the
	##  older version of the survival library, the above will crash S.
	## How to tell -- list the pyears function, and see whether it's
	##  call to pyears2 has weights in the argument list.
	##
	time <- y[, ny - 1L]
	status <- y[, ny]
	ilength <- diff(itable)         #lengths of intervals
	ngct <- length(ilength)         #number of intervals

	## The code below is as opaque as any I've written, all in the
	##  service of "no for loops".
	## First, 'index' gives the time interval (as defined by itable)
	##  in which the end of each observation's follow-up (time) lies.
	##  Then 'itime' will be the amount of time spent in that last
	##  interval, which is of course somewhat less than ilength.
	index <- unclass(cut(time, itable, include.lowest = TRUE))
	itime <- time - itable[index]
	if (ny == 3L) {
	    ## there is both a start time and a stop time
	    ##  compute the amount of time NOT spent in the interval that
	    ##  the start time lies in.
	    stime <- y[, 1L]             #start time for each interval
	    index2 <- unclass(cut(stime, itable, include.lowest = TRUE))
	    itime2 <- stime - itable[index2]
        }

	## Compute the amount of person-years in each of the intervals
	##   This is:  (width of interval) * (number of "time" elements that
	##                                     end in an interval farther right)
	##            + (ending times in this interval)
	## By construction, I know that there is at least 1 obs in each of the
	##  intervals, so tab1 is of a determined length
	tab1 <- table(index)
	temp <- rev(cumsum(rev(tab1))) # cumsum, counting from the right
	pyears <- ilength * c(temp[-1L], 0) + tapply(itime, index, sum)
	if (ny == 3L) {
	    ## subtract off the time before "start"
	    tab2 <- table(index2, levels = 1:ngct) # force the length of tab2
	    temp <- rev(cumsum(rev(tab2)))
	    py2 <- ilength * c(0, temp[-ngct]) + tapply(itime2, index2, sum)
	    pyears <- pyears - py2
        }

	deaths <- tapply(status, index, sum)
	rate <- deaths/pyears           # hazard rate in each interval
	rate
    }

    ##
    ## Now, compute the "new y" for each observation.
    ##  This is a stretching of the time axis
    ## The cumulative hazard over each interval is rate*length(interval),
    ##  and is the basis of the rescaling.
    rate <- drate2(n, ny, y, wt, itable)
    cumhaz <- cumsum(c(0, rate * diff(itable)))
    newy <- approx(itable, cumhaz, time)$y
    if (ny == 3L) newy <- newy - approx(itable, cumhaz, y[, 1L])$y

    if (length(offset) == n)  newy <- newy * exp(offset)

    if (missing(parms)) parms <- c(shrink = 1L, method = 1L)
    else {
	parms <- as.list(parms)
        if (is.null(names(parms))) stop("You must input a named list for parms")
        parmsNames <- c("method", "shrink")
        indx <- pmatch(names(parms), parmsNames, 0L)
        if (any(indx == 0L))
            stop(gettextf("'parms' component not matched: %s",
                          names(parms)[indx == 0L]), domain = NA)
	else names(parms) <- parmsNames[indx]

	if (is.null(parms$method)) method <- 1L
	else method <- pmatch(parms$method, c("deviance", "sqrt"))
	if (is.na(method)) stop("Invalid error method for Poisson")

	if (is.null(parms$shrink)) shrink <- 2L - method
	else shrink <- parms$shrink
	if (!is.numeric(shrink) || shrink < 0L)
            stop("Invalid shrinkage value")
        parms <- c(shrink = shrink, method = method)
    }
    list(y = cbind(newy, y[, 2L]), parms = parms, numresp = 2L, numy = 2L,
	 summary = function(yval, dev, wt, ylevel, digits) {
	     paste0("  events=", formatg(yval[, 2L]),
                   ",  estimated rate=" , formatg(yval[, 1L], digits),
                   " , mean deviance=", formatg(dev/wt, digits))
         },
	 text = function(yval, dev, wt, ylevel, digits, n, use.n) {
	     if (use.n) paste0(formatg(yval[, 1L], digits), "\n",
                              formatg(yval[, 2L]) , "/", n)
             else paste(formatg(yval[, 1L], digits))
         })
}
