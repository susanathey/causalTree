causalTree.control <-
    function(minsplit = 20L, minbucket = round(minsplit/3), cp = 0,
             maxcompete = 4L, maxsurrogate = 5L, usesurrogate = 2L, xval = 10L,
             surrogatestyle = 0L, maxdepth = 30L, ...)
    {
        
        if (maxcompete < 0L) {
            warning("The value of 'maxcompete' supplied is < 0; the value 0 was used instead")
            maxcompete <- 0L
        }
        if (any(xval < 0L)) {
            warning("The value of 'xval' supplied is < 0; the value 0 was used instead")
            xval <- 0L
        }
        if (maxdepth > 30L) stop("Maximum depth is 30")
        if (maxdepth < 1L)  stop("Maximum depth must be at least 1")
        if (missing(minsplit) && !missing(minbucket)) minsplit <- minbucket * 3L
        if ((usesurrogate < 0L) || (usesurrogate > 2L)) {
            warning("The value of 'usesurrogate' supplied was out of range, the default value of 2 is used instead.")
            usesurrogate <- 2L
        }
        if ((surrogatestyle < 0L) || (surrogatestyle > 1L)) {
            warning("The value of 'surrogatestyle' supplied was out of range, the default value of 0 is used instead.")
            surrogatestyle <- 0L
        }
        
        ## Because xval can be of length either 1 or n, and the C code
        ##   refers to parameters by number, i.e., "opt[5]" in rpart.c,
        ##   the xval parameter should always be last on the list.
        list(minsplit = minsplit, minbucket = minbucket, cp = cp,
             maxcompete = maxcompete, maxsurrogate = maxsurrogate,
             usesurrogate = usesurrogate,
             surrogatestyle = surrogatestyle, maxdepth = maxdepth, xval = xval)
    }
