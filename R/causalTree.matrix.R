## Rewritten 11/2012; the prior version of the code failed with the
##   formula y ~ . - zed;  zed got included in the result
## Rather than be more sophisticated in how I understand the terms
##  structure (because there might be yet another special case we've not
##  thought of) use model.matrix on modified data.
##
causalTree.matrix <- function(frame)
{
    ## First line is just a failsafe: this should always be called with
    ##   a model frame.
    if (!inherits(frame, "data.frame") ||
       is.null(attr(frame, "terms")))  return(as.matrix(frame))

    ## turn other classes into numerics.
    ## We replace columns in frame rather than making a new object, since
    ## model.matrix wants a model.frame object as it's argument.
    for (i in 1:ncol(frame)) {
        if (is.character(frame[[i]])) frame[[i]] <- as.numeric(factor(frame[[i]]))
        else if (!is.numeric(frame[[i]])) frame[[i]] <- as.numeric(frame[[i]])
    }
    ## Toss the intercept term when done (column 1)
    X <- model.matrix(attr(frame, "terms"), frame)[, -1L, drop = FALSE]
    ## model.matrix labels columns with backticks, and causalTree.matrix did not.
    colnames(X) <- sub("^`(.*)`", "\\1", colnames(X))
    class(X) <- c("causalTree.matrix", class(X)) # ipred package expects this class
    X
}

