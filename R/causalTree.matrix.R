causalTree.matrix <- function(frame)
{

    if (!inherits(frame, "data.frame") ||
       is.null(attr(frame, "terms")))  return(as.matrix(frame))

    for (i in 1:ncol(frame)) {
        if (is.character(frame[[i]])) frame[[i]] <- as.numeric(factor(frame[[i]]))
        else if (!is.numeric(frame[[i]])) frame[[i]] <- as.numeric(frame[[i]])
    }

    X <- model.matrix(attr(frame, "terms"), frame)[, -1L, drop = FALSE]

    colnames(X) <- sub("^`(.*)`", "\\1", colnames(X))
    class(X) <- c("causalTree.matrix", class(X)) 
    X
}

