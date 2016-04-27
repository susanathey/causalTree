##
## Do Rpart estimations given a tree and a matrix of predictors
honest.est.causalTree <- function(fit, x, wt, treatment, y)
{
    frame <- fit$frame
    if (nrow(frame) == 1L) {              # root only tree
        #print ("here end?")
        #return(structure(rep(1, nrow(x), names = rownames(x))))
    }
    
    nc <- frame[, c("ncompete", "nsurrogate")]
    frame$index <- 1L + c(0L, cumsum((frame$var != "<leaf>") +
                                         nc[[1L]] + nc[[2L]]))[-(nrow(frame) + 1L)]
    frame$index[frame$var == "<leaf>"] <- 0L
    vnum <- match(rownames(fit$split), colnames(x))
    if (any(is.na(vnum)))
        stop("Tree has variables not found in new data")
    temp <- .Call(C_honest_estimate_causalTree,
                  as.integer(dim(x)),
                  as.integer(dim(frame)[1L]),
                  as.integer(dim(fit$splits)),
                  as.integer(if (is.null(fit$csplit)) rep(0L, 2L) else dim(fit$csplit)),
                  as.integer(row.names(frame)),
                  as.integer(unlist(frame[, c("ncompete", "nsurrogate", "index")])),
                  as.integer(frame[, "n"]),
                  as.double(frame[, "wt"]), 
                  as.double(frame[, "dev"]),
                  as.double(frame[, "yval"]),
                  as.integer(vnum),
                  as.double(fit$splits),
                  as.integer(fit$csplit - 2L), # we added 2 in causalTree.R
                  as.integer((fit$control)$usesurrogate),
                  as.double(x),
                  as.double(wt),
                  as.double(treatment),
                  as.double(y),
                  as.integer(is.na(x)))
    return (fit)
    #names(temp) <- rownames(x)
    #temp
}
