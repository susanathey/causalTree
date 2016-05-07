#
# Run down the built tree and get the final leaf ids for estimation sample
#
est.causalTree <- function(fit, x)
{
    frame <- fit$frame
    if (nrow(frame) == 1L)                # root only tree
        return(structure(rep(1, nrow(x), names = rownames(x))))

    nc <- frame[, c("ncompete", "nsurrogate")]
    frame$index <- 1L + c(0L, cumsum((frame$var != "<leaf>") +
                                     nc[[1L]] + nc[[2L]]))[-(nrow(frame) + 1L)]
    frame$index[frame$var == "<leaf>"] <- 0L
    vnum <- match(rownames(fit$split), colnames(x))
    if (any(is.na(vnum)))
        stop("Tree has variables not found in new data")
    temp <- .Call(C_estimate_causalTree,
                  as.integer(dim(x)),
                  as.integer(dim(frame)[1L]),
                  as.integer(dim(fit$splits)),
                  as.integer(if (is.null(fit$csplit)) rep(0L, 2L) else dim(fit$csplit)),
                  as.integer(row.names(frame)),
                  as.integer(unlist(frame[, c("n", "ncompete", "nsurrogate", "index")])),
                  as.integer(vnum),
                  as.double(fit$splits),
                  as.integer(fit$csplit - 2L),
                  as.integer((fit$control)$usesurrogate),
                  as.double(x),
                  as.integer(is.na(x)))
    names(temp) <- rownames(x)
    temp
}
