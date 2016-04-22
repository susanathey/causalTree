#
# Caclulate variable importance
# Each primary split is credited with the value of splits$improve
# Each surrogate split gets split$adj times the primary split's value
#
# Called only internally by causalTree
#
importance <- function(fit)
{
    ff <- fit$frame
    fpri <- which(ff$var != "<leaf>")  # points to primary splits in ff
    spri <- 1 + cumsum(c(0, 1 + ff$ncompete[fpri] + ff$nsurrogate[fpri]))
    spri <- spri[seq_along(fpri)] # points to primaries in the splits matrix
    nsurr <- ff$nsurrogate[fpri]  # number of surrogates each has

    sname <- vector("list", length(fpri))
    sval <- sname

    ## The importance for primary splits needs to be scaled
    ## It was a printout choice for the anova method to list % improvement in
    ##  the sum of squares, an importance calculation needs the total SS.
    ## All the other methods report an unscaled change.
     scaled.imp <- if (fit$method == "anova")
        fit$splits[spri, "improve"] * ff$dev[fpri]
    else fit$splits[spri, "improve"]

    sdim <- rownames(fit$splits)
    for (i in seq_along(fpri)) {
        ## points to surrogates
        if (nsurr[i] > 0L) {
            indx <- spri[i] + ff$ncompete[fpri[i]] + seq_len(nsurr[i])
            sname[[i]] <- sdim[indx]
            sval[[i]] <- scaled.imp[i] * fit$splits[indx, "adj"]
        }
    }

    import <- tapply(c(scaled.imp, unlist(sval)),
                     c(as.character(ff$var[fpri]), unlist(sname)),
                     sum)
    sort(c(import), decreasing = TRUE) # a named vector
}
