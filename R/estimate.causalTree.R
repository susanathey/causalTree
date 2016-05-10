estimate.causalTree <- function(object, data, weights, treatment, na.action = na.causalTree)
{
    if (!inherits(object, "rpart")) stop("Not a legitimate \"rpart\" object")
    # get the leaf of the object
    leaf <- as.numeric(row.names(object$frame)[which(object$frame$var == "<leaf>")])
    Terms <- object$terms
    data$tr <- treatment
    if (missing(weights)) {
        Terms <- object$terms
        m <- model.frame(Terms, data = data, na.action = na.action, treatment = tr,
                         xlev = attr(object, "xlevels"))
    } else {
        attr(Terms, "dataClasses")["weights"] <- "numeric"
        data$w <- weights
        m <- model.frame(Terms, data = data, na.action = na.action, treatment = tr,
                         weights = w, xlev = attr(object, "xlevels"))
    }
    
    if (!is.null(cl <- attr(Terms, "dataClasses")))
        .checkMFClasses(cl, m, TRUE)
    
    treatment <- m$`(treatment)`
    n <- nrow(m)
    Y <- model.response(m)
    X <- causalTree.matrix(m)
    if (missing(weights))
        wts <- rep(1, nrow(m))
    else 
        wts <- model.weights(m)
    new_object <- data.table::copy(object)
    ans <- honest.est.causalTree(new_object, X, wts, treatment, Y)
    return(ans)
}
