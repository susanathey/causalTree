# get model frame of causalTree, same as rpart
model.frame.causalTree <- function(formula, ...)
{
    m <- formula$model
    if (!is.null(m)) return(m)
    oc <- formula$call
    if (substring(deparse(oc[[1L]]), 1L, 7L) == "predict") {
        m <- eval(oc$newdata)
        if (is.null(attr(m, "terms"))) {
            object <- eval(oc$object)
            m <- model.frame(object$terms, m, na.causalTree)
        }
        return(m)
    }
    while(!deparse(oc[[1L]]) %in%  c("causalTree", "causalTree::causalTree", "causalTree:::causalTree"))
        oc <- eval(oc[[2L]])$call
    oc$subset <- names(formula$where)
    oc$method <- formula$method
    eval(oc)
}

