#
#  Plot best fit causal tree and p-values at leaf nodes.
#
# 

# Helper functions for plot.causalTree

# Concatenates p-values to leaf node labels
node.pvals <- function(x, labs, digits, varlen) {
  ifelse(is.na(x$frame$p.value), labs,
         paste(labs, "\np=", round(x$frame$p.value, 4)))
}

# Uses the model and response (x and y) values from a rpart object to construct
# a data frame from which p-values between control and treatment within a 
# leaf node is calculated. TODO: multiple testing corrections?
add.pvals <- function(tree) {
  dat <- cbind(treatment=tree$x[,1], outcome=tree$y, node=tree$where)
  dat <- aggregate(outcome ~ treatment + node, data=dat, FUN=c)
  merged <- merge(dat[dat$treatment == 0,], dat[dat$treatment == 1,], by="node",
                  suffixes=c(".ctl", ".trt"))
  p.values <- do.call(
    rbind, apply(merged, 1, FUN=function(x) {
      data.frame(node=x$node,
                 p.value=t.test(x$outcome.ctl, x$outcome.trt)$p.value) }))
  tree$frame$p.value[p.values$node] <- p.values$p.value
  return(tree)
}

# Takes the optimally pruned causal tree and adds pvalues.  Then plots.
plot.causalTree <- function(tree, ...) {
  if (is.null(tree$x) || is.null(tree$y))
    stop("Must build causalTree with x=TRUE, y=TRUE")
  opCp <- tree$cptable[,1][which.min(tree$cptable[,4])]
  opFit <- prune(tree, opCp)
  opFit <- add.pvals(opFit) 
  rpart.plot(opFit, node.fun=node.pvals, ...)
}
