print.rpart <- function(x, minlength = 0L, spaces = 2L, cp,
               digits = getOption("digits"), ...)
{
    if (!inherits(x, "rpart")) stop("Not a legitimate \"rpart\" object")

    if (!missing(cp)) x <- prune.rpart(x, cp = cp)
    frame <- x$frame
    ylevel <- attr(x, "ylevels")
    node <- as.numeric(row.names(frame))
    depth <- tree.depth(node)
    indent <- paste(rep(" ", spaces * 32L), collapse = "")
    ## 32 is the maximal depth
    indent <- if (length(node) > 1L) {
        indent <- substring(indent, 1L, spaces * seq(depth))
        paste0(c("", indent[depth]), format(node), ")")
    } else paste0(format(node), ")")

    tfun <- (x$functions)$print
    yval <- if (!is.null(tfun)) {
	if (is.null(frame$yval2)) tfun(frame$yval, ylevel, digits)
	else tfun(frame$yval2, ylevel, digits)
    } else format(signif(frame$yval, digits))
    term <- rep(" ", length(depth))
    term[frame$var == "<leaf>"] <- "*"
    z <- labels(x, digits = digits, minlength = minlength, ...)
    n <- frame$n
    z <- paste(indent, z, n, format(signif(frame$dev, digits)), yval, term)

    omit <- x$na.action
    if (length(omit)) cat("n=", n[1L], " (", naprint(omit), ")\n\n", sep = "")
    else cat("n=", n[1L], "\n\n")

    ## This is stolen, unabashedly, from print.tree
    if (x$method == "class") cat("node), split, n, loss, yval, (yprob)\n")
    else cat("node), split, n, deviance, yval\n")
    cat("      * denotes terminal node\n\n")

    cat(z, sep = "\n")
    invisible(x)
}
