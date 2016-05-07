causalTree.anova <- function(y, offset, wt)
{
    if (!is.null(offset)) y <- y - offset
      list(y = y, numresp = 1L, numy = 1L,
	 summary = function(yval, dev, wt, ylevel, digits) {
     
	   paste0("  causal effect=", formatg(yval, digits),
	          ", error=" , formatg(dev, digits))
         },
	 text = function(yval, dev, wt, ylevel, digits, n, use.n ) {
	     if (use.n) paste0(formatg(yval, digits), "\nn=", n) else
             formatg(yval, digits)
         })
}
