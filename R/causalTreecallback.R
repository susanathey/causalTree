## This routine sets up the callback code for user-written split
##  routines in causalTree
##
causalTreecallback <- function(mlist, nobs, init)
{
    if (length(mlist) < 3L)
        stop("User written methods must have 3 functions")
    if (!is.function(mlist$init))
        stop("User written method does not contain an 'init' function")
    if (!is.function(mlist$split))
        stop("User written method does not contain a 'split' function")
    if (!is.function(mlist$eval))
        stop("User written method does not contain an 'eval' function")

    user.eval <- mlist$eval
    user.split <- mlist$split

    numresp <- init$numresp
    numy <- init$numy
    parms <- init$parms

    ##
    ## expr2 is an expression that will call the user "evaluation"
    ##   function, and check that what comes back is valid
    ## expr1 does the same for the user "split" function
    ##
    ## For speed in the C interface, yback, xback, and wback are
    ##  fixed S vectors of a fixed size, and nback tells us how
    ##  much of the vector is actually being used on this particular
    ##  callback.
    ##
    if (numy == 1L) {
        expr2 <- quote({
            temp <- user.eval(yback[1:nback], wback[1:nback], parms)
            if (length(temp$label) != numresp)
                stop("User 'eval' function returned invalid label")
            if (length(temp$deviance) != 1L)
                stop("User 'eval' function returned invalid deviance")
            as.numeric(as.vector(c(temp$deviance, temp$label)))
        })
        expr1 <- quote({
            if (nback < 0L) { #categorical variable
                n2 <- -nback
                temp <- user.split(yback[1L:n2], wback[1L:n2],
                                   xback[1L:n2], parms, FALSE)
                ncat <- length(unique(xback[1L:n2]))
                if (length(temp$goodness) != ncat - 1L ||
                    length(temp$direction) != ncat)
                    stop("Invalid return from categorical 'split' function")
            } else {
                temp <- user.split(yback[1L:nback], wback[1L:nback],
                                   xback[1L:nback], parms, TRUE)
                if (length(temp$goodness) != (nback - 1L))
                    stop("User 'split' function returned invalid goodness")
                if (length(temp$direction) != (nback - 1L))
                    stop("User 'split' function returned invalid direction")
            }
            as.numeric(as.vector(c(temp$goodness, temp$direction)))
        })
    } else {
        expr2 <- quote({
            tempy <- matrix(yback[1L:(nback * numy)], ncol = numy)
            temp <- user.eval(tempy, wback[1L:nback], parms)
            if (length(temp$label) != numresp)
                stop("User 'eval' function returned invalid label")
            if (length(temp$deviance) != 1L)
                stop("User 'eval' function returned invalid deviance")
            as.numeric(as.vector(c(temp$deviance, temp$label)))
        })
        expr1 <- quote({
            if (nback < 0L) { #categorical variable
                n2 <- -nback
                tempy <- matrix(yback[1L:(n2 * numy)], ncol = numy)
                temp <- user.split(tempy, wback[1L:n2], xback[1L:n2],
                                   parms, FALSE)
                ncat <- length(unique(xback[1L:n2]))
                if (length(temp$goodness) != ncat - 1L ||
                    length(temp$direction) != ncat)
                    stop("Invalid return from categorical 'split' function")
            } else {
                tempy <- matrix(yback[1L:(nback * numy)], ncol = numy)
                temp <- user.split(tempy, wback[1:nback], xback[1L:nback],
                                   parms, TRUE)
                if (length(temp$goodness) != (nback - 1L))
                    stop("User 'split' function returned invalid goodness")
                if (length(temp$direction) != (nback - 1L))
                    stop("User 'split' function returned invalid direction")
            }
            as.numeric(as.vector(c(temp$goodness, temp$direction)))
        })
    }
    ##
    ##  The vectors nback, wback, xback and yback will have their
    ##  contents constantly re-inserted by C code.  It's one way to make
    ##  things very fast.  It is dangerous to do this, so they
    ##  are tossed into a separate frame to isolate them.  Evaluations of
    ##  the above expressions occur in that frame.
    ##
    rho <- new.env()
    assign("nback", integer(1L), envir = rho)
    assign("wback", double(nobs), envir = rho)
    assign("xback", double(nobs), envir = rho)
    assign("yback", double(nobs * numy), envir = rho)
    assign("user.eval", user.eval, envir = rho)
    assign("user.split", user.split, envir = rho)
    assign("numy", numy, envir = rho)
    assign("numresp", numresp, envir = rho)
    assign("parms", parms, envir = rho)
    .Call(C_init_ctcallback, rho, as.integer(numy), as.integer(numresp),
          expr1, expr2)
    list(expr1 = expr1, expr2 = expr2, rho = rho)
}
