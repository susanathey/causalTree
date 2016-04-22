/*
 *   Cut down a list of death times, to avoid ones that differ by
 *    very, very little.
 *	n	number of death times
 *	y	list of death times, sorted
 *	eps     machine precision
 * output
 *      keep    1=keep this one, 0=don't
 */
#include "causalTree.h"
#include <Rinternals.h>

static void
Rpartexp2(int n, double *y, double eps, int *keep)
{
    double delta;
    int i, j;
    double lasty;

    /* let delta = eps * interquartile range */

    i = n / 4;
    j = (3 * n) / 4;
    delta = eps * (y[j] - y[i]);


    /*
     * make sure that each y is at least "delta" greater than
     * the last y that we have decided to keep
     */
    lasty = y[0];
    keep[0] = 1;
    for (i = 1; i < n; i++) {
	if ((y[i] - lasty) <= delta)
	    keep[i] = 0;
	else {
	    keep[i] = 1;
	    lasty = y[i];
	}
    }
}

SEXP
causalTreeexp2(SEXP dtimes, SEXP eps)
{
    int n = LENGTH(dtimes);
    SEXP keep = PROTECT(allocVector(INTSXP, n));
    Rpartexp2(n, REAL(dtimes), asReal(eps), INTEGER(keep));
    UNPROTECT(1);
    return keep;
}
