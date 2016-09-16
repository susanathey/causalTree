/*
 *  Do rpart predictions given the matrix form of the tree.
 *
 *  Input
 *      dimx        : # of rows and columns in the new data
 *      nnode       : # of nodes in the tree
 *      nsplit      : # of split structures
 *      dimc        : dimension of the categorical splits matrix
 *      nnum        : node number for each row of 'nodes'
 *      nodes       : matrix of node info
 *                      row 0= count, 1=index of primary, 2=#competitors,
 *                          3= number of surrogates
 *      vnum        : variable number of each split
 *      split       : matrix of split info
 *                  :   row 0=usage count, 1= #categories if >1, otherwise
 *                         the split parity, 2= utility, 3= index to csplit
 *                         or numeric split point
 *      csplit      : matrix of categorical split info
 *      usesur      : at what level to use surrogates
 *      xdata       : the new data
 *      xmiss       : shows missings in the new data
 *
 *  Output
 *      where       : the "final" row in nodes for each observation
 */
//#include "rpart.h"
//#include "rpartproto.h"

#include "causalTree.h"
#include "causalTreeproto.h"

static void
honest_estimate_rparttree0(const int *dimx, int nnode, int nsplit, const int *dimc,
	    const int *nnum, const int *nodes2, const int *vnum,
	    const double *split2, const int *csplit2, const int *usesur,
	    int *n1, double *wt1, double *dev1, double *yval1,
	    const double *xdata2,const double *wt2,const double *y2, const int *xmiss2, int *where)
{
    int i, j;
    int n;
    int ncat;
    int node, nspl, var, dir;
    int lcount, rcount;
    int npos;
    double temp;
    const int *nodes[4];
    const double *split[4];
    const int **csplit = NULL, **xmiss;
    const double **xdata;
///
Rprintf("inside honest C code fn: 1\n");
double *trs = NULL;
double *cons = NULL; 
double *trsums = NULL; 
double *consums = NULL;
double *trsqrsums = NULL;
double *consqrsums = NULL;
int nnodemax = -1;
int *invertdx = NULL;

trs = (double *) ALLOC(nnode, sizeof(double));
cons = (double *) ALLOC(nnode, sizeof(double));
trsums = (double *) ALLOC(nnode, sizeof(double));
consums = (double *) ALLOC(nnode, sizeof(double));
trsqrsums = (double *) ALLOC(nnode, sizeof(double));
consqrsums = (double *) ALLOC(nnode, sizeof(double));
Rprintf("inside honest C code fn: 2\n");

// initialize:
for (i = 0; i < nnode; i++) {
  trs[i] = 0.;
  cons[i] = 0.;
  trsums[i] = 0.;
  consums[i] = 0.;
  trsqrsums[i] = 0.;
  consqrsums[i] = 0.;
  n1[i] = 0;
  wt1[i] = 0.;
  if (nnum[i] > nnodemax) {
    nnodemax = nnum[i]; 
  }
}
Rprintf("inside honest C code fn: 3\n");

invertdx = (int *) ALLOC(nnodemax + 1, sizeof(int));
// construct an invert index:
for (i = 0; i <= nnodemax + 1; i++) {
  invertdx[i] = -1;
}
for (i = 0; i != nnode; i++) {
  invertdx[nnum[i]] = i;
}
Rprintf("inside honest C code fn: 4\n");


/*n = dimx[0]; // n = # of obs
for (i = 0; i < 3; i++) {
  nodes[i] = &(nodes2[nnode * i]);
}*/

///
    n = dimx[0];
    for (i = 0; i < 4; i++) {
	nodes[i] = &(nodes2[nnode * i]); //tbd: change this to above from causaltree?
	split[i] = &(split2[nsplit * i]);
    }
    Rprintf("inside honest C code fn: 5\n");
    
    if (dimc[1] > 0) {
	csplit = (const int **) ALLOC((int) dimc[1], sizeof(int *));
	for (i = 0; i < dimc[1]; i++)
	    csplit[i] = &(csplit2[i * dimc[0]]);
    }
    xmiss = (const int **) ALLOC((int) dimx[1], sizeof(int *));
    xdata = (const double **) ALLOC((int) dimx[1], sizeof(double *));
    for (i = 0; i < dimx[1]; i++) {
	xmiss[i] = &(xmiss2[i * dimx[0]]);
	xdata[i] = &(xdata2[i * dimx[0]]);
    }
    Rprintf("inside honest C code fn: 6\n");
    
    for (i = 0; i < n; i++) {
	node = 1;               /* current node of the tree */
next:
	for (npos = 0; nnum[npos] != node; npos++);  /* position of the node */
  
  n1[npos]++;
  wt1[npos] += wt2[i];
  trs[npos] += wt2[i];
  cons[npos] += wt2[i];
  trsums[npos] += wt2[i]* y2[i];
  consums[npos] += wt2[i]* y2[i];
  trsqrsums[npos] +=  wt2[i]* y2[i] * y2[i];
  consqrsums[npos] += wt2[i] * y2[i] * y2[i];
  
  Rprintf("inside honest C code fn: 7\n");
  
       /* walk down the tree */
	nspl = nodes[3][npos] - 1;      /* index of primary split */
	if (nspl >= 0) {        /* not a leaf node */
	    var = vnum[nspl] - 1;
	    if (xmiss[var][i] == 0) {   /* primary var not missing */
		ncat = (int) split[1][nspl];
		temp = split[3][nspl];
		if (ncat >= 2)
		    dir = csplit[(int) xdata[var][i] - 1][(int) temp - 1];
		else if (xdata[var][i] < temp)
		    dir = ncat;
		else
		    dir = -ncat;
		if (dir) {
		    if (dir == -1)
			node = 2 * node;
		    else
			node = 2 * node + 1;
		    goto next;
		}
	    }
	    if (*usesur > 0) {
	      Rprintf("inside honest C code fn: 8\n");
	      
		for (j = 0; j < nodes[2][npos]; j++) {
		    nspl = nodes[1][npos] + nodes[3][npos] + j;
		    var = vnum[nspl] - 1;
		    if (xmiss[var][i] == 0) {   /* surrogate not missing */
			ncat = (int) split[1][nspl];
			temp = split[3][nspl];
			if (ncat >= 2)
			    dir = csplit[(int)xdata[var][i] - 1][(int)temp - 1];
			else if (xdata[var][i] < temp)
			    dir = ncat;
			else
			    dir = -ncat;
			if (dir) {
			    if (dir == -1)
				node = 2 * node;
			    else
				node = 2 * node + 1;
			    goto next;
			}
		    }
		}
	    }
	    if (*usesur > 1) {  /* go with the majority */
       Rprintf("inside honest C code fn: 9\n");
	      
		for (j = 0; nnum[j] != (2 * node); j++);
		//lcount = nodes[0][j];
		lcount=n1[j];
		for (j = 0; nnum[j] != (1 + 2 * node); j++);
		//rcount = nodes[0][j];
		rcount=n1[j];
		if (lcount != rcount) {
		    if (lcount > rcount)
			node = 2 * node;
		    else
			node = 2 * node + 1;
		    goto next;
		}
	    }
	}
	where[i] = node; //npos + 1;
    }
    Rprintf("inside honest C code fn: 10\n");
    
    //change yval1,dev1 to rpart terms
    for (i = 0; i <= nnodemax; i++) {
      if (invertdx[i] == -1)
        continue;
      int origindx = invertdx[i];
      //base case
      if (trs[origindx] != 0 && cons[origindx] != 0) {
        double tr_mean = trsums[origindx] * 1.0 / trs[origindx];
        double con_mean = consums[origindx] * 1.0 / cons[origindx];
        yval1[origindx] = tr_mean; // - con_mean;
        dev1[origindx] = trsqrsums[origindx] - trs[origindx] * tr_mean * tr_mean; 
          //+ consqrsums[origindx] - cons[origindx] * con_mean * con_mean;
      } else {
        int parentdx = invertdx[i / 2];
        yval1[origindx] = yval1[parentdx];
        dev1[origindx] = yval1[parentdx];
      }
    }
    Rprintf("inside honest C code fn: 11\n");
    
    
}

#include <Rinternals.h>

SEXP
honest_estimate_rparttree(SEXP dimx, SEXP nnode, SEXP nsplit, SEXP dimc,
	   SEXP nnum, SEXP nodes2, SEXP n1, SEXP wt1, SEXP dev1, SEXP yval1,
	   SEXP vnum, SEXP split2,
	   SEXP csplit2, SEXP usesur, SEXP xdata2, SEXP wt2, SEXP y2,
	   SEXP xmiss2)
{
    Rprintf("inside honest C code wrapper fn\n");
    int n = asInteger(dimx);
    SEXP where = PROTECT(allocVector(INTSXP, n));
    honest_estimate_rparttree0(INTEGER(dimx), asInteger(nnode), asInteger(nsplit),
		INTEGER(dimc), INTEGER(nnum), INTEGER(nodes2),
		INTEGER(vnum), REAL(split2), INTEGER(csplit2),
		INTEGER(usesur),
		INTEGER(n1), REAL(wt1), REAL(dev1), REAL(yval1),
		REAL(xdata2), REAL(wt2),REAL(y2),
		INTEGER(xmiss2),
		INTEGER(where));
    UNPROTECT(1);
    return where;
}
