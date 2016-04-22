/*
 *  Do causalTree predictions given the matrix form of the tree.
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
#include "causalTree.h"
#include "causalTreeproto.h"

    static void
estimate_causalTree0(const int *dimx, int nnode, 
                     int nsplit, const int *dimc, const int *nnum, const int *nodes2,
                     const int *vnum,
        const double *split2, const int *csplit2, const int *usesur,
        const double *xdata2, const int *xmiss2, int *where)
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

    n = dimx[0];
    for (i = 0; i < 4; i++) {
        nodes[i] = &(nodes2[nnode * i]);
        split[i] = &(split2[nsplit * i]);
    }

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

    for (i = 0; i < n; i++) {
        node = 1;               /* current node of the tree */
next:
        for (npos = 0; nnum[npos] != node; npos++);  /* position of the node */
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
                for (j = 0; nnum[j] != (2 * node); j++);
                lcount = nodes[0][j];
                for (j = 0; nnum[j] != (1 + 2 * node); j++);
                rcount = nodes[0][j];
                if (lcount != rcount) {
                    if (lcount > rcount)
                        node = 2 * node;
                    else
                        node = 2 * node + 1;
                    goto next;
                }
            }
        }
        where[i] = node;
    }
}

#include <Rinternals.h>

    SEXP
estimate_causalTree(SEXP dimx, SEXP nnode, SEXP nsplit, SEXP dimc,
        SEXP nnum, SEXP nodes2, SEXP vnum, SEXP split2,
        SEXP csplit2, SEXP usesur, SEXP xdata2, SEXP xmiss2)
{
    int n = asInteger(dimx);
    SEXP where = PROTECT(allocVector(INTSXP, n));
    estimate_causalTree0(INTEGER(dimx), asInteger(nnode), asInteger(nsplit),
            INTEGER(dimc), INTEGER(nnum), INTEGER(nodes2),
            INTEGER(vnum), REAL(split2), INTEGER(csplit2),
            INTEGER(usesur), REAL(xdata2), INTEGER(xmiss2),
            INTEGER(where));
    UNPROTECT(1);
    return where;
}
