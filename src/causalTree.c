/*
 * The main entry point for recursive partitioning routines.
 *
 * Input variables:
 *      ncat    = # categories for each var, 0 for continuous variables.
 *      split_Rule = 1 - TOT
 *                   2 - CT
 *                   3 - fit
 *                   4 - tstats
 *      Numbuckets = 0 - no discrete
 *                   o.w. - discrete
 *      
 *      crossmeth = 1 - TOT
 *                  2 - matching
 *                  3 - fitH
 *                  4 - fitA
 *                  5 - CTH
 *                  6 - CTA
 *                  
 *      opt  =  vector of options.  Same order as causalTree.control, as a vector
 *                   of doubles.
 *      minsize = minimum number of treated observations, control observations in a leaf
 *      p = propensity score
 *      xvals  = number of cross-validations to do
 *      xgrp  = indices for the cross-validations
 *      ymat  = vector of response variables
 *      xmat  = matrix of continuous variables
 *      wt      = vector of case weights
 *      treatment = vector of case treatment status: 1= treated, 0 = control
 *      ny   = number of columns of the y matrix (it is passed in as a vector)
 *      cost = nost
 *      xvar = vector of x features variance in column
 *      alpha = weight parameter for error function 
 *
 * Returned: a list with elements
 *      which = vector of final node numbers for each input obs
 *      cptable = the complexity table
 *      dsplit = for each split, numeric variables (doubles)
 *      isplit = for each split, integer variables
 *      dnode =  for each node, numeric variables
 *      inode =  for each node, integer variables
 *
 *   Naming convention: ncat = pointer to an integer vector, ncat2 = the
 *   input R object (SEXP) containing that vector, ncat3 = an output S object
 *   containing that vector.
 */

#define MAINRP
#include <math.h>
#include<stdio.h>
#include "causalTree.h"
#include "node.h"
#include "func_table.h"
#include "causalTreeproto.h"

SEXP
causalTree(SEXP ncat2, SEXP split_Rule2, SEXP bucketnum2, SEXP bucketMax2, SEXP method2, 
           SEXP crossmeth2, SEXP crosshonest2, SEXP opt2,
           SEXP minsize2, SEXP p2, SEXP xvals2, SEXP xgrp2,
        SEXP ymat2, SEXP xmat2, SEXP wt2, SEXP treatment2, SEXP ny2, SEXP cost2, 
        SEXP xvar2, SEXP split_alpha2, SEXP cv_alpha2, SEXP NumHonest2, SEXP gamma2)
{
    pNode tree;          /* top node of the tree */
    char *errmsg;
    int i, j, k, n;
    int maxcat;
    double temp, temp2;
    int *savesort = NULL /* -Wall */ ;
    double *dptr;               /* temp */
    int *iptr;
    int method;
    int split_Rule;
    int bucketnum;  
    int bucketMax;
    int crossmeth;
    int crosshonest;
    /*
     * pointers to R objects
     */
    int *ncat, *xgrp;
    int xvals;

    double *wt;
    double *treatment;
    int minsize;
    /* add propensity score: */
    double propensity;
    double split_alpha, cv_alpha;
    double gamma;
    int NumHonest;
    double train_to_est_ratio = 0.;
    
    /*
     * Return objects for R -- end in "3" to avoid overlap with internal names
     */
    SEXP which3, cptable3, dsplit3, isplit3, csplit3 = R_NilValue, /* -Wall */
         dnode3, inode3;

    /* work arrays for the return process */
    int nodecount, catcount, splitcount;
    double **ddnode, *ddsplit[3];
    int *iinode[6], *iisplit[3];
    int **ccsplit;
    double scale;
    
    CpTable cp;

    ncat = INTEGER(ncat2);
    xgrp = INTEGER(xgrp2);
    xvals = asInteger(xvals2);
    wt = REAL(wt2);
    treatment = REAL(treatment2);
    minsize = asInteger(minsize2);
    propensity = asReal(p2);
    split_alpha = asReal(split_alpha2);
    cv_alpha = asReal(cv_alpha2);
    gamma=asReal(gamma2);
    method = asInteger(method2); 
    crossmeth = asInteger(crossmeth2);
    crosshonest = asInteger(crosshonest2);
    split_Rule = asInteger(split_Rule2);
    bucketnum  = asInteger(bucketnum2);
    bucketMax = asInteger(bucketMax2);
    NumHonest = asInteger(NumHonest2);

    int split_id, cv_id;
    char getchar1;
    
    R_FlushConsole();
    //R_Process();
    //Rprintf("test print\n");
    //printf("test print2\n");
    //getchar1=getchar();
    /*
     * initialize the splitting functions from the function table
     */
    
    if (split_Rule <= NUM_SPLIT_RULE && crossmeth <= NUM_CROSSMETH) {
        split_id = split_Rule - 1;
        cv_id = crossmeth - 1;
        ct_init = split_func_table[split_id].init_split;
        ct_choose = split_func_table[split_id].choose_split;
        ct_eval = split_func_table[split_id].eval;
        ct_xeval = cv_func_table[cv_id].xeval;
        ct.num_y = asInteger(ny2);
    } else {
        error(_("Invalid value for 'split.Rule' or 'cv.option' "));
    }
    
    
    
    /*
     * set some other parameters
     */
    dptr = REAL(opt2);
    ct.min_node = (int) dptr[1];
    ct.min_split = (int) dptr[0];
    ct.complexity = dptr[2];
    ct.maxpri = (int) dptr[3] + 1;      /* max primary splits =
                                           max competitors + 1 */
    if (ct.maxpri < 1) ct.maxpri = 1;
    ct.maxsur = (int) dptr[4];
    ct.usesurrogate = (int) dptr[5];
    ct.sur_agree = (int) dptr[6];
    ct.maxnode = (int) pow((double) 2.0, (double) dptr[7]) - 1;
    ct.n = nrows(xmat2);
    ct.NumHonest = NumHonest;
    
    n = ct.n;                   
    ct.nvar = ncols(xmat2);
    ct.numcat = INTEGER(ncat2);
    ct.wt = wt;
    ct.treatment = treatment;
    ct.iscale = 0.0;
    ct.vcost = REAL(cost2);
    
    ct.xvar = REAL(xvar2);
    ct.NumXval = xvals;
    
    
    dptr = REAL(xmat2);
    ct.xdata = (double **) ALLOC(ct.nvar, sizeof(double *));
    for (i = 0; i < ct.nvar; i++) {
        ct.xdata[i] = dptr;
        dptr += n;
    }
    
    ct.ydata = (double **) ALLOC(n, sizeof(double *));
    
    dptr = REAL(ymat2);
    temp2 = 0;
    for (i = 0; i < n; i++) {
        ct.ydata[i] = dptr;
        for (j = 0; j < ct.num_y; j++) {
            if (fabs(ct.ydata[i][j]) > temp2) temp2 = fabs(ct.ydata[i][j]);        
        }
        dptr += ct.num_y;
    }
    ct.max_y = temp2;
    ct.propensity = propensity;
    
    /*
     * allocate some scratch
     */
    
    ct.tempvec = (int *) ALLOC(n, sizeof(int));
    ct.xtemp = (double *) ALLOC(n, sizeof(double));
    ct.ytemp = (double **) ALLOC(n, sizeof(double *));
    ct.wtemp = (double *) ALLOC(n, sizeof(double));
    ct.trtemp = (double *) ALLOC(n, sizeof(double));
    
    /*
     * create a matrix of sort indices, one for each continuous variable
     *   This sort is "once and for all".
     * I don't have to sort the categoricals.
     */
    ct.sorts = (int **) ALLOC(ct.nvar, sizeof(int *));
    ct.sorts[0] = (int *) ALLOC(n * ct.nvar, sizeof(int));
    maxcat = 0;
    for (i = 0; i < ct.nvar; i++) {
        ct.sorts[i] = ct.sorts[0] + i * n;
        for (k = 0; k < n; k++) {
            if (!R_FINITE(ct.xdata[i][k])) {
                ct.tempvec[k] = -(k + 1);       /* this variable is missing */
                ct.xtemp[k] = 0;        /* avoid weird numerics in S's NA */
            } else {
                ct.tempvec[k] = k;
                ct.xtemp[k] = ct.xdata[i][k];
            }
        }
        if (ncat[i] == 0)
            mysort(0, n - 1, ct.xtemp, ct.tempvec);
        else if (ncat[i] > maxcat)
            maxcat = ncat[i];
        for (k = 0; k < n; k++)
            ct.sorts[i][k] = ct.tempvec[k];
    }

    /*
     * save away a copy of the ct.sorts, if needed for xval
     */

    if (xvals > 1) {
        savesort = (int *) ALLOC(n * ct.nvar, sizeof(int));
        memcpy(savesort, ct.sorts[0], n * ct.nvar * sizeof(int));
    }

    /*
     * And now the last of my scratch space
     */
    
    if (maxcat > 0) {
        ct.csplit = (int *) ALLOC(3 * maxcat, sizeof(int));
        ct.left = ct.csplit + maxcat;
        ct.right = ct.left + maxcat;
        ct.lwt = (double *) ALLOC(2 * maxcat, sizeof(double));
        ct.rwt = ct.lwt + maxcat;
        ct.ltr = (double *) ALLOC(2 * maxcat, sizeof(double));
        ct.rtr = ct.ltr + maxcat;
    } else
        ct.csplit = (int *) ALLOC(1, sizeof(int));

    /*
     * initialize the top node of the tree
     */
    errmsg = _("unknown error");
    which3 = PROTECT(allocVector(INTSXP, n));
    ct.which = INTEGER(which3);
    temp = 0;
    temp2 = 0; 

    for (i = 0; i < n; i++) {
        ct.which[i] = 1;
        temp += wt[i];
        temp2 += treatment[i];
    }
    
    train_to_est_ratio = 100;
    i = (*ct_init) (n, ct.ydata, maxcat, &errmsg, &ct.num_resp, 1, wt, treatment,
         bucketnum, bucketMax, &train_to_est_ratio);
    

    
    if (i > 0)
        error(errmsg);
   

    nodesize = sizeof(Node) + (ct.num_resp - 20) * sizeof(double);
    tree = (pNode) ALLOC(1, nodesize);
    memset(tree, 0, nodesize);
    tree->num_obs = n;
    tree->sum_wt = temp;
    tree->sum_tr = temp2;
    tree->parent = NULL;


    if (split_Rule == 1) {
        //tot:
        (*ct_eval) (n, ct.ydata, tree->response_est, tree->controlMean, tree->treatMean, 
         &(tree->risk), wt, treatment, ct.max_y, ct.propensity);
    } else if (split_Rule == 2) {
        // ct:
        (*ct_eval) (n, ct.ydata, tree->response_est, tree->controlMean, tree->treatMean, 
         &(tree->risk), wt, treatment, ct.max_y, split_alpha, train_to_est_ratio);
    } else if (split_Rule == 3) {
        //fit
        (*ct_eval) (n, ct.ydata, tree->response_est, tree->controlMean, tree->treatMean, 
         &(tree->risk), wt, treatment, ct.max_y, split_alpha, train_to_est_ratio);
    } else if (split_Rule == 4) {
        // tstats
        (*ct_eval) (n, ct.ydata, tree->response_est, tree->controlMean, tree->treatMean,
         &(tree->risk), wt, treatment, ct.max_y, split_alpha, train_to_est_ratio);
    } else if (split_Rule == 5) {
        // totD
        (*ct_eval) (n, ct.ydata, tree->response_est, tree->controlMean, tree->treatMean,
         &(tree->risk), wt, treatment, ct.max_y, ct.propensity);
    } else if (split_Rule == 6) {
        //CTD
        (*ct_eval) (n, ct.ydata, tree->response_est, tree->controlMean, tree->treatMean, 
         &(tree->risk), wt, treatment, ct.max_y, split_alpha, train_to_est_ratio);
    } else if (split_Rule == 7) {
        // fitD
        (*ct_eval) (n, ct.ydata, tree->response_est, tree->controlMean, tree->treatMean, 
         &(tree->risk), wt, treatment, ct.max_y, split_alpha, train_to_est_ratio);
    } else if (split_Rule == 8) {
        //tstatsD:
        (*ct_eval) (n, ct.ydata, tree->response_est, tree->controlMean, tree->treatMean, 
         &(tree->risk), wt, treatment, ct.max_y, split_alpha, train_to_est_ratio);
    } else if (split_Rule == 9) {
        // user (temporarily set as CT)
        (*ct_eval) (n, ct.ydata, tree->response_est, tree->controlMean, tree->treatMean, 
         &(tree->risk), wt, treatment, ct.max_y, split_alpha, train_to_est_ratio);
    } else if (split_Rule == 10) {
        // userD (temporarily set as CTD)
        (*ct_eval) (n, ct.ydata, tree->response_est, tree->controlMean, tree->treatMean, 
         &(tree->risk), wt, treatment, ct.max_y, split_alpha, train_to_est_ratio);
    }else if (split_Rule == 11) {
      // policy
      (*ct_eval) (n, ct.ydata, tree->response_est, tree->controlMean, tree->treatMean, 
       &(tree->risk), wt, treatment, ct.max_y, split_alpha, train_to_est_ratio);
    } else if (split_Rule == 12) {
      // policyD
      (*ct_eval) (n, ct.ydata, tree->response_est, tree->controlMean, tree->treatMean, 
       &(tree->risk), wt, treatment, ct.max_y, split_alpha, train_to_est_ratio);
    }
    tree->complexity = tree->risk;
    ct.alpha = ct.complexity * tree->risk;

    /*
     * Do the basic tree
     */
    
    partition(1, tree, &temp, 0, n, minsize, split_Rule, split_alpha, bucketnum, bucketMax,
              train_to_est_ratio); // temp store sumrisk
  
    CpTable cptable = (CpTable) ALLOC(1, sizeof(cpTable));

    cptable->cp = tree->complexity;
    cptable->risk = tree->risk;
    cptable->nsplit = 0;
    cptable->forward = 0;
    cptable->xrisk = 0;
    cptable->xstd = 0;
    ct.num_unique_cp = 1;
    
    if (tree->rightson) {
        make_cp_list(tree, tree->complexity, cptable);
        make_cp_table(tree, tree->complexity, 0);
        
        if (xvals > 1) {
            myxval(xvals, cptable, xgrp, maxcat, &errmsg, minsize, savesort, split_Rule,
                   crossmeth, split_alpha, cv_alpha, bucketnum, bucketMax,gamma);
        }
    }
    /*
     * all done, create the return list for R
     * first the cp table
     */
    scale = 1 / tree->risk;
    i = 0;
    cptable3 = PROTECT(allocMatrix(REALSXP, xvals > 1 ? 5 : 3,
                ct.num_unique_cp));
    dptr = REAL(cptable3);
    for (cp = cptable; cp; cp = cp->forward) {
        dptr[i++] = cp->cp * scale;
        dptr[i++] = cp->nsplit;
        dptr[i++] = cp->risk * scale;
        if (xvals > 1) {
            dptr[i++] = cp->xrisk * scale ;
            dptr[i++] = cp->xstd * scale ;
        }
    }

    /*
     * Return the body of the tree
     *  For each component we first create a vector to hold the
     *  result, then a ragged array index into the vector.
     * The ctmatrix routine then fills everything in.
     */
    ctcountup(tree, &nodecount, &splitcount, &catcount);
    dnode3 = PROTECT(allocMatrix(REALSXP, nodecount, (3 + ct.num_resp)));
    ddnode = (double **) ALLOC(3 + ct.num_resp, sizeof(double *));
    dptr = REAL(dnode3);
    for (i = 0; i < 3 + ct.num_resp; i++) {
        ddnode[i] = dptr;
        dptr += nodecount;
    }

    dsplit3 = PROTECT(allocMatrix(REALSXP, splitcount, 3));
    dptr = REAL(dsplit3);
    for (i = 0; i < 3; i++) {
        ddsplit[i] = dptr;
        dptr += splitcount;
        for (j = 0; j < splitcount; j++)
            ddsplit[i][j] = 0.0;
    }

    inode3 = PROTECT(allocMatrix(INTSXP, nodecount, 6));
    iptr = INTEGER(inode3);
    for (i = 0; i < 6; i++) {
        iinode[i] = iptr;
        iptr += nodecount;
    }

    isplit3 = PROTECT(allocMatrix(INTSXP, splitcount, 3));
    iptr = INTEGER(isplit3);
    for (i = 0; i < 3; i++) {
        iisplit[i] = iptr;
        iptr += splitcount;
    }

    if (catcount > 0) {
        csplit3 = PROTECT(allocMatrix(INTSXP, catcount, maxcat));
        ccsplit = (int **) ALLOC(maxcat, sizeof(int *));
        iptr = INTEGER(csplit3);
        for (i = 0; i < maxcat; i++) {
            ccsplit[i] = iptr;
            iptr += catcount;
            for (j = 0; j < catcount; j++)
                ccsplit[i][j] = 0;      /* zero it out */
        }
    } else {
        ccsplit = NULL;
    }

    ctmatrix(tree, ct.numcat, ddsplit, iisplit, ccsplit, ddnode, iinode, 1);
    free_tree(tree, 0);         /* let the memory go */

    /*
     * Fix up the 'which' array
     *  Nodes are sometimes trimmed during the
     *  tree building, and 'which' is not updated in that case
     */
    for (i = 0; i < n; i++) {
        k = ct.which[i];
        do {
            for (j = 0; j < nodecount; j++)
                if (iinode[0][j] == k) {
                    ct.which[i] = j + 1;
                    break;
                }
            k /= 2;
        } while (j >= nodecount);
    }

    /* Create the output list */
    int nout = catcount > 0 ? 7 : 6;
    SEXP rlist = PROTECT(allocVector(VECSXP, nout));
    SEXP rname = allocVector(STRSXP, nout);
    setAttrib(rlist, R_NamesSymbol, rname);
    SET_VECTOR_ELT(rlist, 0, which3);
    SET_STRING_ELT(rname, 0, mkChar("which"));
    SET_VECTOR_ELT(rlist, 1, cptable3);
    SET_STRING_ELT(rname, 1, mkChar("cptable"));
    SET_VECTOR_ELT(rlist, 2, dsplit3);
    SET_STRING_ELT(rname, 2, mkChar("dsplit"));
    SET_VECTOR_ELT(rlist, 3, isplit3);
    SET_STRING_ELT(rname, 3, mkChar("isplit"));
    SET_VECTOR_ELT(rlist, 4, dnode3);
    SET_STRING_ELT(rname, 4, mkChar("dnode"));
    SET_VECTOR_ELT(rlist, 5, inode3);
    SET_STRING_ELT(rname, 5, mkChar("inode"));
    if (catcount > 0) {
        SET_VECTOR_ELT(rlist, 6, csplit3);
        SET_STRING_ELT(rname, 6, mkChar("csplit"));
    }

    UNPROTECT(1 + nout);
    return rlist;
}
