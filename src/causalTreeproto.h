/*
 * prototypes for all of the causalTree functions
 *   This helps the ansi compiler do tight checking.
 *
 */
#include "node.h"

pNode branch(pNode tree, int obs);

//void bsplit(pNode me, int n1, int n2);
void bsplit(pNode me, int n1, int n2, int minsize, int split_Rule, double alpha,
            int bucketnum, int bucketMax, double train_to_est_ratio);

void choose_surg(int n1, int n2, int *y, double *x, int *order,
		 int ncat, double *agreement, double *split, int *csplit,
		 double ltot, double rtot, double *adj);

void fix_cp(pNode me, double parent_cp);

void free_tree(pNode node, int freenode);

void graycode_init0(int maxcat);
void graycode_init1(int numcat, int *count);
void graycode_init2(int numcat, int *count, double *val);
int graycode(void);

pSplit insert_split(pSplit *listhead, int ncat, double improve, int max);

void make_cp_list(pNode me, double parent, CpTable cptable_head);

CpTable make_cp_table(pNode me, double parent, int nsplit);

void mysort(int start, int stop, double *x, int *cvec);

void nodesplit(pNode me, int nodenum, int n1, int n2, int *nleft, int *nright);

//int partition(int nodenum, pNode splitnode, double *sumrisk, int n1, int n2);
int partition(int nodenum, pNode splitnode, double *sumrisk, int n1, int n2, 
              int minsize, int split_Rule, double alpha, int bucketnum, int bucketMax, 
              double train_to_est_ratio);

int print_tree(pNode me, int maxdepth);


SEXP causalTree(SEXP ncat2, SEXP split_Rule2, SEXP bucketnum2, SEXP bucketMax2, SEXP method2, 
                SEXP crossmeth2, SEXP crosshonest2, SEXP opt2, SEXP minsize2, SEXP p2, 
                SEXP ymat2, SEXP xmat2, SEXP xvals2, SEXP xgrp2, SEXP wt2, SEXP treatment2, SEXP ny2,
                SEXP cost2, SEXP xvar2, SEXP split_alpha2, SEXP cv_alpha2, SEXP NumHonest2);

void causalTree_callback0(int *nr);
void causalTree_callback1(int n, double *y[], double *wt, double *z);
void causalTree_callback2(int n, int ncat, double *y[], double *wt,
		     double *x, double *good);
void ctcountup(pNode me, int *nnode, int *nsplit, int *ncat);

void ctmatrix(pNode me, int *numcat, double **dsplit, int **isplit,
	      int **csplit, double **dnode, int **inode, int id);


void totrundown(pNode tree, int obs, double *cp, double *xpred, double *xtemp);
void matching_rundown(pNode tree, int obs, int neighbor, double *cp, double *xpred, 
                     double *xpred2, double *xtemp);
void fitH_rundown(pNode tree, int obs, double *cp, double *xpred, double *xtemp, int k,
                  double alpha, double xtrain_to_est_ratio);
void fitA_rundown(pNode tree, int obs, double *cp, double *xpred, double *xtemp, int k);

void CTH_rundown(pNode tree, int obs, double *cp, double *xpred, double *xtemp, int k, 
                 double alpha, double xtrain_to_est_ratio, double propensity);
void CTA_rundown(pNode tree, int obs, double *cp, double *xpred, double *xtemp, int k, double alpha);
void surrogate(pNode me, int n1, int n2);


void myxval(int n_xval, CpTable cptable_head, int *x_grp, int maxcat, char **errmsg, 
           int minsize, int *savesort, int split_Rule,
           int crossmeth, double split_alpha, double cv_alpha, int bucketnum, int bucketMax);
/* ---------------------- for xvalHelper --------------------- */
int findNeighbor(int obs, int k);
double measureDistance(int i, int j);
