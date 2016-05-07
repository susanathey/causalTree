/*
 * This rundown function for CTH.
 *
 */
#include "causalTree.h"
#include "node.h"
#include "causalTreeproto.h"

#ifdef NAN
/* NAN is supported */
#endif

void
CTH_rundown(pNode tree, int obs, double *cp, double *xpred, double *xtemp, int k, double alpha, 
            double xtrain_to_est_ratio, double propensity)
{
    int i, obs2 = (obs < 0) ? -(1 + obs) : obs;
    int my_leaf_id;
    pNode otree =  tree;
    pNode otree_tmp = tree;
    pNode tree_tmp = tree;
    
    int opnumber = 0;
    int j, s;
    int tmp_obs, tmp_id;
    double tr_mean, con_mean;
    double tr_sqr_sum, con_sqr_sum;
    double consums, trsums, cons, trs;
    double tr_var, con_var;

    /*
     * Now, repeat the following: for the cp of interest, run down the tree
     *   until I find a node with smaller complexity.  The parent node will
     *   not have collapsed, but this split will have, so this is my
     *   predictor.
     */
    for (i = 0; i < ct.num_unique_cp; i++) {
        cons = 0.;
        trs = 0.;
        consums = 0.;
        trsums = 0.;
        tr_sqr_sum = 0.;
        con_sqr_sum = 0.;
        
        while (cp[i] < tree->complexity) {
	        tree = branch(tree, obs);
	        if (tree == 0)
		        goto oops;
	        otree = tree;
	    }
	    xpred[i] = tree->response_est[0];
        my_leaf_id = tree->id;
        
        for (s = k; s < ct.n; s++) {
            tree_tmp = otree_tmp;
            j = ct.sorts[0][s];
            tmp_obs = (j < 0) ? -(1 + j) : j;
            while (cp[i] < tree_tmp->complexity) {
                tree_tmp = branch(tree_tmp, tmp_obs);
            }
            tmp_id = tree_tmp->id;
            if (tmp_id == my_leaf_id) {
                if (ct.treatment[tmp_obs] == 0) {
                    cons += ct.wt[tmp_obs];
                    consums += *ct.ydata[tmp_obs] * ct.wt[tmp_obs];
                    con_sqr_sum += (*ct.ydata[tmp_obs]) * (*ct.ydata[tmp_obs]) * ct.wt[tmp_obs];
                } else {
                    trs += ct.wt[tmp_obs];
                    trsums += *ct.ydata[tmp_obs] * ct.wt[tmp_obs];
                    tr_sqr_sum += (*ct.ydata[tmp_obs]) * (*ct.ydata[tmp_obs]) * ct.wt[tmp_obs];
                }
            }
        }

        if (trs == 0) {
            tr_mean = tree->parent->xtreatMean[0];
            tr_var = 0;
        } else {
            tr_mean = trsums / trs;
            tree->xtreatMean[0] = tr_mean;
            tr_var = tr_sqr_sum / trs - tr_mean * tr_mean;
        }
        
        if (cons == 0) {
            con_mean = tree->parent->xcontrolMean[0];
            con_var = 0;
        } else {
            con_mean = consums / cons;
            tree->xcontrolMean[0] = con_mean;
            con_var = con_sqr_sum / cons - con_mean * con_mean;
        }
        
        xtemp[i] = (*ct_xeval)(ct.ydata[obs2], ct.wt[obs2], ct.treatment[obs2], tr_mean, 
                    con_mean, trs, cons, alpha, xtrain_to_est_ratio, propensity);
    }
    return;

oops:;
    if (ct.usesurrogate < 2) {  /* must have hit a missing value */
	for (i = 0; i < ct.num_unique_cp; i++)
	    xpred[i] = otree->response_est[0];

	xtemp[i] = (*ct_xeval)(ct.ydata[obs2], ct.wt[obs2], ct.treatment[obs2], tr_mean, con_mean);
	Rprintf("oops number %d.\n", opnumber++);
  return;
    }
    warning("Warning message--see rundown.c");
}
