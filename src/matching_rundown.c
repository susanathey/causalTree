/*
 * This rundown function is typically for TOT method. You may change it to make it 
 * compatibel with other splitting funcitons.
 *
 */
#include "causalTree.h"
#include "node.h"
#include "causalTreeproto.h"

void
matching_rundown(pNode tree, int obs, int neighbor, double *cp, double *xpred, 
                 double *xpred2, double *xtemp)
{
    int i, obs2 = (obs < 0) ? -(1 + obs) : obs;
    pNode otree =  tree;
    pNode tree2 = tree;
    pNode otree2 = tree;
    int neighbor2 = (neighbor < 0)? -( 1 + neighbor) : neighbor;
    // for debug only:
    int opnumber = 0;

    /*
     * Now, repeat the following: for the cp of interest, run down the tree
     *   until I find a node with smaller complexity.  The parent node will
     *   not have collapsed, but this split will have, so this is my
     *   predictor.
     */
    for (i = 0; i < ct.num_unique_cp; i++) {
        while (cp[i] < tree->complexity) {
	        tree = branch(tree, obs);
	        if (tree == 0)
		        goto oops;
	        otree = tree;
	    }
        xpred[i] = tree->response_est[0];
        
        while (cp[i] < tree2->complexity) {
            tree2 = branch(tree2, neighbor);
            if (tree2 == 0) 
                goto oops;
            otree2 = tree2;
        }
        
        xpred2[i] = tree2->response_est[0]; 
        // temp = (2 * ct.treatment[obs2] - 1) * (*ct.ydata[obs2] - *ct.ydata[neighbor2]);
        // xtemp[i] = (temp - 0.5 * (xpred[i] + xpred2[i])) * (temp - 0.5 * (xpred[i] + xpred2[i]));
        xtemp[i] = (*ct_xeval)(ct.ydata[obs2], ct.ydata[neighbor2], ct.wt[obs2], ct.wt[neighbor2],
                    ct.treatment[obs2], ct.treatment[neighbor2], xpred[i], xpred2[i]);
    }
    return;

oops:;
    if (ct.usesurrogate < 2) {  /* must have hit a missing value */
	for (i = 0; i < ct.num_unique_cp; i++)
	    xpred[i] = otree->response_est[0];
        xpred2[i] = otree2->response_est[0];

	xtemp[i] = (*ct_xeval)(ct.ydata[obs2], ct.ydata[neighbor2], ct.wt[obs2], ct.wt[neighbor2],
             ct.treatment[obs2], ct.treatment[neighbor2], xpred[i], xpred2[i]);
	Rprintf("oops number %d.\n", opnumber++);
  return;
    }
    warning("Warning message--see rundown.c");
}
