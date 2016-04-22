/*
 * Walk an observation 'one more split' down the tree.  If there are no
 *   more splits, return 0, otherwise return the address of the new node.
 * A return of zero also comes about if surrogates aren't being used, and I
 *   hit a missing value.
 *
 * tree :  the current node
 * obs  :  the observation number.  This will be negative iff the primary
 *           split is missing.
 *
 * This is the one routine that accesses the tree by observation number,
 *  without the mediation of the ct.sorts array.  For missing value
 *  information it thus has to look at X directly using a macro.
 */
#include "causalTree.h"
#include "node.h"
#include "causalTreeproto.h"

pNode
branch(pNode tree, int obs)
{
    int j, dir;
    int category;               /* for categorical variables */
    pNode me;
    pSplit tsplit;
    double **xdata;

    if (!tree->leftson) return NULL;

    me = tree;
    xdata = ct.xdata;

   /*
    * choose left or right son
    *   this may use lots of surrogates before we're done
    */
    tsplit = me->primary;
    j = tsplit->var_num;
    if (R_FINITE(xdata[j][obs])) {
	if (ct.numcat[j] == 0) {        /* continuous */
	    dir = (xdata[j][obs] < tsplit->spoint) ?
		tsplit->csplit[0] : -tsplit->csplit[0];
	    goto down;
	} else {                /* categorical predictor */
	    category = (int) xdata[j][obs];     /* factor predictor -- which
						 * level? */
	    dir = (tsplit->csplit)[category - 1];
	    if (dir)
		goto down;
	}
    }
    if (ct.usesurrogate == 0)
	return NULL;
   /*
    * use the surrogates
    */
    for (tsplit = me->surrogate; tsplit; tsplit = tsplit->nextsplit) {
	j = tsplit->var_num;
	if (R_FINITE(xdata[j][obs])) {  /* not missing */
	    if (ct.numcat[j] == 0) {
		dir = (ct.xdata[j][obs] < tsplit->spoint) ?
		    tsplit->csplit[0] : -tsplit->csplit[0];
		goto down;
	    } else {
		category = (int) xdata[j][obs]; /* factor predictor -- which
						 * level */
		dir = (tsplit->csplit)[category - 1];
		if (dir)
		    goto down;
	    }
	}
    }


    if (ct.usesurrogate < 2)
	return NULL;
    /*
     * split it by default
     */
    dir = me->lastsurrogate;

down:
    return (dir == LEFT) ? me->leftson : me->rightson;
}
