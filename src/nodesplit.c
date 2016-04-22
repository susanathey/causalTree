/*
 * nodesplit -- Split the node in two, and keep a count as we do of how
 *  many splits are determined by each surrogate variable.
 *
 * me	  : pointer to the node of the tree
 * nodenum: the node number of the current node, used to update "which"
 * n1, n2 : starting and ending indices for the observation numbers
 * nnleft, nnright: at the end, how many obs were sent to the left and
 *            to the right.  Beware-  for certain settings of the
 *            usesurrogate option, some obs go neither left nor right
 */
#include "causalTree.h"
#include "node.h"
#include "causalTreeproto.h"
#include <stdio.h>

void
nodesplit(pNode me, int nodenum, int n1, int n2, int *nnleft, int *nnright)
{
    int i, j, k;
    pSplit tsplit;
    int var, extra, lastisleft, someleft;
    int i1, i2, i3;
    int leftson, rightson;
    int pvar;
    double psplit;
    int *index;
    int *which;
    int **sorts;
    int *sindex;                /* sindex[i] is a shorthand for sorts[var][i] */
    double **xdata;
    int nleft, nright;
    //int           dummy;      /* debugging */

    which = ct.which;
    sorts = ct.sorts;
    xdata = ct.xdata;
    leftson = 2 * nodenum;      /* the label that will go with the left son */
    rightson = leftson + 1;

    /*
     * Walk through the variables (primary, then surrogate 1, then surr 2...)
     *   and reassign "which"
     */
    tsplit = me->primary;
    pvar = tsplit->var_num;     /* primary variable */
    someleft = 0;
    nleft = 0;
    nright = 0;

    if (ct.numcat[pvar] > 0) {  /* categorical primary variable */
	index = tsplit->csplit;
	for (i = n1; i < n2; i++) {
	    j = sorts[pvar][i];
	    if (j < 0)
		someleft++;     /* missing value */
	    else
		switch (index[(int) xdata[pvar][j] - 1]) {
		case LEFT:
		    which[j] = leftson;
		    nleft++;
		    break;
		case RIGHT:
		    which[j] = leftson + 1;
		    nright++;
		    break;
		}
	}
    } else {
	psplit = tsplit->spoint;        /* value of split point */
	extra = tsplit->csplit[0];
	for (i = n1; i < n2; i++) {
	    j = sorts[pvar][i];
	    if (j < 0)
		someleft++;
	    else {
		if (xdata[pvar][j] < psplit)
		    k = extra;
		else
		    k = -extra;
		if (k == LEFT) {
		    which[j] = leftson;
		    nleft++;
		} else {
		    which[j] = leftson + 1;
		    nright++;
		}
	    }
	}
    }

    /*
     * Now the surrogates
     *   Usually, there aren't a lot of observations that need to
     *   be split.  So it is more efficient to make one 1:n pass,
     *   with multiple runs through the surrogate list.
     */
    if (someleft > 0 && ct.usesurrogate > 0) {
	for (i = n1; i < n2; i++) {
	    j = ct.sorts[pvar][i];
	    if (j >= 0)
		continue;       /* already split */

	    j = -(j + 1);       /* obs number - search for surrogates */
	    for (tsplit = me->surrogate; tsplit; tsplit = tsplit->nextsplit) {
		var = tsplit->var_num;
		if (!R_FINITE(xdata[var][j]))
		    continue;
		/* surrogate not missing - process it */

		if (ct.numcat[var] > 0) {       /* categorical surrogate */
		    index = tsplit->csplit;
		    k = (int) xdata[var][j];    /* the value of the surrogate  */
		    /*
		     * The need for the if stmt below may not be obvious.
		     * The surrogate's value must not be missing, AND there
		     * must have been at least 1 person with both this
		     * level of the surrogate and a primary split value
		     * somewhere in the node.  If everyone in this node
		     * with level k of the surrogate also had a missing
		     * value of the primary variable, then index[k-1] will
		     * be zero.
		     */
		    if (index[k - 1]) {
			tsplit->count++;
			if (index[k - 1] == LEFT) {
			    which[j] = leftson;
			    nleft++;
			} else {
			    which[j] = leftson + 1;
			    nright++;
			}
			someleft--;
			break;
		    }
		} else {
		    psplit = tsplit->spoint;    /* continuous surrogate */
		    extra = tsplit->csplit[0];
		    tsplit->count++;
		    if (xdata[var][j] < psplit)
			k = extra;
		    else
			k = -extra;
		    if (k == LEFT) {
			which[j] = leftson;
			nleft++;
		    } else {
			which[j] = rightson;
			nright++;
		    }
		    someleft--;
		    break;
		}
	    }
	}
    }
    if (someleft > 0 && ct.usesurrogate == 2) {
	/* all surrogates missing, use the default */
	i = me->lastsurrogate;
	if (i) {           /* 50-50 splits are possible - there is no
			    * "default" */
	    if (i < 0) {
		lastisleft = leftson;
		nleft += someleft;
	    } else {
		lastisleft = rightson;
		nright += someleft;
	    }

	    for (i = n1; i < n2; i++) {
		j = sorts[pvar][i];
		/*
		 * only those who weren't split by the primary (j < 0) and
		 * weren't split by a surrogate (which == nodenum) need to be
		 * assigned
		 */
		if (j < 0) {
		    j = -(j + 1);
		    if (which[j] == nodenum)
			which[j] = lastisleft;
		}
	    }
	}
    }
    /*
     * Last part of the work is to update the sorts matrix
     *
     * Say that n1=5, n2=12, 4 go left, 3 go right, and one obs
     *   stays home, and the data looks like this:
     *
     *   sorts[var][5] = 17    which[17]= 17 = 2*nodenum +1 = rightson
     *       "            4    which[4] = 16 = 2*nodenum    = leftson
     *                   21          21   16
     *                    6           6   17
     *      "		  7           7   16
     *                   30          30   17
     *                    8           8   16
     *	sorts[var][12]=	-11    which[11] = 8 = nodenum  (X = missing)
     *
     *  Now, every one of the rows of the sorts contains these same
     *    elements -- 4,6,7,8,11,17,21,30 -- in some order, which
     *    represents both how they are sorted and any missings via
     *    negative numbers.
     *  We need to reorder this as {{goes left}, {goes right}, {stays
     *    here}}, preserving order within the first two groups.  The
     *    order within the "stay here" group doesn't matter since they
     *    won't be looked at again for splitting.
     *  So the result in this case should be
     *       4, 21, 7, 8,   17, 6, 30,  -11
     *  The algorithm is the opposite of a merge sort.
     *
     *  Footnote: if no surrogate variables were used, then one could
     *   skip this process for the primary split variable, as that
     *   portion of "sorts" would remain unchanged.  It's not worth
     *   the bother of checking, however.
     */
    for (k = 0; k < ct.nvar; k++) {
	sindex = ct.sorts[k];   /* point to variable k */
	i1 = n1;
	i2 = i1 + nleft;
	i3 = i2 + nright;
	for (i = n1; i < n2; i++) {
	    j = sindex[i];
	    if (j < 0)
		j = -(j + 1);
	    if (which[j] == leftson)
		sindex[i1++] = sindex[i];
	    else {
		if (which[j] == rightson)
		    ct.tempvec[i2++] = sindex[i];
		else
		    ct.tempvec[i3++] = sindex[i];       /* went nowhere */
	    }
	}
	for (i = n1 + nleft; i < n2; i++)
	    sindex[i] = ct.tempvec[i];
    }

    *nnleft = nleft;
    *nnright = nright;
}
