/*
 *  For S's usage, convert the linked list data into matrix form
 */
#include "causalTree.h"
#include "node.h"
#include "causalTreeproto.h"

/* These four preserve from call to call */
static int ncnt, scnt, ccnt;
static double cp_scale;

void
ctmatrix(pNode me, int *numcat, double **dsplit,
	 int **isplit, int **csplit, double **dnode, int **inode, int id)
{
    /*
     * dsplit  0: improvement
     *         1: split point if continuous; index into csplit if not
     *         2: surrogate: adjusted agreement,  primary: nothing
     * isplit  0: variable #
     *         1: count
     *         2: if continuous: direction -1=left, 1=right
     *            if categorical: # of categories
     * csplit[i]: -1=left, 0=missing category, 1=right
     * dnode   0: risk
     *         1: complexity threshold
     *         2: sum of weights
     *         3, ...: response estimate
     * inode   0: node number
     *         1: index of the first primary, in the split list
     *         2: #primary    ==0 if this is a terminal node
     *         3: #surrogates
     *         4: # observations
     *         5: # obs for which this is the final resting place
     */

    int i, j, k;
    pSplit spl;

    if (id == 1) {              /* this is the top node */
	cp_scale = 1 / me->risk;
	scnt = 0;
	ncnt = 0;
	ccnt = 0;
    }
    dnode[0][ncnt] = me->risk;
    dnode[1][ncnt] = me->complexity * cp_scale;
    dnode[2][ncnt] = me->sum_wt;
    for (i = 0; i < ct.num_resp; i++)
	dnode[3 + i][ncnt] = me->response_est[i];
    inode[0][ncnt] = id;
    inode[4][ncnt] = me->num_obs;

    if (me->complexity <= ct.alpha || me->leftson == 0) {       /* no kids */
	inode[1][ncnt] = 0;
	inode[2][ncnt] = 0;
	inode[3][ncnt] = 0;
	inode[5][ncnt] = me->num_obs;
	ncnt++;
    } else {
	inode[1][ncnt] = scnt + 1;      /* S has 1 based, not 0 based
					 * subscripts */
	i = 0;
	for (spl = me->primary; spl; spl = spl->nextsplit) {
	    i++;
	    j = spl->var_num;
	    dsplit[0][scnt] = spl->improve;
	    if (numcat[j] == 0) {
		dsplit[1][scnt] = spl->spoint;
		isplit[2][scnt] = spl->csplit[0];
	    } else {            /* categorical */
		dsplit[1][scnt] = ccnt + 1;
		isplit[2][scnt] = numcat[j];
		for (k = 0; k < numcat[j]; k++)
		    csplit[k][ccnt] = spl->csplit[k];
		ccnt++;
	    }
	    isplit[0][scnt] = j + 1;    /* use "1" based subscripts */
	    isplit[1][scnt] = spl->count;
	    scnt++;
	}
	inode[2][ncnt] = i;

	i = 0;
	for (spl = me->surrogate; spl; spl = spl->nextsplit) {
	    i++;
	    j = spl->var_num;
	    dsplit[0][scnt] = spl->improve;
	    dsplit[2][scnt] = spl->adj;
	    if (numcat[j] == 0) {
		dsplit[1][scnt] = spl->spoint;
		isplit[2][scnt] = spl->csplit[0];
	    } else {
		dsplit[1][scnt] = ccnt + 1;
		isplit[2][scnt] = numcat[j];
		for (k = 0; k < numcat[j]; k++)
		    csplit[k][ccnt] = spl->csplit[k];
		ccnt++;
	    }
	    isplit[0][scnt] = j + 1;
	    isplit[1][scnt] = spl->count;
	    scnt++;
	}
	inode[3][ncnt] = i;
	inode[5][ncnt] = me->num_obs -
	    ((me->leftson)->num_obs + (me->rightson)->num_obs);

	ncnt++;

	ctmatrix(me->leftson, numcat,
		 dsplit, isplit, csplit, dnode, inode, 2 * id);
	ctmatrix(me->rightson, numcat,
		 dsplit, isplit, csplit, dnode, inode, 2 * id + 1);
    }
}
