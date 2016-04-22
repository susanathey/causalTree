/*
 * count up the number of nodes and splits in the final result
 *
 *  Gather the counts for myself, add in those of my children, and
 *    pass the total back to my parent
 */
#include "causalTree.h"
#include "node.h"
#include "causalTreeproto.h"

void
ctcountup(pNode me, int *nnode, int *nsplit, int *ncat)
{
    int node2, split2;
    int cat2;
    int i, j, k;
    pSplit ss;

    if (me->complexity <= ct.alpha || me->leftson == 0) {       /* no kids */
	*nnode = 1;
	*nsplit = 0;
	*ncat = 0;
    } else {
	i = 0;
	j = 0;
	k = 0;
	for (ss = me->primary; ss; ss = ss->nextsplit) {
	    i++;
	    if (ct.numcat[ss->var_num] > 0)
		k++;
	}
	for (ss = me->surrogate; ss; ss = ss->nextsplit) {
	    j++;
	    if (ct.numcat[ss->var_num] > 0)
		k++;
	}

	ctcountup(me->leftson, nnode, nsplit, ncat);
	ctcountup(me->rightson, &node2, &split2, &cat2);
	*nnode += 1 + node2;
	*nsplit += i + j + split2;
	*ncat += k + cat2;
    }
}
