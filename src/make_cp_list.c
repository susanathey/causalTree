/*
 * This routine creates the list of unique complexity parameters.
 * The list is maintained in sorted order.  If two parameters are within
 * "cplist_epsilon" of one another, then only the larger of them is
 * retained.
 *       CHANGE: 7/2000, the "cplist-epsilon" logic moved to S code
 *   I want the list sorted with the largest cp at the top of the list, since
 * that is the order that the CP table will be printed in.
 *
 *    After the partition routine is done, each node is labeled with the
 * complexity parameter appropriate if that node were the top of the tree.
 * However, if there is a more vunerable node further up, the node in
 * question will actually have the smaller complexity parameter; it will
 * be removed when its parent collapses. So this routine also adjusts each
 * C.P. to = minimum(my C.P., parent's C.P.).
 *
 *   This routine is called at the top level by causalTree, after causalTree has
 * initialized the first member of the linked cp-list, set its number of
 * splits to zero, and its risk to that for no splits at all.  This routine
 * allocates and links in the rest of the cp-list.  The make_cp_table
 * routine then fills in the rest of the variables in the list.
 *
 *  node *me;         pointer to my node structure
 *  double parent;    complexity of my parent node
 *
 *   When it comes time to cross-validate, we fill in xrisk and xstd
 */
#include <math.h>
#include "causalTree.h"
#include "node.h"
#include "causalTreeproto.h"

void
make_cp_list(pNode me, double parent, CpTable cptable_head)
{
    double me_cp;
    CpTable cplist, cptemp = NULL;

    if (me->complexity > parent)
      me->complexity = parent;
    me_cp = me->complexity;
    if (me_cp < ct.alpha)
	    me_cp = ct.alpha;       /* table should go no lower */
    if (me->leftson) {
      make_cp_list(me->leftson, me_cp, cptable_head);
	    make_cp_list(me->rightson, me_cp, cptable_head);
    }
    if (me_cp < parent) {       /* if not, then it can't be unique */
      for (cplist = cptable_head; cplist; cplist = cplist->forward) {
	   /* am I tied? */
	    if (me_cp == cplist->cp) /* exact ties */
        return;         

	    if (me_cp > cplist->cp)
		    break;
	    cptemp = cplist;
	  }

       /* insert new stuff after cptemp */
       /* was CALLOC and not cleaned up */
	cplist = (CpTable) ALLOC(1, sizeof(cpTable));
	cplist->cp = me_cp;
	cplist->risk = cplist->xrisk = cplist->xstd = 0;
	cplist->nsplit = 0;
	cplist->back = cptemp;
	cplist->forward = cptemp->forward;
	if (cptemp->forward)
	    (cptemp->forward)->back = cplist;
	else
	    cptable_tail = cplist;
	cptemp->forward = cplist;
	ct.num_unique_cp++;
	return;
    }
}
