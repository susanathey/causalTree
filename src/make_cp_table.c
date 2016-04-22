/*
 * Given a cptable list already initialized with the unique cp's in it,
 *  fill in the columns for risk and number of splits.
 *
 * The basic logic is : for each terminal node on the tree, start myself
 *  down at the bottom of the list of complexity parameters.  For each
 *  unique C.P. until my parent collapses, the node I'm in adds into that
 *  line of the CP table.  So walk up the CP list, adding in, until my
 *  parent would collapse; then report my position in the cp list to the
 *  parent and quit.
 *
 *  parent: complexity of the parent node
*/
#include "causalTree.h"
#include "node.h"
#include "causalTreeproto.h"

CpTable
make_cp_table(pNode me, double parent, int nsplit)
{
    CpTable cplist;

    if (me->leftson) {          /* if there are splits below */
	/*
	 * The 2 lines below are perhaps devious
	 *  1) Since the return value depends on ones parent, both calls will
	 *       return the same thing.
	 *  2) I send 0 to the left to keep the current split (me) from
	 *       being counted twice, once by each child.
	 */
	make_cp_table(me->leftson, me->complexity, 0);
	cplist = make_cp_table(me->rightson, me->complexity, nsplit + 1);
    } else
	cplist = cptable_tail;

    while (cplist->cp < parent) {
	cplist->risk += me->risk;
	cplist->nsplit += nsplit;
	cplist = cplist->back;
    }

    return cplist;
}
