
/*
 * When partition is done, each node is labeled with the complexity
 *  appropriate if it were the top of the tree.  Actually, the complexity
 *  should be min(me, any-node-above-me).  This routine fixes that.
 */
#include "causalTree.h"
#include "node.h"
#include "causalTreeproto.h"

void
fix_cp(pNode me, double parent_cp)
{
    if (me->complexity > parent_cp)
	me->complexity = parent_cp;

    if (me->leftson) {
	fix_cp(me->leftson, me->complexity);
	fix_cp(me->rightson, me->complexity);
    }
}
