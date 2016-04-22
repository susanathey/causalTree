/*
 * print out the tree, in all it's glory
 *
 * This routine exists in the sources only for debugging puctoses -
 *   you will see occasional commented out calls to it in the code
 *
 * It makes for the nicest printout if the tree is printed out by depth,
 *  rather than by the usual recursive search.  That is, print it out the
 *  way it would be drawn: top node, then the two children of the top node,
 *  then the four grandchildren of the top node, etc.  In order to do
 *  this we run up and down the tree like crazy.
 *
 * maxdepth = max depth to print, 1= top node only, etc.
 */
#include <stdio.h>
#include "node.h"
#include "causalTree.h"

static void printme(pNode me, int id);
static void print_tree2(pNode me, int id, int mydepth, int target);

void
print_tree(pNode me, int maxdepth)
{
    int i;

    printme(me, 1);
    for (i = 2; i <= maxdepth; i++) {
	if (me->leftson)
	    print_tree2(me->leftson, 2, 2, i);
	if (me->rightson)
	    print_tree2(me->rightson, 3, 2, i);
    }
}

/*
 * Recursively run down and find children of the right depth
 */
static void
print_tree2(pNode me, int id, int mydepth, int target)
{
    if (mydepth == target)
	printme(me, id);
    else {
	if (me->leftson)
	    print_tree2(me->leftson, 2 * id, mydepth + 1, target);
	if (me->rightson)
	    print_tree2(me->rightson, 2 * id + 1, mydepth + 1, target);
    }
}

/* the actual work routine */
static void
printme(pNode me, int id)
{
    int i, j, k;
    pSplit ss;

    Rprintf("\n\nNode number %d: %d observations", id, me->num_obs);
    Rprintf("\t   Complexity param= %f\n", me->complexity);
    Rprintf("  response estimate=%f,  risk/n= %f\n", *(me->response_est),
	    me->risk / me->num_obs);

    if (me->leftson)
	Rprintf("  left son=%d (%d obs)", 2 * id, (me->leftson)->num_obs);
    if (me->rightson)
	Rprintf(" right son=%d (%d obs)", 2 * id + 1, (me->rightson)->num_obs);
    if (me->leftson && me->rightson) {
	i = me->num_obs - ((me->leftson)->num_obs + (me->rightson)->num_obs);
	if (i == 0)
	    Rprintf("\n");
	else
	    Rprintf(", %d obs do not split\n", i);
    } else
	Rprintf("\n");

    Rprintf("  Primary splits:\n");
    for (ss = me->primary; ss; ss = ss->nextsplit) {
	j = ss->var_num;
	if (ct.numcat[j] == 0) {
	    if (ss->csplit[0] == LEFT)
		Rprintf
		    ("\tvar%d < %5g to the left, improve=%5.3f,  (%d missing)\n",
		     j, ss->spoint, ss->improve, me->num_obs - ss->count);
	    else
		Rprintf
		    ("\tvar%d > %5g to the left, improve=%5.3f, (%d missing)\n",
		     j, ss->spoint, ss->improve, me->num_obs - ss->count);
	} else {
	    Rprintf("\tvar%d splits as ", j);
	    for (k = 0; k < ct.numcat[j]; k++) {
		switch (ss->csplit[k]) {
		case LEFT:
		    Rprintf("L");
		    break;
		case RIGHT:
		    Rprintf("R");
		    break;
		case 0:
		    Rprintf("-");
		}
	    }
	    if (ct.numcat[j] < 7)
		Rprintf(",\timprove=%5.3f, (%d missing)\n",
			ss->improve, (me->num_obs - ss->count));
	    else
		Rprintf(", improve=%5.3f, (%d missing)\n",
			ss->improve, (me->num_obs - ss->count));
	}
    }

   /*
    * Now print the surrogate splits.
    */
    if (me->surrogate)
	Rprintf("  Surrogate splits:\n");
    for (ss = me->surrogate; ss; ss = ss->nextsplit) {
	j = ss->var_num;
	if (ct.numcat[j] == 0) {
	    if (ss->csplit[0] == LEFT)
		Rprintf
		    ("\tvar%d < %5g to the left, agree=%5.3f, (%d split)\n",
		     j, ss->spoint, ss->improve, ss->count);
	    else
		Rprintf
		    ("\tvar%d > %5g to the left, agree=%5.3f, (%d split)\n",
		     j, ss->spoint, ss->improve, ss->count);
	} else {
	    Rprintf("\tvar%d splits as ", j);
	    for (k = 0; k < ct.numcat[j]; k++) {
		switch (ss->csplit[k]) {
		case LEFT:
		    Rprintf("L");
		    break;
		case RIGHT:
		    Rprintf("R");
		    break;
		case 0:
		    Rprintf("-");
		}
	    }
	    if (ct.numcat[j] < 7)
		Rprintf(",\tagree=%5.3f, (%d split)\n",
			ss->improve, ss->count);
	    else
		Rprintf(", agree=%5.3f, (%d split)\n",
			ss->improve, ss->count);
	}
    }
}
