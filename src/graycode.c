  /*
 * Walk through subsets in an ordered way.
 *   For all subsets, this is the classic gray code.
 *
 * Graycode_init0 is called once at the very beginning with
 *   the maximum number of categories.  It allocates a scratch vector.
 *
 * Graycode_init1 is called once for each unordered variable.
 * Graycode_init2 is called once for each orderable variable,
 *   the second argument is a vector that will be used to rank the variables
 */
#include "causalTree.h"
#include "causalTreeproto.h"
static int *gray;
static int maxc, gsave;

void
graycode_init0(int maxcat)
{
    gray = (int *) ALLOC(maxcat, sizeof(int));
}

void
graycode_init1(int numcat, int *count)
{
    int i;

    maxc = numcat;
    for (i = 0; i < maxc; i++)
	gray[i] = count[i] ? 1 : 0;
    gsave = -2;
}


void
graycode_init2(int numcat, int *count, double *val)
{
    int i, j, k;
    double temp;
    maxc = numcat;
    
    /*
     *   sort categories with no members first
     *   then order by val
     */
     
    gray[0] = 0;
    if (count[0] == 0)
	    k = 1;
    else
	    k = 0;
    for (i = 1; i < maxc; i++) {
      if (count[i] == 0) {
	      for (j = i - 1; j >= k; j--) {
		        gray[j + 1] = gray[j];
		        val[j + 1] = val[j];
	      }
	      gray[k++] = i;
	    } else {
        // bubble sort
	      temp = val[i];
	      for (j = i - 1; j >= k && val[j] > temp; j--) {
		      gray[j + 1] = gray[j];
		      val[j + 1] = val[j];
	      }
	      val[j + 1] = temp;
	      gray[j + 1] = i;
	    }
    }
    gsave = k - 1;
}

/*
* Everyone starts in the right hand group
* This routine returns the next subject who needs to
*  change allegiance.
* A value of maxc means that we're done.
*/
int
graycode(void)
{
    int i;

    if (gsave > -2) {           /* ordered data */
	gsave++;
	if (gsave < maxc)
	    return gray[gsave];
	else
	    return maxc;
    } else {
	/*
	 * Form next subgroup.  We do this using the classic Gray code.
	 *  The initial subset has everyone in the right group.  Each
	 *  subset varies from the prior by only one member -- the
	 *  following item changes groups: 1,2,1,4,1,2,1,8,1,2,1,4,1,...
	 *  The outer loop only goes up to maxc-1: we know for causalTree that
	 *    changing the allegiance of the last subject is never necessary
	 */
	for (i = 0; i < maxc - 1; i++) {
	    if (gray[i] == 1) {
		gray[i] = 2;
		return i;
	    } else if (gray[i] == 2)
		gray[i] = 1;
	}
	return maxc;            /* signal "done" */
    }
}
