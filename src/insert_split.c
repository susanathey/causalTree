/*
 * sort a new split into a linked list, based on its "improvement"
 *
 *  allocates new memory as needed
 *   returns 0 if the new element isn't good enough,
 *   the address of the new element otherwise
 */
#include "causalTree.h"
#include "node.h"
#include "causalTreeproto.h"

pSplit
insert_split(pSplit *listhead, int ncat, double improve, int max)
{
    int nlist;
    pSplit s1, s2, s3 = NULL, s4;

    // csplit[0] gets used even for continuous splits.
    if (ncat == 0) ncat = 1; 
    int splitsize = sizeof(Split) + (ncat - 20) * sizeof(int);

    // The Split structure is sized for 2 categpries.
    if (*listhead == 0) {
    /* first call to a new list */
	s3 = (pSplit) CALLOC(1, splitsize);
	s3->nextsplit = NULL;
	*listhead = s3;
	return s3;
    }
    if (max < 2) {
       /* user asked for only 1 to be retained! */
	s3 = *listhead;
	if (improve <= s3->improve)
	    return NULL;
	if (ncat > 1) {
	    Free(s3);
	    s3 = (pSplit) CALLOC(1, splitsize);
	    s3->nextsplit = NULL;
	    *listhead = s3;
	}
	return s3;
    }
   /* set up --- nlist = length of list, s4=last element, s3=next to last */
    nlist = 1;
    for (s4 = *listhead; s4->nextsplit; s4 = s4->nextsplit) {
	s3 = s4;
	nlist++;
    }

   /* now set up so that the "to be added" is between s1 and s2 */
    s1 = *listhead;
    for (s2 = *listhead; s2; s2 = s2->nextsplit) {
	if (improve > s2->improve)
	    break;
	s1 = s2;
    }

    if (nlist == max) {
	if (s2 == 0)
	    return NULL;        /* not good enough */
	if (ncat > 1) {
	   // FIXME: use Realloc
	    Free(s4);           /* get new memory -- this chunk may be too
				 * small */
	    s4 = (pSplit) CALLOC(1, splitsize);
	}
	if (s1 == s3)
	    s4->nextsplit = NULL;
	else {
	    s3->nextsplit = NULL;
	    s4->nextsplit = s2;
	}
    } else {
	s4 = (pSplit) CALLOC(1, splitsize);
	s4->nextsplit = s2;
    }
    if (s2 == *listhead)
	*listhead = s4;
    else
	s1->nextsplit = s4;
    return s4;
}
