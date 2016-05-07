/*
 * quick sort routine : sort a vector of floats, and carry along an int
 *
 *  x:     vector to sort on
 *  start: first element of x to sort
 *  stop:  last element of x to sort
 *  cvec:  a vector to carry along
 */
#include "causalTree.h"
#include "causalTreeproto.h"

void
mysort(int start, int stop, double *x, int *cvec)
{
    int i, j, k;
    double temp, median;
    int tempd;

    while (start < stop) {
       /*
    	* first-- if the list is short, do an ordinary insertion sort
	    */
	if ((stop - start) < 11) {
	    for (i = start + 1; i <= stop; i++) {
		temp = x[i];
		tempd = cvec[i];
		j = i - 1;

		while (j >= start && x[j] > temp) {
		    x[j + 1] = x[j];
		    cvec[j + 1] = cvec[j];
		    j--;
		}
		x[j + 1] = temp;
		cvec[j + 1] = tempd;
	    }
	    return;
	}
	/*
	 * list is longer -- split it into two
	 *  I use the median of 3 values as the split point
	 */
	i = start;
	j = stop;
	k = (start + stop) / 2;

	median = x[k];
	if (x[i] >= x[k]) {     /* one of j or k is smallest */
	    if (x[j] > x[k]) {  /* k is smallest */
		if (x[i] > x[j])
		    median = x[j];
		else
		    median = x[i];
	    }
	} else {
	    if (x[j] < x[k]) {
		if (x[i] > x[j])
		    median = x[i];
		else
		    median = x[j];
	    }
	}

	/*
	 *  Now actually do the partitioning
	 *   Because we must have at least one element >= median, "i"
	 *   will never run over the end of the array.  Similar logic
	 *   applies to j.
	 * A note on the use of "<" rather than "<=".  If a list has lots
	 *   of identical elements, e.g. 80/100 are "3.5", then we will
	 *   often go to the swap step with x[i]=x[j]=median.  But we will
	 *   get the pointers i and j to meet approximately in the middle of
	 *   the list, and that is THE important condition for speed in a
	 *   quicksort.
	 *
	 */
	while (i < j) {
	   /*
	    * top pointer down till it points at something too large
	    */
	    while (x[i] < median)
		i++;

	   /*
	    * bottom pointer up until it points at something too small
	    */
	    while (x[j] > median)
		j--;

	    if (i < j) {
		if (x[i] > x[j]) {      /* swap */
		    temp = x[i];
		    x[i] = x[j];
		    x[j] = temp;
		    tempd = cvec[i];
		    cvec[i] = cvec[j];
		    cvec[j] = tempd;
		}
		i++;
		j--;
	    }
	}

	/*
	 * The while() step helps if there are lots of ties.  It will break
	 *  the list into 3 parts: < median, ==median, >=median, of which only
	 *  the top and bottom ones need further attention.
	 * The ">=" is needed because i may be  == to j
	 */
	while (x[i] >= median && i > start)
	    i--;
	while (x[j] <= median && j < stop)
	    j++;

	/*
	 * list has been split, now do a recursive call
	 *   always recur on the shorter list, as this keeps the total
	 *       depth of nested calls to less than log_base2(n).
	 */
	if ((i - start) < (stop - j)) { /* top list is shorter */
	    if ((i - start) > 0)
		mysort(start, i, x, cvec);
	    start = j;
	} else {                /* bottom list is shorter */
	    if ((stop - j) > 0)
		mysort(j, stop, x, cvec);
	    stop = i;
	}
    }
}
