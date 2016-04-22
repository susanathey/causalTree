/*
 *  rescale an Survival time so that it is essentially exponential
 *   each interval between deaths has the same number of person-years
 *
 *	n	number of observations
 *	y	survival object, sorted by death times
 *	wt      vector of weights, sorted
 * output
 *      newy    rescaled time vector
 * scratch
 *      wtemp(n)
 */
#include "causalTree.h"

void
causalTreeexp(int *n2, double *y, double *wt, double *newy, double *wtemp)
{
    int n;
    double *stop, *event;
    int i, j;
    double tsum, dsum;          /* weighted sums of times and deaths */
    double time, ltime, rtime;
    double temp;
    double psum, scale;
    int last;

    n = *n2;
    stop = y;
    event = y + n;

    temp = 0;
    for (i = n - 1; i >= 0; i--) {
	temp += wt[i];
	wtemp[i] = temp;        /* sum of weights */
    }

    last = 0;
    ltime = 0;                  /* last event time */
    rtime = 0;                  /* rescaled time, cumulative */
    while (last < n) {
       /*
	* look ahead to find the next death
	*/
	psum = 0;
	for (i = last; i < n && event[i] == 0; i++)
	    psum += wt[i] * (stop[i] - ltime);  /* partial intervals */

       /*
	* Found it (or the end of the data)
	*/
	if (i > n) {            /* no more deaths */
	    for (i = last; i < n; i++)
		newy[i] = rtime;
	    last = n;
	} else {                /* rescale this interval */
	   /*
	    * count up the sum of the weighted deaths
	    */
	    dsum = 0;
	    time = stop[i];
	    for (; i < n && event[i] == 1 && stop[i] == time; i++)
		dsum += wt[i];  /* tied deaths */

	    tsum = (wtemp[i] + dsum) * (time - ltime) + psum;
	    scale = dsum / tsum;        /* scaling factor */

	    for (j = last; j < i; j++)
		newy[j] = rtime + (stop[j] - ltime) * scale;
	    rtime += (time - ltime) * scale;
	    last = i;
	    ltime = time;
	}
    }
}
