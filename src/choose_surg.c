/*
 * A particular split routine, optimized for the surrogate variable
 *  search.  The "goodness" of a split is the total weights of concordant
 *  observations between the surrogate and the primary split.
 *  Note that the CART folks use the %concordance, which factors missing
 *  values into the equations somewhat differently.
 *
 *  y is coded as  +1=left, -1=right, 0=missing
 *
*/
#include "causalTree.h"
#include "causalTreeproto.h"

void
choose_surg(int n1, int n2, int *y, double *x, int *order,
	    int ncat, double *agreement, double *split, int *csplit,
	    double tleft, double tright, double *adj)
{
    int *left = ct.left, *right = ct.right;
    double *lwt = ct.lwt, *rwt = ct.rwt;
    double llwt, lrwt, rrwt, rlwt;      /* sum of weights for each */
    double agree, majority, total_wt;
    int success = 0;  // set to 1 when something worthwhile is found

    /*
     * I enforce that at least 2 obs must go each way, to avoid having an
     *  uncorrelated surrogate beat the "null" surrogate too easily
     * Observations with 0 weight don't count in this total
     */
    if (ncat == 0) {            /* continuous case */
	/*
	 * ll = y's that go left that are also sent left by my split
	 * lr = y's that go left that I send right
	 * rl = y's that go right that I send to the left
	 * rr = y's that go right that I send to the right
	 *
	 * The agreement is max(ll+rr, lr+rl), if weights were = 1;
	 *   actually max(llwt + rrwt, lrwt + rlwt) / denominator
	 */
	double lastx = 0.0;
	int ll, lr, rr, rl;
	ll = rl = 0;
	llwt = 0;
	rlwt = 0;
	for (int i = n2 - 1; i >= n1; i--) {  /* start with me sending all to the left */
	    int j = order[i];
	    if (j >= 0) {
		lastx = x[j];   /*this is why I run the loop backwards */
		switch (y[j]) {
		case LEFT:
		    if (ct.wt[j] > 0)
			ll++;
		    llwt += ct.wt[j];
		    break;
		case RIGHT:
		    if (ct.wt[j] > 0)
			rl++;
		    rlwt += ct.wt[j];
		    break;
		default:;
		}
	    }
	}

	if (llwt > rlwt)
	    agree = llwt;
	else
	    agree = rlwt;

	majority = agree;       /* worst possible agreement */
	total_wt = llwt + rlwt;

	lr = rr = 0;
	lrwt = 0;
	rrwt = 0;
	/*
	 *  March across, moving things from the right to the left
	 *    the "lastx" code is caring for ties in the x var
	 *    (The loop above sets it to the first unique x value).
	 */
	/* NB: might never set csplit[0] or split */
	csplit[0] = LEFT;
	*split = lastx; // a valid splitting value
	for (int i = n1; (ll + rl) >= 2; i++) {
	    int j = order[i];
	    if (j >= 0) {       /* not a missing value */
		if ((lr + rr) >= 2 && x[j] != lastx) {
		   /* new x found, evaluate the split */
		    if ((llwt + rrwt) > agree) {
			success = 1;
			agree = llwt + rrwt;
			csplit[0] = RIGHT;      /* < goes to the right */
			*split = (x[j] + lastx) / 2;
		    } else if ((lrwt + rlwt) > agree) {
			success = 1;
			agree = lrwt + rlwt;
			csplit[0] = LEFT;
			*split = (x[j] + lastx) / 2;
		    }
		}

		switch (y[j]) { /* update numbers */
		case LEFT:
		    if (ct.wt[j] > 0) {
			ll--;
			lr++;
		    }
		    llwt -= ct.wt[j];
		    lrwt += ct.wt[j];
		    break;
		case RIGHT:
		    if (ct.wt[j] > 0) {
			rl--;
			rr++;
		    }
		    rlwt -= ct.wt[j];
		    rrwt += ct.wt[j];
		    break;
		default:;      /* ignore missing y's */
		}
		lastx = x[j];
	    }
	}
    } else {                      /* categorical predictor */
	int defdir;
	int lcount = 0, rcount = 0;
	for (int i = 0; i < ncat; i++) {
	    left[i] = 0;
	    right[i] = 0;
	    lwt[i] = 0;
	    rwt[i] = 0;
	}

	/* First step:
	 *  left = table(x[y goes left]), right= table(x[y goes right])
	 *  so left[2] will be the number of x == 2's that went left,
	 *  and lwt[2] the sum of the weights for those observations.
	 * Only those with weight > 0 count in the totals
	 */
	for (int i = n1; i < n2; i++) {
	    int j = order[i];
	    if (j >= 0) {
		int k = (int) x[j] - 1;
		switch (y[j]) {
		case LEFT:
		    if (ct.wt[j] > 0)
			left[k]++;
		    lwt[k] += ct.wt[j];
		    break;
		case RIGHT:
		    if (ct.wt[j] > 0)
			right[k]++;
		    rwt[k] += ct.wt[j];
		    break;
		default:;
		}
	    }
	}

       /*
	*  Compute which is better: everyone to the right or left
	*/
	llwt = 0;
	rrwt = 0;
	for (int i = 0; i < ncat; i++) {
	    llwt += lwt[i];
	    rrwt += rwt[i];
	}
	if (llwt > rrwt) {
	    defdir = LEFT;
	    majority = llwt;
	} else {
	    defdir = RIGHT;
	    majority = rrwt;
	}
	total_wt = llwt + rrwt;

	/*
	 * We can calculate the best split category by category--- send each
	 *  x value individually to its better direction
	 */
	agree = 0.0;
	for (int i = 0; i < ncat; i++) {
	    if (left[i] == 0 && right[i] == 0)
		csplit[i] = 0;
	    else {
		if (lwt[i] < rwt[i] || (lwt[i] == rwt[i] && defdir == RIGHT)) {
		    agree += rwt[i];
		    csplit[i] = RIGHT;
		    lcount += left[i];
		    rcount += right[i];
		} else {
		    agree += lwt[i];
		    csplit[i] = LEFT;
		    lcount += right[i];
		    rcount += left[i];
		}
	    }
	}
	success = lcount > 1 && rcount > 1; /* sends at least 2 each way */
    }

    /*
     * success = 0 means no split was found that had at least 2 sent each
     *   way (not counting weights of zero), and for continuous splits also had
     *   an improvement in agreement.
     *  Due to round-off error such a split could still appear to have adj>0
     *  in the calculation futher below; avoid this.
     */

    if (!success) {
	*agreement = 0.0;
	*adj = 0.0;
	return;
    }

    /*
     *  Now we have the total agreement.  Calculate the %agreement and
     *    the adjusted agreement
     *  For both, do I use the total y vector as my denominator (my
     *    preference), or only the y's for non-missing x (CART book)?
     *    If the former, need to reset some totals.
     */
    if (ct.sur_agree == 0) {    /* use total table */
	total_wt = tleft + tright;
	if (tleft > tright)
	    majority = tleft;
	else
	    majority = tright;
    }

    *agreement = agree / total_wt;
    majority /= total_wt;
    *adj = (*agreement - majority) / (1. - majority);
}
