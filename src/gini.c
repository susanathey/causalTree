/*
 * The routines for gini-classification
 */
#include <math.h>
#include "causalTree.h"
#include "causalTreeproto.h"

static int numclass;
static double *left,            /* left branch n (weighted) */
 *right, **ccnt;
static double *prior, *aprior,  /* altered priors */
 *freq, *loss;                  /* loss matrix */
static int *tsplit, *countn;
static double *awt, *rate;
static double (*impurity) ();

static double
gini_impure1(double p)
{
    return p * (1.0 - p);
}

static double
gini_impure2(double p)
{
    if (p == 0) return 0.0; else return -p * log(p);
}

int
giniinit(int n, double **y, int maxcat, char **error,
	 double *parm, int *size, int who, double *wt)
{
    int i, j, k;
    double temp;

   /* allocate memory  and setup losses */
    if (who == 1) {
	numclass = 0;           /* number of classes */
	for (i = 0; i < n; i++)
	    if (*y[i] > numclass)
		numclass = (int) *y[i];

	if (parm[numclass + numclass * numclass] == 2)
	    impurity = gini_impure2;
	else
	    impurity = gini_impure1;

	left = (double *) ALLOC(numclass * 2, sizeof(double));
	right = left + numclass;

	tsplit = (int *) ALLOC(maxcat * 2, sizeof(int));
	countn = tsplit + maxcat;

	awt = (double *) ALLOC(maxcat * 2, sizeof(double));
	rate = awt + maxcat;

	if (maxcat > 0) {
	    graycode_init0(maxcat);
	    ccnt = (double **) ALLOC(numclass, sizeof(double *));
	    ccnt[0] = (double *) ALLOC(numclass * maxcat, sizeof(double));
	    for (i = 1; i < numclass; i++)
		ccnt[i] = ccnt[i - 1] + maxcat;
	}
	i = 3 * numclass + numclass * numclass;
	prior = (double *) ALLOC(i, sizeof(double));
	aprior = prior + numclass;
	freq = aprior + numclass;
	loss = freq + numclass;

	for (i = 0; i < numclass; i++)
	    freq[i] = 0;
	temp = 0;
	for (i = 0; i < n; i++) {
	    j = (int) *y[i] - 1;
	    freq[j] += wt[i];
	    temp += wt[i];      /* sum total of weights */
	}
	for (i = 0; i < numclass; i++)
	    freq[i] /= temp;    /* relative frequency */

	temp = 0;
	for (i = 0; i < numclass; i++) {
	    prior[i] = parm[i];
	    aprior[i] = 0;
	    for (j = 0; j < numclass; j++) {
		k = numclass * j + i;
		loss[k] = parm[numclass + k];
		temp += loss[k] * prior[i];
		aprior[i] += loss[k] * prior[i];
	    }
	}
	for (i = 0; i < numclass; i++) {
	    if (freq[i] > 0) {  /* watch out for a missing class */
		prior[i] /= freq[i];
		aprior[i] /= (temp * freq[i]);  /* pi_i / n_i */
	    }
	}
    }
    *size = 2 + numclass;
    return 0;
}

/*
 * Compute the predicted response and the classification error
 *   This is R(T) (this node's contribution) in the paper.
 */
void
ginidev(int n, double **y, double *value, double *risk, double *wt)
{
    int i, j, max = 0;
    double temp, dev = 0;
    double prob;

   /*
    * count up number in each class,
    *   and P(T), the probability of reaching this branch of the tree
    */
    for (i = 0; i < numclass; i++)
	freq[i] = 0;
    temp = 0;
    for (i = 0; i < n; i++) {
	j = (int) y[i][0] - 1;
	freq[j] += wt[i];
	temp += wt[i] * prior[j];
    }
    prob = temp;                /* this is actually P(T)*n; R code will fix
				 * it up */

   /*
    * Now compute best class and its error
    */
    for (i = 0; i < numclass; i++) {    /* assume class i were the prediction */
	temp = 0;
	for (j = 0; j < numclass; j++)
	    temp += freq[j] * loss[i * numclass + j] * prior[j];
	if (i == 0 || temp < dev) {
	    max = i;
	    dev = temp;
	}
    }

    value[0] = max + 1;         /* remember: external groups start at 1 */
    for (i = 0; i < numclass; i++)
	value[i + 1] = freq[i];
    value[numclass + 1] = prob;
    *risk = dev;
}


/*
* return the error for a particular point
*/
double
ginipred(double *y, double *pred)
{
    int i, j;
    double temp;
    i = (int) y[0] - 1;
    j = (int) *pred - 1;
    temp = prior[i] * loss[i * numclass + j];
    return temp;
}


/*
 * The gini splitting function.  Find that split point in x such that
 *  the rss within the two groups is decreased as much
 *  as possible.
 */
void
gini(int n, double *y[], double *x, int numcat,
     int edge, double *improve, double *split, int *csplit, double my_risk,
     double *wt)
{
    int i, j, k;
    double lwt, rwt;
    int rtot, ltot;
    int direction = LEFT, where = 0;
    double total_ss, best, temp, p;
    double lmean, rmean;        /* used to decide direction */

    for (i = 0; i < numclass; i++) {
	left[i] = 0;
	right[i] = 0;
    }
    lwt = 0;
    rwt = 0;
    rtot = 0;
    ltot = 0;
    for (i = 0; i < n; i++) {
	j = (int) *y[i] - 1;
	rwt += aprior[j] * wt[i];  /* altered weight = prior * case_weight */
	right[j] += wt[i];
	rtot++;
    }
    total_ss = 0;
    for (i = 0; i < numclass; i++) {
	temp = aprior[i] * right[i] / rwt;      /* p(class=i, given node A) */
	total_ss += rwt * (*impurity) (temp);   /* p(A) * I(A) */
    }
    best = total_ss;  /* total weight of right * impurity of right + 0 *0 */

    /*
     * at this point we split into 2 disjoint paths
     */
    if (numcat > 0)
	goto categorical;

    for (i = 0; rtot > edge; i++) {
	j = (int) *y[i] - 1;
	rwt -= aprior[j] * wt[i];
	lwt += aprior[j] * wt[i];
	rtot--;
	ltot++;
	right[j] -= wt[i];
	left[j] += wt[i];

	if (ltot >= edge && x[i + 1] != x[i]) {
	    temp = 0;
	    lmean = 0;
	    rmean = 0;
	    for (j = 0; j < numclass; j++) {
		p = aprior[j] * left[j] / lwt;  /* p(j | left) */
		temp += lwt * (*impurity) (p);  /* p(left) * I(left) */
		lmean += p * j;
		p = aprior[j] * right[j] / rwt; /* p(j | right) */
		temp += rwt * (*impurity) (p);  /* p(right) * I(right) */
		rmean += p * j;
	    }
	    if (temp < best) {
		best = temp;
		where = i;
		direction = lmean < rmean ? LEFT : RIGHT;
	    }
	}
    }

    *improve = total_ss - best;
    if (*improve > 0) {         /* found something */
	csplit[0] = direction;
	*split = (x[where] + x[where + 1]) / 2;
    }
    return;

categorical:;
    /*
     * First collapse the data into a numclass x numcat array
     *  ccnt[i][j] = number of class i obs, category j of the predictor
     */
    for (j = 0; j < numcat; j++) {
	awt[j] = 0;
	countn[j] = 0;
	for (i = 0; i < numclass; i++)
	    ccnt[i][j] = 0;
    }
    for (i = 0; i < n; i++) {
	j = (int) *y[i] - 1;
	k = (int) x[i] - 1;
	awt[k] += aprior[j] * wt[i];
	countn[k]++;
	ccnt[j][k] += wt[i];
    }

    for (i = 0; i < numcat; i++) {
	if (awt[i] == 0)
	    tsplit[i] = 0;
	else {
	    rate[i] = ccnt[0][i] / awt[i];      /* a scratch array */
	    tsplit[i] = RIGHT;
	}
    }

    if (numclass == 2)
	graycode_init2(numcat, countn, rate);
    else
	graycode_init1(numcat, countn);

    while ((i = graycode()) < numcat) {
       /* item i changes groups */
	if (tsplit[i] == LEFT) {
	    tsplit[i] = RIGHT;
	    rwt += awt[i];
	    lwt -= awt[i];
	    rtot += countn[i];
	    ltot -= countn[i];
	    for (j = 0; j < numclass; j++) {
		right[j] += ccnt[j][i];
		left[j] -= ccnt[j][i];
	    }
	} else {
	    tsplit[i] = LEFT;
	    rwt -= awt[i];
	    lwt += awt[i];
	    rtot -= countn[i];
	    ltot += countn[i];
	    for (j = 0; j < numclass; j++) {
		right[j] -= ccnt[j][i];
		left[j] += ccnt[j][i];
	    }
	}

	if (ltot >= edge && rtot >= edge) {
	    temp = 0;
	    lmean = 0;
	    rmean = 0;
	    for (j = 0; j < numclass; j++) {
		p = aprior[j] * left[j] / lwt;
		temp += lwt * (*impurity) (p);
		lmean += p * j;
		p = aprior[j] * right[j] / rwt; /* p(j | right) */
		temp += rwt * (*impurity) (p);  /* p(right) * I(right) */
		rmean += p * j;
	    }
	    if (temp < best) {
		best = temp;
		if (lmean < rmean)
		    for (j = 0; j < numcat; j++) csplit[j] = tsplit[j];
		else
		    for (j = 0; j < numcat; j++) csplit[j] = -tsplit[j];
	    }
	}
    }
    *improve = total_ss - best;
}
