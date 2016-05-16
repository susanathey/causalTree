/*
 * The discrte version for user (set temporarily as CTD)
 */

#include "causalTree.h"
#include "causalTreeproto.h"

#ifndef min
#define min(a,b) (((a) < (b)) ? (a) : (b))
#endif

#ifndef max
#define max(a,b)  (((a) > (b)) ? (a) : (b))
#endif

static double *sums, *wtsums, *treatment_effect;
static double *wts, *trs, *trsums;
static int *countn;
static int *tsplit;
static double *wtsqrsums, *wttrsqrsums;

// for discrete version:
static int *n_bucket, *n_tr_bucket, *n_con_bucket;
static double *wts_bucket, *trs_bucket;
static double *wtsums_bucket, *trsums_bucket;
static double *wtsqrsums_bucket, *trsqrsums_bucket; 
static double *tr_end_bucket, *con_end_bucket;


int
userDinit(int n, double *y[], int maxcat, char **error,
		int *size, int who, double *wt, double *treatment, 
		int bucketnum, int bucketMax, double *train_to_est_ratio)
{
	if (who == 1 && maxcat > 0) {
		graycode_init0(maxcat);
		countn = (int *) ALLOC(2 * maxcat, sizeof(int));
		tsplit = countn + maxcat;
		treatment_effect = (double *) ALLOC(8 * maxcat, sizeof(double));
		wts = treatment_effect + maxcat;
		trs = wts + maxcat;
		sums = trs + maxcat;
		wtsums = sums + maxcat;
		trsums = wtsums + maxcat;
		wtsqrsums = trsums + maxcat;
		wttrsqrsums = wtsqrsums + maxcat;
	}
	*size = 1;
	*train_to_est_ratio = n * 1.0 / ct.NumHonest;
	if (bucketnum == 0) 
		Rprintf("ERROR for buket!\n");
	return 0;
}


void
userDss(int n, double *y[], double *value, double *con_mean, double *tr_mean, double *risk, double *wt, double *treatment, 
		double max_y, double alpha, double train_to_est_ratio)
{
	int i;
	double temp0 = 0., temp1 = 0., twt = 0.; /* sum of the weights */ 
	double ttreat = 0.;
	double effect;
	double tr_var, con_var;
	double con_sqr_sum = 0., tr_sqr_sum = 0.;
	for (i = 0; i < n; i++) {
		temp1 += *y[i] * wt[i] * treatment[i];
		temp0 += *y[i] * wt[i] * (1 - treatment[i]);
		twt += wt[i];
		ttreat += wt[i] * treatment[i];
		tr_sqr_sum += (*y[i]) * (*y[i]) * wt[i] * treatment[i];
		con_sqr_sum += (*y[i]) * (*y[i]) * wt[i] * (1- treatment[i]);
	}

	effect = temp1 / ttreat - temp0 / (twt - ttreat);
	tr_var = tr_sqr_sum / ttreat - temp1 * temp1 / (ttreat * ttreat);
	con_var = con_sqr_sum / (twt - ttreat) - temp0 * temp0 / ((twt - ttreat) * (twt - ttreat));

	*tr_mean = temp1 / ttreat;
	*con_mean = temp0 / (twt - ttreat);
	*value = effect;
	*risk = 4 * twt * max_y * max_y - alpha * twt * effect * effect + 
		(1 - alpha) * (1 + train_to_est_ratio) * twt * (tr_var /ttreat  + con_var / (twt - ttreat));
}

void
userD(int n, double *y[], double *x, int nclass,
		int edge, double *improve, double *split, int *csplit,
		double myrisk, double *wt, double *treatment, int minsize, 
		double alpha, int bucketnum, int bucketMax, double train_to_est_ratio)
{
	int i, j;
	double temp;
	double left_sum, right_sum;
	double left_tr_sum, right_tr_sum;
	double left_tr_sqr_sum, right_tr_sqr_sum;
	double left_sqr_sum, right_sqr_sum;
	double tr_var, con_var;
	double left_tr_var, left_con_var, right_tr_var, right_con_var;
	double left_tr, right_tr;
	double left_wt, right_wt;
	int left_n, right_n;
	double best;
	int direction = LEFT;
	int where = 0;
	double node_effect, left_effect, right_effect;
	double left_temp, right_temp;
	int min_node_size = minsize;
	int bucketTmp;
	double trsum = 0.;
	int Numbuckets;

	double *cum_wt, *tmp_wt, *fake_x;
	double tr_wt_sum, con_wt_sum, con_cum_wt, tr_cum_wt;

	// for overlap:
	double tr_min, tr_max, con_min, con_max;
	double left_bd, right_bd;
	double cut_point;

	right_wt = 0;
	right_tr = 0;
	right_sum = 0;
	right_tr_sum = 0;
	right_sqr_sum = 0;
	right_tr_sqr_sum = 0;
	right_n = n;
	for (i = 0; i < n; i++) {
		right_wt += wt[i];
		right_tr += wt[i] * treatment[i];
		right_sum += *y[i] * wt[i];
		right_tr_sum += *y[i] * wt[i] * treatment[i];
		right_sqr_sum += (*y[i]) * (*y[i]) * wt[i];
		right_tr_sqr_sum += (*y[i]) * (*y[i]) * wt[i] * treatment[i];
		trsum += treatment[i];
	}

	temp = right_tr_sum / right_tr - (right_sum - right_tr_sum) / (right_wt - right_tr);
	tr_var = right_tr_sqr_sum / right_tr - right_tr_sum * right_tr_sum / (right_tr * right_tr);
	con_var = (right_sqr_sum - right_tr_sqr_sum) / (right_wt - right_tr)
		- (right_sum - right_tr_sum) * (right_sum - right_tr_sum) 
		/ ((right_wt - right_tr) * (right_wt - right_tr));
	node_effect = alpha * temp * temp * right_wt - (1 - alpha) * (1 + train_to_est_ratio) 
		* right_wt * (tr_var / right_tr  + con_var / (right_wt - right_tr));

	if (nclass == 0) {
		/* continuous predictor */
		cum_wt = (double *) ALLOC(n, sizeof(double));
		tmp_wt = (double *) ALLOC(n, sizeof(double));
		fake_x = (double *) ALLOC(n, sizeof(double));

		tr_wt_sum = 0.;
		con_wt_sum = 0.;
		con_cum_wt = 0.;
		tr_cum_wt = 0.;  

		// find the abs max and min of x:
		double max_abs_tmp = fabs(x[0]);
		for (i = 0; i < n; i++) {
			if (max_abs_tmp < fabs(x[i])) {
				max_abs_tmp = fabs(x[i]);
			}
		}

		// set tr_min, con_min, tr_max, con_max to a large/small value
		tr_min = max_abs_tmp;
		tr_max = -max_abs_tmp;
		con_min = max_abs_tmp;
		con_max = -max_abs_tmp;

		for (i = 0; i < n; i++) {
			if (treatment[i] == 0) {
				con_wt_sum += wt[i];
				if (con_min > x[i]) {
					con_min = x[i];
				}
				if (con_max < x[i]) {
					con_max = x[i];
				}
			} else {
				tr_wt_sum += wt[i];
				if (tr_min > x[i]) {
					tr_min = x[i];
				}
				if (tr_max < x[i]) {
					tr_max = x[i];
				}
			}
			cum_wt[i] = 0.;
			tmp_wt[i] = 0.;
			fake_x[i] = 0.;
		}

		// compute the left bound and right bound
		left_bd = max(tr_min, con_min);
		right_bd = min(tr_max, con_max);

		bucketTmp = min(round(trsum / (double)bucketnum), round(((double)n - trsum) / (double)bucketnum));
		Numbuckets = max(minsize, min(bucketTmp, bucketMax));

		for (i = 0; i < n; i++) {
			if (treatment[i] == 0) {
				tmp_wt[i] = wt[i] / con_wt_sum;     
				con_cum_wt += tmp_wt[i];
				cum_wt[i] = con_cum_wt;
				fake_x[i] = (int)floor(Numbuckets * cum_wt[i]);
			} else {
				tmp_wt[i] = wt[i] / tr_wt_sum;
				tr_cum_wt += tmp_wt[i];
				cum_wt[i] = tr_cum_wt;
				fake_x[i] = (int)floor(Numbuckets * cum_wt[i]);
			}
		}

		n_bucket = (int *) ALLOC(Numbuckets + 1,  sizeof(int));
		n_tr_bucket = (int *) ALLOC(Numbuckets + 1, sizeof(int));
		n_con_bucket = (int *) ALLOC(Numbuckets + 1, sizeof(int));
		wts_bucket = (double *) ALLOC(Numbuckets + 1, sizeof(double));
		trs_bucket = (double *) ALLOC(Numbuckets + 1, sizeof(double));
		wtsums_bucket = (double *) ALLOC(Numbuckets + 1, sizeof(double));
		trsums_bucket = (double *) ALLOC(Numbuckets + 1, sizeof(double));
		wtsqrsums_bucket = (double *) ALLOC(Numbuckets + 1, sizeof(double));
		trsqrsums_bucket = (double *) ALLOC(Numbuckets + 1, sizeof(double));
		tr_end_bucket = (double *) ALLOC(Numbuckets + 1, sizeof(double));
		con_end_bucket = (double *) ALLOC (Numbuckets + 1, sizeof(double));


		for (j = 0; j < Numbuckets + 1; j++) {
			n_bucket[j] = 0;
			n_tr_bucket[j] = 0;
			n_con_bucket[j] = 0;
			wts_bucket[j] = 0.;
			trs_bucket[j] = 0.;
			wtsums_bucket[j] = 0.;
			trsums_bucket[j] = 0.;
			wtsqrsums_bucket[j] = 0.;
			trsqrsums_bucket[j] = 0.;
		}

		for (i = 0; i < n; i++) {
			j = fake_x[i];
			n_bucket[j]++;
			wts_bucket[j] += wt[i];
			trs_bucket[j] += wt[i] * treatment[i];
			wtsums_bucket[j] += *y[i] * wt[i];
			trsums_bucket[j] += *y[i] * wt[i] * treatment[i];
			wtsqrsums_bucket[j] += (*y[i]) * (*y[i]) * wt[i];
			trsqrsums_bucket[j] += (*y[i]) * (*y[i]) * wt[i] * treatment[i];
			if (treatment[i] == 1) {
				tr_end_bucket[j] = x[i];
			} else {
				con_end_bucket[j] = x[i];
			}
		}

		left_wt = 0;
		left_tr = 0;
		left_n = 0;
		left_sum = 0;
		left_tr_sum = 0;
		left_sqr_sum = 0;
		left_tr_sqr_sum = 0;
		left_temp = 0.;
		right_temp = 0.;

		best = 0;

		for (j = 0; j < Numbuckets; j++) {
			left_n += n_bucket[j];
			right_n -= n_bucket[j];
			left_wt += wts_bucket[j];
			right_wt -= wts_bucket[j];
			left_tr += trs_bucket[j]; 
			right_tr -= trs_bucket[j];

			left_sum += wtsums_bucket[j];
			right_sum -= wtsums_bucket[j];

			left_tr_sum += trsums_bucket[j];
			right_tr_sum -= trsums_bucket[j];

			left_sqr_sum += wtsqrsums_bucket[j];
			right_sqr_sum -= wtsqrsums_bucket[j];

			left_tr_sqr_sum += trsqrsums_bucket[j];
			right_tr_sqr_sum -= trsqrsums_bucket[j];

			cut_point = (tr_end_bucket[j] + con_end_bucket[j]) / 2.0;
			if (left_n >= edge && right_n >= edge &&
					(int) left_tr >= min_node_size &&
					(int) left_wt - (int) left_tr >= min_node_size &&
					(int) right_tr >= min_node_size &&
					(int) right_wt - (int) right_tr >= min_node_size &&
					cut_point < right_bd && cut_point > left_bd) {
				left_temp = left_tr_sum / left_tr - 
					(left_sum - left_tr_sum) / (left_wt - left_tr);
				left_tr_var = left_tr_sqr_sum / left_tr - 
					left_tr_sum  * left_tr_sum / (left_tr * left_tr);
				left_con_var = (left_sqr_sum - left_tr_sqr_sum) / (left_wt - left_tr)  
					- (left_sum - left_tr_sum) * (left_sum - left_tr_sum)
					/ ((left_wt - left_tr) * (left_wt - left_tr));        
				left_effect = alpha * left_temp * left_temp * left_wt
					- (1 - alpha) * (1 + train_to_est_ratio) * left_wt 
					* (left_tr_var / left_tr + left_con_var / (left_wt - left_tr));

				right_temp = right_tr_sum / right_tr -
					(right_sum - right_tr_sum) / (right_wt - right_tr);
				right_tr_var = right_tr_sqr_sum / right_tr -
					right_tr_sum * right_tr_sum / (right_tr * right_tr);
				right_con_var = (right_sqr_sum - right_tr_sqr_sum) / (right_wt - right_tr)
					- (right_sum - right_tr_sum) * (right_sum - right_tr_sum) 
					/ ((right_wt - right_tr) * (right_wt - right_tr));
				right_effect = alpha * right_temp * right_temp * right_wt
					- (1 - alpha) * (1 + train_to_est_ratio) * right_wt 
					* (right_tr_var / right_tr + right_con_var / (right_wt - right_tr));
				temp = left_effect + right_effect - node_effect;
				if (temp > best) {
					best = temp;                  
					where = j; 
					if (left_temp < right_temp)
						direction = LEFT;
					else
						direction = RIGHT;
				}
			}
		}

		*improve = best;
		if (best > 0) {         /* found something */
			csplit[0] = direction;
			*split = (tr_end_bucket[where] + con_end_bucket[where]) / 2.0;
		}

	} else {
		/*
		 * Categorical predictor
		 */
		for (i = 0; i < nclass; i++) {
			countn[i] = 0;
			wts[i] = 0;
			trs[i] = 0;
			sums[i] = 0;
			wtsums[i] = 0;
			trsums[i] = 0;
			wtsqrsums[i] = 0;
			wttrsqrsums[i] = 0;
		}


		for (i = 0; i < n; i++) {
			j = (int) x[i] - 1;
			countn[j]++;
			wts[j] += wt[i];
			trs[j] += wt[i] * treatment[i];
			sums[j] += *y[i];
			wtsums[j] += *y[i] * wt[i];
			trsums[j] += *y[i] * wt[i] * treatment[i];
			wtsqrsums[j] += (*y[i]) * (*y[i]) * wt[i];
			wttrsqrsums[j] += (*y[i]) * (*y[i]) * wt[i] * treatment[i];
		}

		for (i = 0; i < nclass; i++) {
			if (countn[i] > 0) {
				tsplit[i] = RIGHT;
				treatment_effect[i] = trsums[j] / trs[j] - (wtsums[j] - trsums[j]) / (wts[j] - trs[j]);
			} else
				tsplit[i] = 0;
		}
		graycode_init2(nclass, countn, treatment_effect);

		/*
		 * Now find the split that we want
		 */
		left_wt = 0;
		left_tr = 0;
		left_n = 0;
		left_sum = 0;
		left_tr_sum = 0;
		left_sqr_sum = 0;
		left_tr_sqr_sum = 0;

		best = 0;
		where = 0;
		while ((j = graycode()) < nclass) {
			tsplit[j] = LEFT;
			left_n += countn[j];
			right_n -= countn[j];

			left_wt += wts[j];
			right_wt -= wts[j];

			left_tr += trs[j];
			right_tr -= trs[j];

			left_sum += wtsums[j];
			right_sum -= wtsums[j];

			left_tr_sum += trsums[j];
			right_tr_sum -= trsums[j];

			left_sqr_sum += wtsqrsums[j];
			right_sqr_sum -= wtsqrsums[j];

			left_tr_sqr_sum += wttrsqrsums[j];
			right_tr_sqr_sum -= wttrsqrsums[j];


			if (left_n >= edge && right_n >= edge &&
					(int) left_tr >= min_node_size &&
					(int) left_wt - (int) left_tr >= min_node_size &&
					(int) right_tr >= min_node_size &&
					(int) right_wt - (int) right_tr >= min_node_size) {

				left_temp = left_tr_sum / left_tr - (left_sum - left_tr_sum) / (left_wt - left_tr);
				left_tr_var = left_tr_sqr_sum / left_tr - left_tr_sum  * left_tr_sum / (left_tr * left_tr);
				left_con_var = (left_sqr_sum - left_tr_sqr_sum) / (left_wt - left_tr) 
					- (left_sum - left_tr_sum) * (left_sum - left_tr_sum)
					/ ((left_wt - left_tr) * (left_wt - left_tr));        
				left_effect = alpha * left_temp * left_temp * left_wt
					- (1 - alpha) * (1 + train_to_est_ratio) * left_wt 
					* (left_tr_var / left_tr + left_con_var / (left_wt - left_tr));

				right_temp = right_tr_sum / right_tr - (right_sum - right_tr_sum) / (right_wt - right_tr);
				right_tr_var = right_tr_sqr_sum / right_tr - 
					right_tr_sum * right_tr_sum / (right_tr * right_tr);
				right_con_var = (right_sqr_sum - right_tr_sqr_sum) / (right_wt - right_tr)
					- (right_sum - right_tr_sum) * (right_sum - right_tr_sum) 
					/ ((right_wt - right_tr) * (right_wt - right_tr));
				right_effect = alpha * right_temp * right_temp * right_wt
					- (1 - alpha) * (1 + train_to_est_ratio) * right_wt 
					* (right_tr_var / right_tr + right_con_var / (right_wt - right_tr)); 

				temp = left_effect + right_effect - node_effect;

				if (temp > best) {
					best = temp;
					if (left_temp > right_temp)
						for (i = 0; i < nclass; i++) csplit[i] = -tsplit[i];
					else
						for (i = 0; i < nclass; i++) csplit[i] = tsplit[i];
				}
			}
		}
		*improve = best;
	}
}

double userDpred(double *y, double wt, double treatment, double *yhat, double propensity) // pass in ct.which
{
	double ystar;
	double temp;

	ystar = y[0] * (treatment - propensity) / (propensity * (1 - propensity));
	temp = ystar - *yhat;
	return temp * temp * wt;
}
