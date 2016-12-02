/*
 * split.Rule = user (set temporarily as CT)
 */
#include <math.h>
#include "causalTree.h"
#include "causalTreeproto.h"

static double *sums, *wtsums, *treatment_effect;
static double *wts, *trs, *trsums;
static int *countn;
static int *tsplit;
static double *wtsqrsums, *trsqrsums;

int
policyinit(int n, double *y[], int maxcat, char **error,
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
		trsqrsums = wtsqrsums + maxcat;
	}
	*size = 1;
	*train_to_est_ratio = n * 1.0 / ct.NumHonest;
	return 0;
}


void
policyss(int n, double *y[], double *value,  double *con_mean, double *tr_mean, 
		double *risk, double *wt, double *treatment, double max_y,
		double alpha, double train_to_est_ratio)
{
	int i,j;
	//double temp0 = 0., temp1 = 0., twt = 0.; /* sum of the weights */ 
	//double ttreat = 0.;
	//double effect;
	//double tr_var, con_var;
	//double con_sqr_sum = 0., tr_sqr_sum = 0.;
	
	double *temp0, *temp1, *twt; /* sum of the weights */ 
  double *ttreat;
  double *effect;
  double *tr_var, *con_var;
  double *con_sqr_sum, *tr_sqr_sum;
  double ntreats = ct.ntreats;

	//alloc space: some matrices?
	 temp0 = (double *) ALLOC(ntreats, sizeof(double));
	 temp1 = (double *) ALLOC(ntreats, sizeof(double));
	 twt = (double *) ALLOC(ntreats, sizeof(double));
	 ttreat = (double *) ALLOC(ntreats, sizeof(double));
	 tr_sqr_sum = (double *) ALLOC(ntreats, sizeof(double));
	 con_sqr_sum = (double *) ALLOC(ntreats, sizeof(double));
	 effect = (double *) ALLOC(ntreats, sizeof(double));
	 tr_var = (double *) ALLOC(ntreats, sizeof(double));
	 con_var = (double *) ALLOC(ntreats, sizeof(double));
	
	for (j = 0; j < ntreats; j++) {
	for (i = 0; i < n; i++) {
		//temp1 += *y[i] * wt[i] * treatment[i];
		//temp0 += *y[i] * wt[i] * (1 - treatment[i]);
		//twt += wt[i];
		//ttreat += wt[i] * treatment[i];
		//tr_sqr_sum += (*y[i]) * (*y[i]) * wt[i] * treatment[i];
		//con_sqr_sum += (*y[i]) * (*y[i]) * wt[i] * (1- treatment[i]);

		//all vars. below need to be vectors of length ntreats, and run a loop (1:ntreats) over this
		temp1[j] += *y[i] * wt[i] * (treatment[i]==j);
		temp0[j] += *y[i] * wt[i] * (1 - (treatment[i]==j));
		twt[j] += wt[i];
		ttreat[j] += wt[i] * (treatment[i]==j);
		tr_sqr_sum[j] += (*y[i]) * (*y[i]) * wt[i] * (treatment[i]==j);
		con_sqr_sum[j] += (*y[i]) * (*y[i]) * wt[i] * (1- (treatment[i]==j));
		
		
			}
	}
	//need a set of sums and effects here, one for each treatment (treatment, wt, y external vars.)

	//effect = temp1 / ttreat - temp0 / (twt - ttreat);
	//tr_var = tr_sqr_sum / ttreat - temp1 * temp1 / (ttreat * ttreat);
	//con_var = con_sqr_sum / (twt - ttreat) - temp0 * temp0 / ((twt - ttreat) * (twt - ttreat));
  
  for (j = 0; j < ntreats; j++) {
  effect[j] = temp1[j] / ttreat[j] - temp0[j] / (twt[j] - ttreat[j]);
  tr_var[j] = tr_sqr_sum[j] / ttreat[j] - temp1[j] * temp1[j] / (ttreat[j] * ttreat[j]);
  con_var[j] = con_sqr_sum[j] / (twt[j] - ttreat[j]) - temp0[j] * temp0[j] / ((twt[j] - ttreat[j]) * (twt[j] - ttreat[j]));
  }
  
	//returned values: tr_mean, con_mean, value, risk: need to be vectorized
	//also, handle top level calling functions to use vectors:
	//how to do this without affecting other options?
	//*tr_mean = temp1 / ttreat;
	//*con_mean = temp0 / (twt - ttreat);
	//*value = effect;
	//*risk = 4 * twt * max_y * max_y - alpha * twt * effect * effect + 
	//	(1 - alpha) * (1 + train_to_est_ratio) * twt * (tr_var /ttreat  + con_var / (twt - ttreat));

	for (j = 0; j < ntreats; j++)
	  {
	tr_mean[j] = temp1[j] / ttreat[j];
	con_mean[j] = temp0[j] / (twt[j] - ttreat[j]);
	value[j] = effect[j];
	//risk is risk_multi for the ct object
	risk[j] = 4 * twt[j] * max_y * max_y - alpha * twt[j] * effect[j] * effect[j] + 
	(1 - alpha) * (1 + train_to_est_ratio) * twt[j] * (tr_var[j] /ttreat[j]  + con_var[j] / (twt[j] - ttreat[j]));
	  }
	}


void policy(int n, double *y[], double *x, int nclass, int edge, double *improve, double *split, 
		int *csplit, double myrisk, double *wt, double *treatment, int minsize, double alpha,
		double train_to_est_ratio)
{
	int i, j;
	double *temp;
	double *left_sum, *right_sum;
	double *left_tr_sum, *right_tr_sum;
	double *left_tr, *right_tr;
	double *left_wt, *right_wt;
	int *left_n, *right_n;
	double *best;
	int *direction; // = LEFT;
	int *where; // = 0;
	double *node_effect, *left_effect, *right_effect;
	double *left_temp, *right_temp;
	int min_node_size = minsize;

	double *tr_var, *con_var;
	double *right_sqr_sum, *right_tr_sqr_sum, *left_sqr_sum, *left_tr_sqr_sum;
	double *left_tr_var, *left_con_var, *right_tr_var, *right_con_var;
  
  //alloc the above vectors
  temp = (double *) ALLOC(ntreats, sizeof(double));
  
  right_wt = (double *) ALLOC(ntreats, sizeof(double));
  right_tr = (double *) ALLOC(ntreats, sizeof(double));
  right_sum = (double *) ALLOC(ntreats, sizeof(double));
  right_tr_sum = (double *) ALLOC(ntreats, sizeof(double));
  right_tr_sqr_sum = (double *) ALLOC(ntreats, sizeof(double));
  right_tr_var = (double *) ALLOC(ntreats, sizeof(double));
  right_effect = (double *) ALLOC(ntreats, sizeof(double));
  right_n = (double *) ALLOC(ntreats, sizeof(double));
  right_con_var = (double *) ALLOC(ntreats, sizeof(double));
  
  left_wt = (double *) ALLOC(ntreats, sizeof(double));
  left_tr = (double *) ALLOC(ntreats, sizeof(double));
  left_sum = (double *) ALLOC(ntreats, sizeof(double));
  left_tr_sum = (double *) ALLOC(ntreats, sizeof(double));
  left_tr_sqr_sum = (double *) ALLOC(ntreats, sizeof(double));
  left_tr_var = (double *) ALLOC(ntreats, sizeof(double));
  left_effect = (double *) ALLOC(ntreats, sizeof(double));
  left_n = (double *) ALLOC(ntreats, sizeof(double));
  left_con_var = (double *) ALLOC(ntreats, sizeof(double));
  
  best = (double *) ALLOC(ntreats, sizeof(double));
  right_effect = (double *) ALLOC(ntreats, sizeof(double));
  left_effect = (double *) ALLOC(ntreats, sizeof(double));
  node_effect = (double *) ALLOC(ntreats, sizeof(double));
  direction = (double *) ALLOC(ntreats, sizeof(double));
  where = (double *) ALLOC(ntreats, sizeof(double));
  
  right_temp = (double *) ALLOC(ntreats, sizeof(double));
  left_temp = (double *) ALLOC(ntreats, sizeof(double));
  
  tr_var = (double *) ALLOC(ntreats, sizeof(double));
  con_var = (double *) ALLOC(ntreats, sizeof(double));
  
  double ntreats = ct.ntreats;
  
  //memsets, set direction to LEFT, rest to 0
  memset(right_wt, 0, ntreats);
  memset(right_tr, 0, ntreats);
  memset(right_sum, 0, ntreats);
  memset(right_tr_sum, 0, ntreats);
  memset(right_sqr_sum, 0, ntreats);
  memset(right_tr_sqr_sum, 0, ntreats);
  memset(right_n, n, ntreats);
  
  memset(direction, LEFT, ntreats);
  memset(where, 0, ntreats);
  
    /*
	right_wt = 0.;
	right_tr = 0.;
	right_sum = 0.;
	right_tr_sum = 0.;
	right_sqr_sum = 0.;
	right_tr_sqr_sum = 0.;
	right_n = n;
  */
  //run loop instead to set the above to 0: memset(tree, 0, nodesize);
  
  
	
	for (j = 0; j < ntreats; j++){
		for (i = 0; i < n; i++) {
		right_wt[j] += wt[i];
		right_tr[j] += wt[i] * (treatment[i]==j);
		right_sum[j] += *y[i] * wt[i];
		right_tr_sum[j] += *y[i] * wt[i] * (treatment[i]==j);
		right_sqr_sum[j] += (*y[i]) * (*y[i]) * wt[i];
		right_tr_sqr_sum[j] += (*y[i]) * (*y[i]) * wt[i] * (treatment[i]==j);
	}
	

	temp[j] = right_tr_sum[j] / right_tr[j] - (right_sum[j] - right_tr_sum[j]) / (right_wt[j] - right_tr[j]);
	tr_var[j] = right_tr_sqr_sum[j] / right_tr[j] - right_tr_sum[j] * right_tr_sum[j] / (right_tr[j] * right_tr[j]);
	con_var[j] = (right_sqr_sum[j] - right_tr_sqr_sum[j]) / (right_wt[j] - right_tr[j])
		- (right_sum[j] - right_tr_sum[j]) * (right_sum[j] - right_tr_sum[j]) 
		/ ((right_wt[j] - right_tr[j]) * (right_wt[j] - right_tr[j]));
	node_effect[j] = alpha * temp[j] * temp[j] * right_wt[j] - (1 - alpha) * (1 + train_to_est_ratio) 
		* right_wt[j] * (tr_var[j] / right_tr[j]  + con_var[j] / (right_wt[j] - right_tr[j]));
	}
	for (j = 0; j < ntreats; j++)
 {
	if (nclass == 0) {
		/* continuous predictor */
		/*
		left_wt = 0;
		left_tr = 0;
		left_n = 0;
		left_sum = 0;
		left_tr_sum = 0;
		left_sqr_sum = 0;
		left_tr_sqr_sum = 0;
		best = 0;
    */
	
		//run loop instead to set the above to zero: memset to 0
	  
	  memset(left_wt, 0, ntreats);
	  memset(left_tr, 0, ntreats);
	  memset(left_sum, 0, ntreats);
	  memset(left_tr_sum, 0, ntreats);
	  memset(left_sqr_sum, 0, ntreats);
	  memset(left_tr_sqr_sum, 0, ntreats);
	  memset(left_n, 0, ntreats);
	  memset(best, 0, ntreats);
	  	
		for (i = 0; right_n > edge; i++) {
			left_wt[j] += wt[i];
			right_wt[j] -= wt[i];
			left_tr[j] += wt[i] * (treatment[i]==j);
			right_tr[j] -= wt[i] * (treatment[i]==j);
			left_n[j]++;
			right_n[j]--;
			temp[j] = *y[i] * wt[i] * (treatment[i]==j);
			left_tr_sum[j] += temp;
			right_tr_sum[j] -= temp;
			left_sum[j] += *y[i] * wt[i];
			right_sum[j] -= *y[i] * wt[i];
			temp[j] = (*y[i]) *  (*y[i]) * wt[i];
			left_sqr_sum[j] += temp[j];
			right_sqr_sum[j] -= temp[j];
			temp[j] = (*y[i]) * (*y[i]) * wt[i] * (treatment[i]==j);
			left_tr_sqr_sum[j] += temp[j];
			right_tr_sqr_sum[j] -= temp[j];


			if (x[i + 1] != x[i] && left_n[j] >= edge &&
					(int) left_tr[j] >= min_node_size &&
					(int) left_wt[j] - (int) left_tr[j] >= min_node_size &&
					(int) right_tr[j] >= min_node_size &&
					(int) right_wt[j] - (int) right_tr[j] >= min_node_size) {

				left_temp[j] = left_tr_sum[j] / left_tr[j] - 
					(left_sum[j] - left_tr_sum[j]) / (left_wt[j] - left_tr[j]);
				left_tr_var[j] = left_tr_sqr_sum[j] / left_tr[j] - 
					left_tr_sum[j]  * left_tr_sum[j] / (left_tr[j] * left_tr[j]);
				left_con_var[j] = (left_sqr_sum[j] - left_tr_sqr_sum[j]) / (left_wt[j] - left_tr[j])  
					- (left_sum[j] - left_tr_sum[j]) * (left_sum[j] - left_tr_sum[j])
					/ ((left_wt[j] - left_tr[j]) * (left_wt[j] - left_tr[j]));        
				left_effect[j] = alpha * left_temp[j] * left_temp[j] * left_wt[j]
					- (1 - alpha) * (1 + train_to_est_ratio) * left_wt[j] 
					* (left_tr_var[j] / left_tr[j] + left_con_var[j] / (left_wt[j] - left_tr[j]));

				right_temp[j] = right_tr_sum[j] / right_tr[j] -
					(right_sum[j] - right_tr_sum[j]) / (right_wt[j] - right_tr[j]);
				right_tr_var[j] = right_tr_sqr_sum[j] / right_tr[j] -
					right_tr_sum[j] * right_tr_sum[j] / (right_tr[j] * right_tr[j]);
				right_con_var[j] = (right_sqr_sum[j] - right_tr_sqr_sum[j]) / (right_wt[j] - right_tr[j])
					- (right_sum[j] - right_tr_sum[j]) * (right_sum[j] - right_tr_sum[j]) 
					/ ((right_wt[j] - right_tr[j]) * (right_wt[j] - right_tr[j]));
				right_effect[j] = alpha * right_temp[j] * right_temp[j] * right_wt[j]
					- (1 - alpha) * (1 + train_to_est_ratio) * right_wt[j] * 
					(right_tr_var[j] / right_tr[j] + right_con_var[j] / (right_wt[j] - right_tr[j]));

				temp[j] = left_effect[j] + right_effect[j] - node_effect[j];
				if (temp[j] > best[j]) {
					best[j] = temp[j];
					where[j] = i;               
					if (left_temp[j] < right_temp[j])
						direction[j] = LEFT;
					else
						direction[j] = RIGHT;
				}             
			}
		}

		//*improve = best;
	  improve[j] = best;
		if (best[j] > 0) {         /* found something */
			csplit[0][j] = direction[j];
			//*split = (x[where] + x[where + 1]) / 2; 
			split[j] = (x[where[j]] + x[where[j] + 1]) / 2; 
		}
	}

	/*
	 * Categorical predictor
	 */
	else {
		for (i = 0; i < nclass; i++) {
			countn[i] = 0;
			wts[i] = 0;
			trs[i] = 0;
			sums[i] = 0;
			wtsums[i] = 0;
			trsums[i] = 0;
			wtsqrsums[i] = 0;
			trsqrsums[i] = 0;
		}

		/* rank the classes by treatment effect */
		for (i = 0; i < n; i++) {
			j = (int) x[i] - 1;
			countn[j]++;
			wts[j] += wt[i];
			trs[j] += wt[i] * treatment[i];
			sums[j] += *y[i];
			wtsums[j] += *y[i] * wt[i];
			trsums[j] += *y[i] * wt[i] * treatment[i];
			wtsqrsums[j] += (*y[i]) * (*y[i]) * wt[i];
			trsqrsums[j] +=  (*y[i]) * (*y[i]) * wt[i] * treatment[i];
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
		left_sqr_sum = 0.;
		left_tr_sqr_sum = 0.;

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

			left_tr_sqr_sum += trsqrsums[j];
			right_tr_sqr_sum -= trsqrsums[j];

			if (left_n >= edge && right_n >= edge &&
					(int) left_tr >= min_node_size &&
					(int) left_wt - (int) left_tr >= min_node_size &&
					(int) right_tr >= min_node_size &&
					(int) right_wt - (int) right_tr >= min_node_size) {

				left_temp = left_tr_sum / left_tr - (left_sum - left_tr_sum) 
					/ (left_wt - left_tr);

				left_tr_var = left_tr_sqr_sum / left_tr 
					- left_tr_sum  * left_tr_sum / (left_tr * left_tr);
				left_con_var = (left_sqr_sum - left_tr_sqr_sum) / (left_wt - left_tr)  
					- (left_sum - left_tr_sum) * (left_sum - left_tr_sum)
					/ ((left_wt - left_tr) * (left_wt - left_tr));       
				left_effect = alpha * left_temp * left_temp * left_wt
					- (1 - alpha) * (1 + train_to_est_ratio) * left_wt * 
					(left_tr_var / left_tr + left_con_var / (left_wt - left_tr));

				right_temp = right_tr_sum / right_tr - (right_sum - right_tr_sum) 
					/ (right_wt - right_tr);
				right_tr_var = right_tr_sqr_sum / right_tr 
					- right_tr_sum * right_tr_sum / (right_tr * right_tr);
				right_con_var = (right_sqr_sum - right_tr_sqr_sum) / (right_wt - right_tr)
					- (right_sum - right_tr_sum) * (right_sum - right_tr_sum) 
					/ ((right_wt - right_tr) * (right_wt - right_tr));
				right_effect = alpha * right_temp * right_temp * right_wt
					- (1 - alpha) * (1 + train_to_est_ratio) * right_wt *
					(right_tr_var / right_tr + right_con_var / (right_wt - right_tr));
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
}


double
policypred(double *y, double wt, double treatment, double *yhat, double propensity)
{
	double ystar;
	double temp;

	ystar = y[0] * (treatment - propensity) / (propensity * (1 - propensity));
	temp = ystar - *yhat;
	return temp * temp * wt;
}

