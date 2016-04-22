/*
 * The four routines for split.Rule = TOT
 */
#include <math.h>
#include "causalTree.h"
#include "causalTreeproto.h"

/*
 * Warning: need to change to discrete version of TOT
 */
//static double *sums, *wtsums, *treatment_effect;
//static double *wts, *trs, *trsums;
//static int *countn;
//static int *tsplit;
static double *mean, *sums;
static double *wts;
static int *countn;
static int *tsplit;

int
totinit(int n, double *y[], int maxcat, char **error,
        int *size, int who, double *wt, double *treatment, 
        int bucketnum, int bucketMax, double *train_to_est_ratio)
{
    if (who == 1 && maxcat > 0) {
        graycode_init0(maxcat);
        countn = (int *) ALLOC(2 * maxcat, sizeof(int));
        tsplit = countn + maxcat;
        mean = (double *) ALLOC(3 * maxcat, sizeof(double));
        wts = mean + maxcat;
        sums = wts + maxcat;

    }
    *size = 1;
    *train_to_est_ratio = n * 1.0 / ct.NumHonest;
    return 0;
}


void
totss(int n, double *y[], double *value,  double *con_mean, double *tr_mean, double *risk, double *wt, 
      double *treatment, double max_y, double propensity)
{
    int i;
    double temp = 0., twt = 0.;
    double temp0, temp1; // control; treatment.
    double trs, cons;
    double mean, ss;
    //double temp0 = 0., temp1 = 0., twt = 0.; /* sum of the weights */ 
    double ystar;
    temp0 = 0.;
    temp1 = 0.;
    trs = 0.;
    cons = 0.;

    //Rprintf("propensity score = %f", propensity);
    
    for (i = 0; i < n; i++) {
        //Rprintf("treatment[%d] = %f\n", i, treatment[i]);
        ystar = *y[i] * (treatment[i] - propensity) / (propensity * (1 - propensity));
        //Rprintf("y[%d] = %f, ystar[%d] = %f\n", i, *y[i], i, ystar);
        temp += ystar * wt[i];
        twt += wt[i];
        if (treatment[i] == 0) {
            //con
            temp0 += *y[i] * wt[i];
            cons += wt[i];
        } else {
            //tr
            temp1 += *y[i] * wt[i];
            trs += wt[i];
        }
    }
    mean = temp / twt;
    //Rprintf("mean = %f\n", mean);
    
    ss = 0.;
    for (i = 0; i < n; i++) {
        ystar = *y[i] * (treatment[i] - propensity) / (propensity * (1 - propensity));
        temp = ystar - mean;
        ss += temp * temp * wt[i];
    }
    *con_mean  = temp0 / cons;
    *tr_mean = temp1 / trs;
    *value = temp1 /trs - temp0 / cons;
    //*value = mean;
    *risk = ss;
}

/*
 * The anova splitting function.  Find that split point in x such that
 *  the sum of squares of y within the two groups is decreased as much
 *  as possible.  It is not necessary to actually calculate the SS, the
 *  improvement involves only means in the two groups.
 */

void tot(int n, double *y[], double *x, int nclass, int edge, double *improve, 
         double *split, int *csplit, double myrisk, double *wt, double *treatment,
         double propensity, int minsize)
{
    int i, j;
    double temp;
    double left_sum, right_sum;
    double left_wt, right_wt;
    int left_n, right_n;
    double grandmean, best;
    int direction = LEFT;
    int where = 0;
    double ystar;
    double left_mean, right_mean;
    double right_tr, left_tr;
    double min_node_size = minsize;
    
    right_wt = 0;
    right_n = n;
    right_sum = 0;
    right_tr = 0;
    for (i = 0; i < n; i++) {
        ystar = *y[i] * (treatment[i] - propensity) / (propensity * (1 - propensity));
        right_sum += ystar * wt[i];
        right_wt += wt[i];
        right_tr += wt[i] * treatment[i];

    }
    grandmean = right_sum / right_wt;

    
    if(nclass == 0) {
        left_sum = 0;
        left_wt = 0;
        left_n = 0;
        left_tr = 0;
        best = 0;
        for (i = 0; right_n > edge; i++) {
            left_wt += wt[i];
            right_wt -= wt[i];
            left_tr += wt[i] * treatment[i];
            right_tr -= wt[i] * treatment[i];
            left_n++;
            right_n--;
            ystar = *y[i] * (treatment[i] - propensity) / (propensity * (1 - propensity));
            left_sum += ystar * wt[i];
            right_sum -= ystar * wt[i];

            if (x[i + 1] != x[i] && left_n >= edge &&
                (int) left_tr >= min_node_size &&
                (int) left_wt - (int) left_tr >= min_node_size &&
                (int) right_tr >= min_node_size &&
                (int) right_wt - (int) right_tr >= min_node_size) {
    
                //if (x[i + 1] != x[i] && left_n >= edge) {
                left_mean = left_sum / left_wt;
                right_mean = right_sum / right_wt;
                temp = left_wt * (grandmean - left_mean) * (grandmean - left_mean) + 
                    right_wt * (grandmean - right_mean) * (grandmean - right_mean);
        
                if (temp > best) {
                    best = temp;
                    where = i;
                    if (left_mean < right_mean)
                        direction = LEFT;
                    else
                        direction = RIGHT;
                }
            }
        }
        
        *improve = best / myrisk;
        if (best > 0) {
            csplit[0] = direction;
            *split = (x[where] + x[where + 1]) / 2;
        }
    } else {
        /*
         * Categorical Predictor
         */
        for (i = 0; i < nclass; i++) {
            sums[i] = 0;
            countn[i] = 0;
            wts[i] = 0;
        }
        /* rank the classes by their mean y value */
        for (i = 0; i < n; i++) {
            j = (int) x[i] - 1;
            countn[j]++;
            wts[j] += wt[i];
            ystar = *y[i] * (treatment[i] - propensity) / (propensity * (1 - propensity));
            sums[j] += (ystar - grandmean) * wt[i];
        }
        for (i = 0; i < nclass; i++) {
            if (countn[i] > 0) {
                tsplit[i] = RIGHT;
                mean[i] = sums[i] / wts[i];
            } else {
                tsplit[i] = 0;
            }
        }
        graycode_init2(nclass, countn, mean);
        
        /*
         * Now find the split that we want
         */
        left_wt = 0;
        left_sum = 0;
        right_sum = 0;
        left_n = 0;
        best = 0;
        where = 0;
        while ((j = graycode()) < nclass) {
            tsplit[j] = LEFT;
            left_n += countn[j];
            right_n -= countn[j];
            left_wt += wts[j];
            right_wt -= wts[j];
            left_sum += sums[j];
            right_sum -= sums[j];
            if (left_n >= edge && right_n >= edge) {
                temp = left_sum * left_sum / left_wt +
                    right_sum * right_sum / right_wt;
                if (temp > best) {
                    best = temp;
                    if ((left_sum / left_wt) > (right_sum / right_wt)) {
                        for (i = 0; i < nclass; i++) csplit[i] = -tsplit[i];
                    } else {
                        for (i = 0; i < nclass; i++) csplit[i] = tsplit[i];
                    }
                }
            }
        }
        
        *improve = best / myrisk;  /* improvement */
    }
}

double
    totpred(double *y, double wt, double treatment, double *yhat, double propensity) // pass in ct.which
    {
        double ystar;
        double temp;
        
        ystar = y[0] * (treatment - propensity) / (propensity * (1 - propensity));
        temp = ystar - *yhat;
        //if (treatment == 1)  temp = y[0] / propensity;
        //else temp = - y[0] / (1 - propensity);
        //double temp = y[0] - *yhat;
        //temp -= *yhat;
        return temp * temp * wt;
    }

double totxeval(int *unique_leaf, int **val_leaf_mat, int cp_id, int t, int *sorts, double *wt,
                double *treatment,  double *y[], double propensity, int k, int nobs, double val_sum_wt, int val_count) {
    /*
    * nobs = ct.n
    * t = # of unique leaves
    * k: start pt of validation set
    * sorts = ct.sorts[0][]
    * wt = ct.wt
    * cp_id = jj
    * treatment = ct.treatment
    * y[] = ct.ydata
    * val_sum_wt = total wts in the validation set:
    */
    
    int s, i, j;
    int jj = cp_id;
    double leaf_total_error;
    double leaf_local_wt;
    double leaf_local_wt_sum;
    double leaf_local_wt_sq_sum;
    double ystar, ystarMean;
    int node;
    
    leaf_total_error = 0.;
    for (s = 0; s < t; s++) {
        // in leaf s:
        leaf_local_wt = 0.;
        leaf_local_wt_sum = 0.;
        leaf_local_wt_sq_sum = 0.;
        node = unique_leaf[s];
        if (node == 0) {
            Rprintf("There is a bug!\n");
        }
        for (i = k; i < nobs; i++) {
            // consider the obs in the node:
            j = sorts[i];
            if(val_leaf_mat[j][jj] == node) {
                ystar = y[0][j] * (treatment[j] - propensity) / (propensity * (1 - propensity));
                leaf_local_wt += wt[j];
                leaf_local_wt_sum += wt[j] * ystar;
                leaf_local_wt_sq_sum += wt[j] * ystar * ystar;
            }
        }
        ystarMean = leaf_local_wt_sum / leaf_local_wt;
        leaf_total_error += leaf_local_wt_sq_sum - leaf_local_wt * ystarMean;
    }
    if (leaf_total_error != leaf_total_error) Rprintf("Nan value happens here!");
    return leaf_total_error;
}



