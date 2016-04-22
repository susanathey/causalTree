/*
 * The four routines for split.Rule = TOT
 */
#include <math.h>
#include "causalTree.h"
#include "causalTreeproto.h"

#ifndef min
#define min(a,b) (((a) < (b)) ? (a) : (b))
#endif

#ifndef max
#define max(a,b)  (((a) > (b)) ? (a) : (b))
#endif

// set temporarily:
//static int bucketnum = 5;
//static int bucketMax = 40;


static int *n_bucket;
static double *wts_bucket, *trs_bucket;
static double *tr_end_bucket, *con_end_bucket;
static double *wtsums_bucket;


/*
 * Warning: need to change to discrete version of TOT
 */

static double *mean, *sums;
static double *wts;
static int *countn;
static int *tsplit;

int
totDinit(int n, double *y[], int maxcat, char **error,
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
    //Rprintf("in totD init, train_to_est_ratio= %f\n", *train_to_est_ratio);
    return 0;
}


void
totDss(int n, double *y[], double *value, double *con_mean, double *tr_mean, double *risk, 
       double *wt, double *treatment, double max_y, double propensity) 
{
    int i;
    double temp = 0., twt = 0.;
    double mean, ss;
    double ystar;
    double temp0, temp1; // control; treatment.
    double trs, cons;
    
    temp0 = 0.;
    temp1 = 0.;
    trs = 0.;
    cons = 0.;

    
    for (i = 0; i < n; i++) {
        ystar = *y[i] * (treatment[i] - propensity) / (propensity * (1 - propensity));
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
    
    ss = 0.;
    for (i = 0; i < n; i++) {
        ystar = *y[i] * (treatment[i] - propensity) / (propensity * (1 - propensity));
        temp = ystar - mean;
        ss += temp * temp * wt[i];
    }
    
    *con_mean  = temp0 / cons;
    *tr_mean = temp1 / trs;
    *value = mean;
    *risk = ss;
}

/*
 * The anova splitting function.  Find that split point in x such that
 *  the sum of squares of y within the two groups is decreased as much
 *  as possible.  It is not necessary to actually calculate the SS, the
 *  improvement involves only means in the two groups.
 */

void totD(int n, double *y[], double *x, int nclass, int edge, double *improve, 
         double *split, int *csplit, double myrisk, double *wt, double *treatment, 
         double propensity, int minsize, int bucketnum, int bucketMax)
{
    int i, j;
    double temp;
    double left_sum, right_sum;
    double left_mean, right_mean;
    double left_wt, right_wt;
    int left_n, right_n;
    double left_tr, right_tr;
    double grandmean, best;
    int direction = LEFT;
    int where = 0;
    double ystar;
    int bucketTmp;
    double trsum = 0.;
    int Numbuckets;
    int min_node_size = minsize;
    
    
    double *cum_wt, *tmp_wt, *fake_x;
    double tr_wt_sum, con_wt_sum, con_cum_wt, tr_cum_wt;
    
    
    right_wt = 0;
    right_n = n;
    right_tr = 0;
    right_sum = 0;
    trsum = 0.;
    for (i = 0; i < n; i++) {
        ystar = *y[i] * (treatment[i] - propensity) / (propensity * (1 - propensity));
        right_sum += ystar * wt[i];
        right_wt += wt[i];
        right_tr += treatment[i] * wt[i];
        trsum += treatment[i];
    }
    grandmean = right_sum / right_wt;
    
    
    if(nclass == 0) {
        cum_wt = (double *) ALLOC(n, sizeof(double));
        tmp_wt = (double *) ALLOC(n, sizeof(double));
        fake_x = (double *) ALLOC(n, sizeof(double));
        tr_wt_sum = 0.;
        con_wt_sum = 0.;
        con_cum_wt = 0.;
        tr_cum_wt = 0.;
        
        for (i = 0; i < n; i ++) {
            if (treatment[i] == 0) {
                con_wt_sum += wt[i];
            } else {
                tr_wt_sum += wt[i];
            }
            cum_wt[i] = 0.;
            tmp_wt[i] = 0.;
            fake_x[i] = 0.;
        }
        
        //bucketTmp = min(round(trsum / (double)bucketnum), round(((double)n - trsum) / (double)bucketnum));
        int test1 = round(trsum / (double)bucketnum);
        int test2 = round(((double)n - trsum) / (double)bucketnum);
        bucketTmp = min(test1, test2);
        Numbuckets = max(minsize, min(bucketTmp, bucketMax));
        
        
        n_bucket = (int *) ALLOC(Numbuckets,  sizeof(int));
        wts_bucket = (double *) ALLOC(Numbuckets, sizeof(double));
        trs_bucket = (double *) ALLOC(Numbuckets, sizeof(double));
        tr_end_bucket = (double *) ALLOC(Numbuckets, sizeof(double));
        con_end_bucket = (double *) ALLOC (Numbuckets, sizeof(double));
        wtsums_bucket = (double *) ALLOC (Numbuckets, sizeof(double));
        
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
        
        for (j = 0; j < Numbuckets; j++) {
            n_bucket[j] = 0;
            wts_bucket[j] = 0.;
            trs_bucket[j] = 0.;
            wtsums_bucket[j]  = 0.;
        }
        
        for (i = 0; i < n; i++) {
            j = fake_x[i];
            n_bucket[j]++;
            wts_bucket[j] += wt[i];
            trs_bucket[j] += wt[i] * treatment[i];
            ystar = *y[i] * (treatment[i] - propensity) / (propensity * (1 - propensity));
            wtsums_bucket[j] += ystar * wt[i];
            if (treatment[i] == 1) {
                tr_end_bucket[j] = x[i];
            } else {
                con_end_bucket[j] = x[i];
            }
        }
        
        left_sum = 0;
        left_wt = 0;
        left_n = 0;
        left_tr = 0.;
        best = 0;
        
        for (j = 0; j < Numbuckets; j++) {
            
            left_n += n_bucket[j];
            right_n -= n_bucket[j];
            left_wt += wts_bucket[j];
            right_wt -= wts_bucket[j];
            left_tr += trs_bucket[j]; // weighted
            right_tr -= trs_bucket[j];
            
            left_sum += wtsums_bucket[j];
            right_sum -= wtsums_bucket[j];
            
            
            if (left_n >= edge && right_n >= edge &&
                (int) left_tr >= min_node_size &&
                (int) left_wt - (int) left_tr >= min_node_size &&
                (int) right_tr >= min_node_size &&
                (int) right_wt - (int) right_tr >= min_node_size) {
                
                left_mean = left_sum / left_wt;
                right_mean = right_sum / right_wt;
                temp = left_wt * (grandmean - left_mean) * (grandmean - left_mean) + 
                    right_wt * (grandmean - right_mean) * (grandmean - right_mean);  
                
                //Rprintf("temp = %f\n", temp);
                
                if (temp > best) {
                    best = temp;
                    where = j;
                    if (left_sum < right_sum)
                        direction = LEFT;
                    else
                        direction = RIGHT;
                }
                
            }
        }
        
        *improve = best / myrisk;
        if (best > 0) {
            csplit[0] = direction;
            *split = (tr_end_bucket[where] + con_end_bucket[where]) / 2;
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
    totDpred(double *y, double wt, double treatment, double *yhat, double propensity) // pass in ct.which
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

