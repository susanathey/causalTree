/*
 * split.Rule = TOT
 */
#include <math.h>
#include "causalTree.h"
#include "causalTreeproto.h"

/*
 * Warning: need to change to discrete version of TOT
 */

static double *sums, *wtsums, *treatment_effect;
static double *wtsqrsums, *wttrsqrsums;
static double *wts, *trs, *trsums;
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
    return 0;
}


void
totss(int n, double *y[], double *value,  double *con_mean, double *tr_mean, double *risk, double *wt, 
      double *treatment, double max_y, double propensity)
{
    int i;
    double temp = 0., twt = 0.;
    double temp0, temp1; 
    double trs, cons;
    double mean, ss;
    double ystar;
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
    *value = temp1 /trs - temp0 / cons;
    *risk = ss;
}

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
      Rprintf("tot: inside cont. split\n");
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
        Rprintf("tot: inside factor split!\n");
        Rprintf("nclass:%d\n",nclass);
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
        left_sum = 0;
        right_sum = 0;
        left_n = 0;
        best = 0;
        where = 0;
        while ((j = graycode()) < nclass) {
          Rprintf("j=%d\n",j);
            tsplit[j] = LEFT;
            left_n += countn[j];
            right_n -= countn[j];
            left_wt += wts[j];
            right_wt -= wts[j];
            left_sum += sums[j];
            right_sum -= sums[j];
            Rprintf("j=%d,sums[j]=%f\n",j,sums[j]);
            Rprintf("left_sum=%f,right_sum=%f\n",left_sum,right_sum);
            if (left_n >= edge && right_n >= edge) {
              Rprintf("tot factor: inside >=edge if \n");
                temp = left_sum * left_sum / left_wt +
                    right_sum * right_sum / right_wt;
                Rprintf("temp=%f\n",temp);
                Rprintf("best=%f\n",best);
                Rprintf("left_sum_fin=%f,left_wt=%f,right_sum_fin=%f,right_wt=%f\n",left_sum,left_wt,right_sum,right_wt);
                if (temp > best) {
                    best = temp;
                  Rprintf("tot factor best:%f\n",best);
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
        return temp * temp * wt;
    }
