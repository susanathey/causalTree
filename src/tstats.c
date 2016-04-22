/*
 * The fifth routines for t - statistic splitting:
 */
#include "causalTree.h"
#include "causalTreeproto.h"

// these only for categorical variables:
static double *mean, *sums, *wtsums;
static double *wts, *trs, *trsums;
static double *wtsqrsums, *wttrsqrsums;
static int *countn;
static int *tsplit;

int
tstatsinit(int n, double *y[], int maxcat, char **error,
        int *size, int who, double *wt, double *treatment, int bucketnum, 
        int bucektMax, double *train_to_est_ratio)
{
    //Rprintf("maxcat = %d\n", maxcat);
    if (who == 1 && maxcat > 0) {
        graycode_init0(maxcat);
        countn = (int *) ALLOC(2 * maxcat, sizeof(int));
        tsplit = countn + maxcat;
        mean = (double *) ALLOC(6 * maxcat, sizeof(double));
        wts = mean + maxcat;
        trs = wts + maxcat;
        sums = trs + maxcat;
        wtsums = sums + maxcat;
        trsums = wtsums + maxcat;
    }
    *size = 1;
    *train_to_est_ratio = n * 1.0 / ct.NumHonest;
    return 0;
}

/*
 * The anova evaluation function.  Return the mean and the ss.
 */
void
tstatsss(int n, double *y[], double *value, double *con_mean,  double *tr_mean, double *risk,
         double *wt, double *treatment, double max_y, double alpha, double train_to_est_ratio)
{
    int i;
    double temp = 0., temp0 = 0., temp1 = 0., twt = 0.; /* sum of the weights */ 
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
        con_sqr_sum += (*y[i]) * (*y[i]) * wt[i] * (1 - treatment[i]);
    }

    effect = temp1 / ttreat - temp0 / (twt - ttreat);
    //Rprintf("effect = %f\n", effect);
    tr_var = tr_sqr_sum / ttreat - temp1 * temp1 / (ttreat * ttreat);
    con_var = con_sqr_sum / (twt - ttreat) - temp0 * temp0 
        / ((twt - ttreat) * (twt - ttreat));
    //Rprintf("tr_var = %f, con_var = %f\n", tr_var, con_var);
    
    *tr_mean = temp1 / ttreat;
    *con_mean = temp0 / (twt - ttreat);
    *value = effect;
    *risk = 4 * n * max_y * max_y - alpha * n * effect * effect 
    + (1 - alpha) * (1 + train_to_est_ratio) * n 
    * (tr_var /ttreat  + con_var / (twt - ttreat));
}

/*
 * The anova splitting function.  Find that split point in x such that
 *  the sum of squares of y within the two groups is decreased as much
 *  as possible.  It is not necessary to actually calculate the SS, the
 *  improvement involves only means in the two groups.
 */

void
tstats(int n, double *y[], double *x, int nclass,
        int edge, double *improve, double *split, int *csplit,
        double myrisk, double *wt, double *treatment, int minsize, double alpha, 
        double train_to_est_ratio)
{
    int i, j;
    double temp;
    double left_sum, right_sum;
    double left_tr_sum, right_tr_sum;
    // add squared sum:
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
    double left_var, right_var, sd;
    double left_temp, right_temp;
    int min_node_size = minsize;
    double improve_temp, improve_best;
    
    //Rprintf("in tstats, train_to_est_ratio = %f\n", train_to_est_ratio);

    right_wt = 0.;
    right_tr = 0.;
    right_sum = 0.;
    right_tr_sum = 0.;
    right_sqr_sum = 0.;
    right_tr_sqr_sum = 0.;
    right_n = n;
    improve_temp = 0.;
    improve_best = 0.;
    for (i = 0; i < n; i++) {
        //Rprintf("wt[%d] = %f\n", i, wt[i]);
        right_wt += wt[i];
        //Rprintf("right_wt = %f\n", right_wt);
        right_tr += wt[i] * treatment[i];
        right_sum += *y[i] * wt[i];
        right_tr_sum += *y[i] * wt[i] * treatment[i];
        right_sqr_sum += (*y[i]) * (*y[i]) * wt[i];
        right_tr_sqr_sum += (*y[i]) * (*y[i]) * wt[i] * treatment[i];
    }
    //Rprintf("right_wt = %f\n", right_wt);
    temp = right_tr_sum / right_tr - (right_sum - right_tr_sum) / (right_wt - right_tr);
    tr_var = right_tr_sqr_sum / right_tr - right_tr_sum * right_tr_sum 
        / (right_tr * right_tr);
    con_var = (right_sqr_sum - right_tr_sqr_sum) / (right_wt - right_tr)
        - (right_sum - right_tr_sum) * (right_sum - right_tr_sum) 
        / ((right_wt - right_tr) * (right_wt - right_tr));
    //now we use t-statistic
    node_effect = alpha * temp * temp * n - (1 - alpha) * (1 + train_to_est_ratio) 
        * n * (tr_var / right_tr  + con_var / (right_wt - right_tr));
    // temporarily set as 0:
    //node_effect = 0
   //Rprintf("n = %d, node_effect = %f\n", n, node_effect);
   
    if (nclass == 0) {
        /* continuous predictor */
        left_wt = 0.;
        left_tr = 0.;
        left_n = 0;
        left_sum = 0.;
        left_tr_sum = 0.;
        left_sqr_sum = 0.;
        left_tr_sqr_sum = 0.;

        best = 0.;

        for (i = 0; right_n > edge; i++) {
            //Rprintf("wt[%d] = %f", i, wt[i]);
            //Rprintf("treatment[%d] = %f", i, treatment[i]);
            left_wt += wt[i];
            right_wt -= wt[i];

            left_tr += wt[i] * treatment[i];
            right_tr -= wt[i] * treatment[i];
            //Rprintf("test = %f\n", right_tr);
            //Rprintf("left_tr = %f\n", left_tr);

            left_n++;
            right_n--;

            temp = *y[i] * wt[i] * treatment[i];
            left_tr_sum += temp;
            right_tr_sum -= temp;

            left_sum += *y[i] * wt[i];
            right_sum -= *y[i] * wt[i];

            temp = (*y[i]) * (*y[i]) * wt[i] * treatment[i];
            left_tr_sqr_sum += temp;
            right_tr_sqr_sum -= temp;

            temp = (*y[i]) * (*y[i]) * wt[i];
            left_sqr_sum += temp;
            right_sqr_sum -= temp;

            //Rprintf("left_n = %d, edge = %d\n", left_n, edge);
            //Rprintf("left_tr = %d, left_wt = %d, right_tr = %d, right_wt = %d\n", left_tr, left_wt, right_tr, right_wt);
            //Rprintf("minsize = %d\n", min_node_size);
            // I have a question about weights here in report:
            if (x[i + 1] != x[i] && left_n >= edge &&
                    (int) left_tr >= min_node_size &&
                    (int) left_wt - (int) left_tr >= min_node_size &&
                    (int) right_tr >= min_node_size &&
                    (int) right_wt - (int) right_tr >= min_node_size) {
                //Rprintf("get in!\n");
                //Rprintf("left_tr = %f\n", left_tr);
                //Rprintf("left_con = %f\n", left_wt - left_tr);
                //Rprintf("right_tr = %f\n", right_tr);
                //Rprintf("right_con = %f\n", right_wt - right_tr);
                //Rprintf("left_n = %d\n", left_n);
                //Rprintf("right_n = %d\n", right_n);

                left_temp = left_tr_sum / left_tr - (left_sum - left_tr_sum) / (left_wt - left_tr);
                left_tr_var = left_tr_sqr_sum / left_tr - left_tr_sum  * left_tr_sum 
                    / (left_tr * left_tr);
                left_con_var = (left_sqr_sum - left_tr_sqr_sum) / (left_wt - left_tr) 
                    - (left_sum - left_tr_sum) * (left_sum - left_tr_sum)
                    / ((left_wt - left_tr) * (left_wt - left_tr));        

                left_var = left_tr_var / left_tr + left_con_var / (left_wt - left_tr);
                left_effect = alpha * left_temp * left_temp * left_n
                  - (1 - alpha) * (1 + train_to_est_ratio) * left_n 
                    * (left_tr_var / left_tr + left_con_var / (left_wt - left_tr));

                //taumean = right_tr / right_wt;
                //temp = (right_tr_sum - right_sum * taumean) /
                //((1 - taumean) * right_tr);
                right_temp = right_tr_sum / right_tr - (right_sum - right_tr_sum) / (right_wt - right_tr);
                right_tr_var = right_tr_sqr_sum / right_tr - right_tr_sum * right_tr_sum / (right_tr * right_tr);
                right_con_var = (right_sqr_sum - right_tr_sqr_sum) / (right_wt - right_tr)
                    - (right_sum - right_tr_sum) * (right_sum - right_tr_sum) / ((right_wt - right_tr) * (right_wt - right_tr));
                right_var = right_tr_var / right_tr + right_con_var / (right_wt - right_tr);
                right_effect = alpha * right_temp * right_temp * right_n
                    - (1 - alpha) * (1 + train_to_est_ratio) * right_n 
                    * (right_tr_var / right_tr + right_con_var / (right_wt - right_tr));    

                sd = sqrt(left_var / left_wt  + right_var / right_wt);
                // very strange !!!!!!
                temp = fabs(left_temp - right_temp) / sd;
                
                //Rprintf("temp = %f\n", temp);

                improve_temp = left_effect + right_effect - node_effect;
                //Rprintf("improve_temp = %f\n", improve_temp);
                //Rprintf("at %f,leftn: %d, lefteffect: %f, rightn: %d, righteffect: %f\n", x[i], left_n, left_effect,right_n, right_effect, node_effect);
                //Rprintf("current best is %f, and temp improv = %f.\n", best, improve_temp);
                //Rprintf("best  = %f\n", best);
                if (temp > best) {
                  
                  //Rprintf("best  = %f\n", best);
                    best = temp;
                    where = i;
                    improve_best = improve_temp;
                    if (left_temp < right_temp)
                        direction = LEFT;
                    else
                        direction = RIGHT;
                }             
            }
        }
       
        // change the value of improvement:
        *improve = improve_best;
        //Rprintf("final improve = %f\n", improve_best);
               
        
        if (improve_best > 0) {         /* found something */
            //Rprintf("best = %f, split = %f\n", improve_best, (x[where] + x[where + 1]) / 2 );
            csplit[0] = direction;
            *split = (x[where] + x[where + 1]) / 2; /* where to split!!!!!!!!! */ 
            //Rprintf("split = %f\n", *split);
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
            wttrsqrsums[i] = 0;
        }

        /* rank the classes by their mean y value */
        /* RANK THE CLASSES BY THEI */
        // Rprintf("nclass = %d ", nclass);
        for (i = 0; i < n; i++) {
            j = (int) x[i] - 1;
            // Rprintf("%d cat, ", j);
            countn[j]++;
            wts[j] += wt[i];
            trs[j] += wt[i] * treatment[i];
            sums[j] += *y[i];
            // adding part
            wtsums[j] += *y[i] * wt[i];
            trsums[j] += *y[i] * wt[i] * treatment[i];
            wtsqrsums[j] += (*y[i]) * (*y[i]) * wt[i];
            wttrsqrsums[j] += (*y[i]) * (*y[i]) * wt[i] * treatment[i];
        }

        for (i = 0; i < nclass; i++) {
            if (countn[i] > 0) {
                tsplit[i] = RIGHT;
                mean[i] = sums[i] / wts[i];
                // mean[i] = sums[i] / countn[i];
                //Rprintf("countn[%d] = %d, mean[%d] = %f\n", i, countn[i], i, mean[i]);
            } else
                tsplit[i] = 0;
        }
        graycode_init2(nclass, countn, mean);

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
            //Rprintf("graycode()= %d\n", j);
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
                    - (left_sum - left_tr_sum) * (left_sum - left_tr_sum)/ ((left_wt - left_tr) * (left_wt - left_tr));   
                left_var = left_tr_var / left_tr + left_con_var / (left_wt - left_tr);
                left_effect = alpha * left_temp * left_temp * left_wt
                  - (1 - alpha) * (1 + train_to_est_ratio) * left_wt
                    * (left_tr_var / left_tr + left_con_var / (left_wt - left_tr));

                //Rprintf("left_sum = %f, left_wt_sum = %f, left_wt = %f, left_n = %d\n", left_sum, left_wt_sum, left_wt, left_n);             
                right_temp = right_tr_sum / right_tr - (right_sum - right_tr_sum) / (right_wt - right_tr);
                right_tr_var = right_tr_sqr_sum / right_tr - right_tr_sum * right_tr_sum / (right_tr * right_tr);
                right_con_var = (right_sqr_sum - right_tr_sqr_sum) / (right_wt - right_tr)
                    - (right_sum - right_tr_sum) * (right_sum - right_tr_sum) / ((right_wt - right_tr) * (right_wt - right_tr));
                right_var = right_tr_var / right_tr + right_con_var / (right_wt - right_tr);
                right_effect = alpha * right_temp * right_temp * right_wt
                - (1 - alpha) * (1 + train_to_est_ratio) 
                    * right_wt * (right_tr_var / right_tr + right_con_var / (right_wt - right_tr));

                sd = sqrt(left_var / left_wt  + right_var / right_wt);
                temp = fabs(left_temp - right_temp) / sd; 
                improve_temp = left_effect + right_effect - node_effect;
                //Rprintf("left_n= %d, lefteffect = %f, right_n = %d, righteffect = %f\n", left_n, left_effect, right_n, right_effect);               
                //temp = left_sum * left_sum / left_wt +
                //right_sum * right_sum / right_wt;
                if (temp > best) {
                    best = temp;
                    improve_best = improve_temp;
                }
            }
        }
        *improve = best;
        if (improve_best > 0) {
            if (left_temp > right_temp)
                for (i = 0; i < nclass; i++) csplit[i] = -tsplit[i];
            else
                for (i = 0; i < nclass; i++) csplit[i] = tsplit[i];
        }
    }
    //Rprintf("for %f variable, improv = %f\n", x[0], *improve);
}

double tstatspred(double *y, double wt, double treatment, double *yhat, double propensity) // pass in ct.which
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
