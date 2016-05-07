/*
 * The anova funciton inherited from rpart package. 
 */
#include "causalTree.h"
#include "causalTreeproto.h"

static double *sums, *wtsums, *treatment_effect;
static double *wts, *trs, *trsums;
static int *countn;
static int *tsplit;


int
anovainit(int n, double *y[], int maxcat, char **error,
        double *parm, int *size, int who, double *wt, double *treatment)
{
    if (who == 1 && maxcat > 0) {
        graycode_init0(maxcat);
        countn = (int *) ALLOC(2 * maxcat, sizeof(int));
        tsplit = countn + maxcat;
        treatment_effect = (double *) ALLOC(6 * maxcat, sizeof(double));
        wts = treatment_effect + maxcat;
        trs = wts + maxcat;
        sums = trs + maxcat;
        wtsums = sums + maxcat;
        trsums = wtsums + maxcat;

    }
    *size = 1;
    return 0;
}


void
anovass(int n, double *y[], double *value, double *risk, double *wt, double *treatment, double max_y)
{
    int i;
    double temp0 = 0., temp1 = 0., twt = 0.; /* sum of the weights */ 
    double ttreat = 0.;
    double effect;

    for (i = 0; i < n; i++) {
        temp1 += *y[i] * wt[i] * treatment[i];
        temp0 += *y[i] * wt[i] * (1 - treatment[i]);
        twt += wt[i];
        ttreat += wt[i] * treatment[i];
    }

    effect = temp1 / ttreat - temp0 / (twt - ttreat);
    *value = effect;
    *risk = 4 * twt * max_y * max_y - twt * effect * effect ;
}

/*
 * The anova splitting function.  Find that split point in x such that
 *  the sum of squares of y within the two groups is decreased as much
 *  as possible.  It is not necessary to actually calculate the SS, the
 *  improvement involves only means in the two groups.
 */

void
anova(int n, double *y[], double *x, int nclass,
        int edge, double *improve, double *split, int *csplit,
        double myrisk, double *wt, double *treatment, int minsize)
{
    int i, j;
    double temp;
    double left_sum, right_sum;
    double left_tr_sum, right_tr_sum;
    double left_tr, right_tr;
    double left_wt, right_wt;
    int left_n, right_n;
    double best;
    int direction = LEFT;
    int where = 0;
    double node_effect, left_effect, right_effect;
    double left_temp, right_temp;
    int min_node_size = minsize;

    right_wt = 0;
    right_tr = 0;
    right_sum = 0;
    right_tr_sum = 0;
    right_n = n;
    for (i = 0; i < n; i++) {
        right_wt += wt[i];
        right_tr += wt[i] * treatment[i];
        right_sum += *y[i] * wt[i];
        right_tr_sum += *y[i] * wt[i] * treatment[i];
    }

    temp = right_tr_sum / right_tr - (right_sum - right_tr_sum) / (right_wt - right_tr);
    node_effect = temp * temp * right_wt;
    
    if (nclass == 0) {
        /* continuous predictor */
        left_wt = 0;
        left_tr = 0;
        left_n = 0;
        left_sum = 0;
        left_tr_sum = 0;
        best = 0;
        for (i = 0; right_n > edge; i++) {
            left_wt += wt[i];
            right_wt -= wt[i];
            left_tr += wt[i] * treatment[i];
            right_tr -= wt[i] * treatment[i];
            left_n++;
            right_n--;
            temp = *y[i] * wt[i] * treatment[i];
            left_tr_sum += temp;
            right_tr_sum -= temp;
            left_sum += *y[i] * wt[i];
            right_sum -= *y[i] * wt[i];

            if (x[i + 1] != x[i] && left_n >= edge &&
                    (int) left_tr >= min_node_size &&
                    (int) left_wt - (int) left_tr >= min_node_size &&
                    (int) right_tr >= min_node_size &&
                    (int) right_wt - (int) right_tr >= min_node_size) {

                left_temp = left_tr_sum / left_tr - (left_sum - left_tr_sum) / (left_wt - left_tr);
                left_effect = left_temp * left_temp * left_wt;

                right_temp = right_tr_sum / right_tr - (right_sum - right_tr_sum) / (right_wt - right_tr);
                right_effect = right_temp * right_temp * right_wt;
                temp = left_effect + right_effect - node_effect;

                if (temp > best) {
                    best = temp;
                    where = i;               
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
            *split = (x[where] + x[where + 1]) / 2; /* where to split!!!!!!!!! */ 
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

            if (left_n >= edge && right_n >= edge &&
                    (int) left_tr >= min_node_size &&
                    (int) left_wt - (int) left_tr >= min_node_size &&
                    (int) right_tr >= min_node_size &&
                    (int) right_wt - (int) right_tr >= min_node_size) {

                left_temp = left_tr_sum / left_tr - (left_sum - left_tr_sum) / (left_wt - left_tr);
                left_effect = left_temp * left_temp * left_wt;

                right_temp = right_tr_sum / right_tr - (right_sum - right_tr_sum) / (right_wt - right_tr);
                right_effect = right_temp * right_temp * right_wt;
    
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
