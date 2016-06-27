/*
 * The cross validation evaluation function
 */
#include<stdio.h>
#include <math.h>
#include "causalTree.h"
#include "causalTreeproto.h"

double
tot_xpred(double *y, double wt, double treatment, double *yhat, double propensity) 
{
    double ystar;
    double temp;
        
    ystar = y[0] * (treatment - propensity) / (propensity * (1 - propensity));
    temp = ystar - yhat[0];
    return temp * temp;
}

double matching_xpred(double *y1, double *y2, double wt1, double wt2, double treatment1,
                     double treatment2, double pred1, double pred2) {
    double temp;
    temp = (2 * treatment1 - 1) * (*y1 - *y2);
    return (temp - 0.5 * (pred1 + pred2)) * (temp - 0.5 * (pred1 + pred2));
}


double fitH_xpred(double *y, double wt, double treatment, double tr_mean, 
                  double con_mean, double trs, double cons, double alpha,
                  double xtrain_to_est_ratio) {
    double res;
    double tr_var;
    double con_var;
    double tmp;
    double tmp_val;

    
    if (treatment == 0) {
        // con
        con_var = wt * (y[0] - con_mean) *  (y[0] - con_mean);
        tmp = con_var / cons;
        tmp_val = con_mean;
    } else {
        // tr
        tr_var = wt * (y[0] - tr_mean) * (y[0] - tr_mean);
        tmp = tr_var / trs;
        tmp_val = tr_mean;
    } 
    res =  4 * ct.max_y * ct.max_y - alpha *  tmp_val * tmp_val +
        (1 + xtrain_to_est_ratio / (ct.NumXval - 1)) * (1 - alpha) *  tmp; 
    return res;
}

double fitA_xpred(double *y, double wt, double treatment, double tree_tr_mean, double tree_con_mean) {
    double res;
    if (treatment == 0) {
        res = (y[0] - tree_con_mean) * (y[0] - tree_con_mean);
    } else {
        res =  (y[0] - tree_tr_mean) * (y[0] - tree_tr_mean);
    }
    return res;
}

double CTH_xpred(double *y, double wt, double treatment, double tr_mean,
                 double con_mean, double trs, double cons, double alpha, 
                 double xtrain_to_est_ratio, double propensity) {
   double res;
   double tr_var;
   double con_var;
   double tmp;
   if (treatment == 0) {
       // con
       con_var = wt * (y[0] - con_mean) *  (y[0] - con_mean);
       tmp = con_var / ((1 - propensity) * cons);
   } else {
       // tr
       tr_var = wt * (y[0] - tr_mean) * (y[0] - tr_mean);
       tmp = tr_var / (propensity * trs);
   } 
   double effect = tr_mean - con_mean;
   
   res = 4 * ct.max_y * ct.max_y - alpha *  effect * effect + (1 + xtrain_to_est_ratio / (ct.NumXval - 1)) 
       * (1 - alpha) *  tmp; 
   
   return res;
}

double CTA_xpred(double *y, double wt, double treatment, double tr_mean, double con_mean,
                     double tree_tr_mean, double tree_con_mean, double alpha) {
    double res;
    double effect_tr = tree_tr_mean - tree_con_mean;
    double effect_te = tr_mean - con_mean;
    res = 2 * ct.max_y * ct.max_y + effect_tr * effect_tr  -  2 *  effect_tr * effect_te;

    return res;
}

// set temporarily as CTH_xpred
double userH_xpred(double *y, double wt, double treatment, double tr_mean,
                 double con_mean, double trs, double cons, double alpha, 
                 double xtrain_to_est_ratio, double propensity) {
    double res;
    double tr_var;
    double con_var;
    double tmp;
    if (treatment == 0) {
        // con
        con_var = wt * (y[0] - con_mean) *  (y[0] - con_mean);
        tmp = con_var / ((1 - propensity) * cons);
    } else {
        // tr
        tr_var = wt * (y[0] - tr_mean) * (y[0] - tr_mean);
        tmp = tr_var / (propensity * trs);
    } 
    double effect = tr_mean - con_mean;
    
    res = 4 * ct.max_y * ct.max_y - alpha *  effect * effect + (1 + xtrain_to_est_ratio / (ct.NumXval - 1)) 
        * (1 - alpha) *  tmp; 
    
    return res;
}

// set temporarily as CTA_xpred
double userA_xpred(double *y, double wt, double treatment, double tr_mean, double con_mean,
                 double tree_tr_mean, double tree_con_mean, double alpha) {
    double res;
    double effect_tr = tree_tr_mean - tree_con_mean;
    double effect_te = tr_mean - con_mean;
//    res = 2 * ct.max_y * ct.max_y + effect_tr * effect_tr  -  2 *  effect_tr * effect_te;
	//res = 2 * ct.max_y * ct.max_y + abs(effect_te)*(1.0 - sign(effect_tr) * sign(effect_te)) / 2.0;
	//(x > 0) ? 1 : ((x < 0) ? -1 : 0)
 //////  printf("userA pred abs fn\n");
   res = 2 * ct.max_y * ct.max_y + .99*(abs(effect_te)*(1.0 - ((effect_tr > 0) ? 1 : ((effect_tr < 0) ? -1 : 0)) * ((effect_te > 0) ? 1 : ((effect_te < 0) ? -1 : 0))) / 2.0) +
    .01*(effect_tr * effect_tr  -  2 *  effect_tr * effect_te);
    return res;
}

// policyH
double policyH_xpred(double *y, double wt, double treatment, double tr_mean,
                   double con_mean, double trs, double cons, double alpha, 
                   double xtrain_to_est_ratio, double propensity) {
  double res;
  double tr_var;
  double con_var;
  double tmp;
  if (treatment == 0) {
    // con
    con_var = wt * (y[0] - con_mean) *  (y[0] - con_mean);
    tmp = con_var / ((1 - propensity) * cons);
  } else {
    // tr
    tr_var = wt * (y[0] - tr_mean) * (y[0] - tr_mean);
    tmp = tr_var / (propensity * trs);
  } 
  double effect = tr_mean - con_mean;
  
  res = 4 * ct.max_y * ct.max_y - alpha *  effect * effect + (1 + xtrain_to_est_ratio / (ct.NumXval - 1)) 
    * (1 - alpha) *  tmp; 
  
  return res;
}

// policyA
double policyA_xpred(double *y, double wt, double treatment, double tr_mean, double con_mean,
                   double tree_tr_mean, double tree_con_mean, double alpha, double gamma) {
  double res;
  double effect_tr = tree_tr_mean - tree_con_mean;
  double effect_te = tr_mean - con_mean;
  //    res = 2 * ct.max_y * ct.max_y + effect_tr * effect_tr  -  2 *  effect_tr * effect_te;
  //res = 2 * ct.max_y * ct.max_y + abs(effect_te)*(1.0 - sign(effect_tr) * sign(effect_te)) / 2.0;
  //(x > 0) ? 1 : ((x < 0) ? -1 : 0)
  //////  printf("userA pred abs fn\n");
  res = 2 * ct.max_y * ct.max_y + gamma*(abs(effect_te)*(1.0 - ((effect_tr > 0) ? 1 : ((effect_tr < 0) ? -1 : 0)) * ((effect_te > 0) ? 1 : ((effect_te < 0) ? -1 : 0))) / 2.0) +
            (1-gamma)*(effect_tr * effect_tr  -  2 *  effect_tr * effect_te);
  return res;
}
