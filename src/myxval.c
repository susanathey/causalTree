#include <math.h>
#include "causalTree.h"
#include "node.h"
#include "causalTreeproto.h"
#include "func_table.h"

void
myxval(int n_xval, CpTable cptable_head, int *x_grp, int maxcat, char **errmsg, 
       int minsize, int *savesort, int split_Rule,
     int crossmeth, double split_alpha, double cv_alpha, int bucketnum, int bucketMax, double gamma)
{
    int i, j, k, ii, jj;
    int last;
    int xgroup;
    double *xtemp, *xpred, *xpred2;
    int *savew;
    double *cp;
    double alphasave;
    pNode xtree;
    CpTable cplist;
    double temp;
    double temp_multi[20];
    double old_wt, total_wt;
    int neighbor; // nearest neighbor number
    alphasave = ct.alpha;
    double xtrain_to_est_ratio = 0.;
    

    /*
     * Allocate a set of temporary arrays
     */
    xtemp = (double *) CALLOC(4 * ct.num_unique_cp, sizeof(double));
    xpred = xtemp + ct.num_unique_cp;
    xpred2 = xpred + ct.num_unique_cp;
    cp = xpred2 + ct.num_unique_cp;
    savew = (int *) CALLOC(ct.n, sizeof(int));
    for (i = 0; i < ct.n; i++)
        savew[i] = ct.which[i]; /* restore at the end */


    // cp[0] = 10 * cptable_head->cp;      /* close enough to infinity */
    cp[0] = 10000 * cptable_head->cp;
    for (cplist = cptable_head, i = 1; i < ct.num_unique_cp;cplist = cplist->forward, i++) {  
        cp[i] = sqrt(cplist->cp * (cplist->forward)->cp);
    }

    total_wt = 0;
    for (i = 0; i < ct.n; i++)
        total_wt += ct.wt[i];
    old_wt = total_wt;

    /*
     * do the validations
     */
    
    k = 0;                      /* -Wall */
    for (xgroup = 0; xgroup < n_xval; xgroup++) {
        /*
         * restore ct.sorts, with the data for this run at the top
         * this requires one pass per variable
         */
        for (j = 0; j < ct.nvar; j++) {
            k = 0;
            for (i = 0; i < ct.n; i++) {
                ii = savesort[j * ct.n + i];
                if (ii < 0)
                    ii = -(1 + ii);     /* missings move too */
                if (x_grp[ii] != xgroup + 1) { 
                    /*
                     * this obs is left in --
                     *  copy to the front half of ct.sorts
                     */
                    ct.sorts[j][k] = savesort[j * ct.n + i]; // the reason to store in savesort
                    k++;
                }
            } 
        }

        /*
         *  Fix up the y vector, and save a list of "left out" obs *   in
         * the tail, unused end of ct.sorts[0][i];
         */
        last = k;

        k = 0;
        temp = 0;
        for (i = 0; i < ct.n; i++) {
            ct.which[i] = 1;    /* everyone starts in group 1 */
            if (x_grp[i] == xgroup + 1) {
                ct.sorts[0][last] = i;
                last++;
            } else {
                ct.ytemp[k] = ct.ydata[i];
                ct.wtemp[k] = ct.wt[i];
                ct.trtemp[k] = ct.treatment[i];
                temp += ct.wt[i];
                k++;
            }
        }

        /* rescale the cp */
        for (j = 0; j < ct.num_unique_cp; j++) 
            cp[j] *= temp / old_wt;
        ct.alpha *= temp / old_wt;
        old_wt = temp;
        
        int iii=0;
        for (iii = 0; iii < ct.ntreats; iii++)
          temp_multi[iii] =temp;
        
        /*
         * partition the new tree
         */
        xtree = (pNode) CALLOC(1, nodesize);
        xtree->num_obs = k;
       
       if (split_Rule == 11){
         ct_init = split_func_table_multi[0].init_split_multi;
         ct_choose = split_func_table_multi[0].choose_split_multi;
         ct_eval = split_func_table_multi[0].eval_multi;
         ct_xeval = cv_func_table_multi[0].xeval_multi;
       }
        
        (*ct_init) (k, ct.ytemp, maxcat, errmsg, &temp, 2, ct.wtemp, ct.trtemp, 
         bucketnum, bucketMax, &xtrain_to_est_ratio);
        
        if (split_Rule == 1) {
            //tot:
            (*ct_eval) (k, ct.ytemp, xtree->response_est, xtree->controlMean, xtree->treatMean,
             &(xtree->risk), ct.wtemp, ct.trtemp, ct.max_y, ct.propensity);
        } else if (split_Rule == 2) {
            // ct:
            (*ct_eval) (k, ct.ytemp, xtree->response_est, xtree->controlMean, xtree->treatMean,
             &(xtree->risk), ct.wtemp, ct.trtemp, ct.max_y, split_alpha, xtrain_to_est_ratio);
        } else if (split_Rule == 3) {
            // fit
            (*ct_eval) (k, ct.ytemp, xtree->response_est, xtree->controlMean, xtree->treatMean, 
             &(xtree->risk), ct.wtemp, ct.trtemp, ct.max_y, split_alpha, xtrain_to_est_ratio);
        } else if (split_Rule == 4) {
            // tstats
            (*ct_eval) (k, ct.ytemp, xtree->response_est, xtree->controlMean, xtree->treatMean, 
             &(xtree->risk), ct.wtemp, ct.trtemp, ct.max_y, split_alpha, xtrain_to_est_ratio);
        } else if (split_Rule == 5) {
            // totD
            (*ct_eval) (k, ct.ytemp, xtree->response_est, xtree->controlMean, xtree->treatMean, 
             &(xtree->risk), ct.wtemp, ct.trtemp, ct.max_y, ct.propensity);
        } else if (split_Rule == 6) {
            // CTD
            (*ct_eval) (k, ct.ytemp, xtree->response_est, xtree->controlMean, xtree->treatMean,
             &(xtree->risk), ct.wtemp, ct.trtemp, ct.max_y, split_alpha, xtrain_to_est_ratio);
        } else if (split_Rule == 7) {
            // fitD:
            (*ct_eval) (k, ct.ytemp, xtree->response_est, xtree->controlMean, xtree->treatMean, 
             &(xtree->risk), ct.wtemp, ct.trtemp, ct.max_y, split_alpha, xtrain_to_est_ratio);
        } else if (split_Rule == 8) {
            //tstatsD:
            (*ct_eval) (k, ct.ytemp, xtree->response_est, xtree->controlMean, xtree->treatMean, 
             &(xtree->risk), ct.wtemp, ct.trtemp, ct.max_y, split_alpha, bucketnum, bucketMax,
             xtrain_to_est_ratio);
        } else if (split_Rule == 9) {
            // user (set temporarily as CT)
            (*ct_eval) (k, ct.ytemp, xtree->response_est, xtree->controlMean, xtree->treatMean,
             &(xtree->risk), ct.wtemp, ct.trtemp, ct.max_y, split_alpha, xtrain_to_est_ratio);
        } else if (split_Rule == 10) {
            // userD (set temporarily as CTD)
            (*ct_eval) (k, ct.ytemp, xtree->response_est, xtree->controlMean, xtree->treatMean,
             &(xtree->risk), ct.wtemp, ct.trtemp, ct.max_y, split_alpha, xtrain_to_est_ratio);
        }else if (split_Rule == 11) {
          // policy
          //(*ct_eval) (k, ct.ytemp, xtree->response_est, xtree->controlMean, xtree->treatMean,
           //&(xtree->risk), ct.wtemp, ct.trtemp, ct.max_y, split_alpha, xtrain_to_est_ratio);
           (*ct_eval_multi) (k, ct.ytemp, xtree->response_est, xtree->controlMean, xtree->treatMean, 
            &(xtree->risk_multi), ct.wtemp, ct.trtemp, ct.max_y, split_alpha, xtrain_to_est_ratio);
        }else if (split_Rule == 12) {
          // policyD
          (*ct_eval) (k, ct.ytemp, xtree->response_est, xtree->controlMean, xtree->treatMean,
           &(xtree->risk), ct.wtemp, ct.trtemp, ct.max_y, split_alpha, xtrain_to_est_ratio);
        }
        
        if (split_Rule == 11){
          xtree->risk = xtree->risk_multi[0];
          xtree->complexity = xtree->risk;
        }
        else
        xtree->complexity = xtree->risk;


        partition(1, xtree, &temp, &temp_multi, 0, k, minsize, split_Rule, split_alpha, bucketnum, bucketMax, 
                  xtrain_to_est_ratio);
        

        fix_cp(xtree, xtree->complexity);

        /*
         * run the extra data down the new tree
         */
      
        for(i = k; i < ct.n; i++) {
            j = ct.sorts[0][i];
            if (crossmeth == 1) {
                // tot
                totrundown(xtree, j, cp, xpred, xtemp);
            } else if (crossmeth == 2) {
                // matching
                neighbor = findNeighbor(j, k);
                matching_rundown(xtree, j, neighbor, cp, xpred, xpred2, xtemp);
            } else if (crossmeth == 3) {
                //fit-honest: alpha need to change to cv.alpha, discriminate with split.alpha
                fitH_rundown(xtree, j, cp, xpred, xtemp, k, cv_alpha, xtrain_to_est_ratio);
            } else if (crossmeth == 4) {
                // fit-adaptive:
                fitA_rundown(xtree, j, cp, xpred, xtemp, k);
                
            } else if (crossmeth == 5) {
                //CT- honest
                CTH_rundown(xtree, j, cp, xpred, xtemp, k, cv_alpha, xtrain_to_est_ratio, ct.propensity);
                
            } else if (crossmeth == 6) {
                //CT- dishonest
                CTA_rundown(xtree, j, cp, xpred, xtemp, k, cv_alpha);
            } else if (crossmeth == 7) {
                // user - honest (set as CT - honest temporarily)
                userH_rundown(xtree, j, cp, xpred, xtemp, k, cv_alpha, xtrain_to_est_ratio, ct.propensity);
            } else if (crossmeth == 8) {
                // user - dishonest (set as CT - dishonest temporarily)
                userA_rundown(xtree, j, cp, xpred, xtemp, k, cv_alpha);
            }else if (crossmeth == 9) {
              // user - honest (set as CT - honest temporarily)
              policyH_rundown(xtree, j, cp, xpred, xtemp, k, cv_alpha, xtrain_to_est_ratio, ct.propensity);
            }else if (crossmeth == 10) {
              // user - dishonest (set as CT - dishonest temporarily)
              //need to pass ** vector for xtemp
              policyA_rundown(xtree, j, cp, xpred, xtemp, k, cv_alpha, gamma);
            }


#if DEBUG > 1
          if (debug > 1) {
                jj = j + 1;
                Rprintf("\nObs %d, y=%f \n", jj, ct.ydata[j][0]);
           }
#endif
            /* add it in to the risk */
            cplist = cptable_head;
            for (jj = 0; jj < ct.num_unique_cp; jj++) {
                cplist->xrisk += xtemp[jj] * ct.wt[j];
                cplist->xstd += xtemp[jj] * xtemp[jj] * ct.wt[j];
                
#if DEBUG > 1
                //	if (debug > 1)
                Rprintf("  cp=%f, pred=%f, xtemp=%f\n",
                        cp[jj] , xpred[jj], xtemp[jj]);
#endif
                cplist = cplist->forward;
            }
        }

        free_tree(xtree, 1);    // Calloc-ed
        R_CheckUserInterrupt();
    }

    for (cplist = cptable_head; cplist; cplist = cplist->forward) {
        cplist->xstd = sqrt(cplist->xstd -
                cplist->xrisk * cplist->xrisk / total_wt);

    }
    ct.alpha = alphasave;
    for (i = 0; i < ct.n; i++)
        ct.which[i] = savew[i];
    Free(savew);
    Free(xtemp);
}
