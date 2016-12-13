/*
 * This rundown functions You may change it to make it compatibel with other splitting funcitons.
 *
 */
#include "causalTree.h"
#include "node.h"
#include "causalTreeproto.h"

#ifdef NAN
/* NAN is supported */
#endif

void
policyA_rundown(pNode tree, int obs, double *cp, double *xpred, double *xtemp, int k, double alpha, double gamma)
{
    int i, obs2 = (obs < 0) ? -(1 + obs) : obs;
    int my_leaf_id;
    pNode otree =  tree;
    pNode otree_tmp = tree;
    pNode tree_tmp = tree;
    
    int opnumber = 0;
    int j, j2, s;
    int tmp_obs, tmp_id;
    //double tr_mean, con_mean;
    double *tr_mean, *con_mean;
    double consums, trsums, cons, trs;
    double *consumsj, *trsumsj, *consj, *trsj;
    
    consumsj = (double *) ALLOC(ct.ntreats, sizeof(double));
    consj = (double *) ALLOC(ct.ntreats, sizeof(double));
    trsumsj = (double *) ALLOC(ct.ntreats, sizeof(double));
    trsj = (double *) ALLOC(ct.ntreats, sizeof(double));
    tr_mean = (double *) ALLOC(ct.ntreats, sizeof(double));
    con_mean = (double *) ALLOC(ct.ntreats, sizeof(double));
    
    /*
     * Now, repeat the following: for the cp of interest, run down the tree
     *   until I find a node with smaller complexity.  The parent node will
     *   not have collapsed, but this split will have, so this is my
     *   predictor.
     */
    for (i = 0; i < ct.num_unique_cp; i++) {
        cons = 0.;
        trs = 0.;
        consums = 0.;
        trsums = 0.;
        
        //consj =0.0;
        //trsj=0.0;
        //consumsj=0.0;
        //trsumsj=0.0;
        memset(consj, 0, ct.ntreats);
        memset(consumsj, 0, ct.ntreats);
        memset(trsj, 0, ct.ntreats);
        memset(trsumsj, 0, ct.ntreats);
        
        while (cp[i] < tree->complexity) {
	        tree = branch(tree, obs);
	        if (tree == 0)
		        goto oops;
	        otree = tree;
	    }
	    xpred[i] = tree->response_est[0];
        // now find other samples in the same leaf;
        my_leaf_id = tree->id;
        
        
        for (s = k; s < ct.n; s++) {
            tree_tmp = otree_tmp;
            j = ct.sorts[0][s];
            // test: 
           
            tmp_obs = (j < 0) ? -(1 + j) : j;
            while (cp[i] < tree_tmp->complexity) {
                tree_tmp = branch(tree_tmp, tmp_obs);
            }
            tmp_id = tree_tmp->id;
          //change the cons, consums, trs, trsums: make multi
            if (tmp_id == my_leaf_id) 
            /*{
                if (ct.treatment[j] == 0) {
                    cons += ct.wt[j];
                    consums += *ct.ydata[j] * ct.wt[j];
                } else {
                    trs += ct.wt[j];
                    trsums += *ct.ydata[j] * ct.wt[j];
                }
            }*/
            
            for(j2=0;j2<ct.ntreats;j2++)
            {
              
              consj[j2]+=(ct.treatment[j]!=j2);
              consumsj[j2]+=*ct.ydata[j]*ct.wt[j]*(ct.treatment[j]!=j2);
              trsj[j2]+=(ct.treatment[j]==j2);
              trsumsj[j2]+=*ct.ydata[j]*ct.wt[j]*(ct.treatment[j]==j2);
            }
        }
        //change: make multi for tr_mean, con_mean
        //calculate tr_mean and con_mean vectors: xtreatmean and xcontrolmean are already vectors in node.h
        
        /*
        if (trs == 0) {
            // want to trace back to tree->parent for tr_mean;
            tr_mean = tree->parent->xtreatMean[0];
        } else {
            tr_mean = trsumsj / trsj;
            tree->xtreatMean[0] = tr_mean;
        }
        
        if (cons == 0) {
            con_mean = tree->parent->xcontrolMean[0];
        } else {
            con_mean = consumsj / consj;
            tree->xcontrolMean[0] = con_mean;
        }
        
        double tree_tr_mean = tree->treatMean[0];
        double tree_con_mean = tree->controlMean[0];
        */
        for(j2=0;j2<ct.ntreats;j2++)
        {
          //revert to parent if not found
          if(consj[j2]==0)
          {
            con_mean[j2] = tree->parent->xcontrolMean[j2];
          }
          else
          {
            con_mean[j2] = consumsj[j2]/consj[j2];
            tree->xcontrolMean[j2] =con_mean[j2];
          }
          
          //revert to parent if not found
          if(trsj[j2]==0)
          {
            tr_mean[j2] = tree->parent->xtreatMean[j2];
          }
          else
          {
            tr_mean[j2] = trsumsj[j2]/trsj[j2];
            tree->xtreatMean[j2] =tr_mean[j2];
          }
          
        }
        // convert tree_tr_mean and tree_con_mean to vectors and pass to fn policyA_pred
      double *tree_tr_mean, *tree_con_mean;
        tree_tr_mean = (double *) ALLOC(ct.ntreats, sizeof(double));
        tree_con_mean = (double *) ALLOC(ct.ntreats, sizeof(double));
        
        for(j2=0;j2<ct.ntreats;j2++)
        {
          tree_tr_mean[j2]=tree->treatMean[j2];
          tree_con_mean[j2]=tree->controlMean[j2];
        }
        
        
        //xtemp[i] = (*ct_xeval)(ct.ydata[obs2], ct.wt[obs2], ct.treatment[obs2], 
        //            tr_mean, con_mean, tree_tr_mean, tree_con_mean, alpha, gamma);
        
        xtemp[i] = (*ct_xeval_multi)(ct.ydata[obs2], ct.wt[obs2], ct.treatment[obs2], 
                    tr_mean, con_mean, tree_tr_mean, tree_con_mean, alpha, gamma);
    }

    return;

oops:;
    if (ct.usesurrogate < 2) {  /* must have hit a missing value */
	for (i = 0; i < ct.num_unique_cp; i++)
	    xpred[i] = otree->response_est[0];

	xtemp[i] = (*ct_xeval)(ct.ydata[obs2], ct.wt[obs2], ct.treatment[obs2], tr_mean, con_mean);
	
	Rprintf("oops number %d.\n", opnumber++);
  return;
    }
    warning("Warning message--see rundown.c");
}
