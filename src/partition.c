/*
 * The main workhorse of the recursive partitioning module.  When called
 *   with a node, it partitions it and then calls itself to partition the
 *   children it has created.
 * If the node is not splittable (too few people, or complexity is too small)
 *   it simply returns.  The routine may not be able to discover that the
 *   complexity is too small until after the children have been partitioned,
 *   so it needs to check this at the end.
 * The vector who[n] indexes which observations are in this node, to speed
 *   up the routine.
 */
#include <R.h>
#include <stdio.h>
#include "causalTree.h"
#include "node.h"
#include "causalTreeproto.h"



int
partition(int nodenum, pNode splitnode, double *sumrisk, double *sumrisk_multi, int n1, int n2,
          int minsize, int split_Rule, double alpha, int bucketnum, int bucketMax,
          double train_to_est_ratio)
{
    pNode me;
    double tempcp;
    int i, j, k;
    double tempcp2;
    double left_risk, right_risk;
    int left_split, right_split;
    double twt, ttr;
    int nleft, nright;
    int n;
    int min_node_size = minsize;
    FILE* fptr;
    double tempcp_multi[20];
    
    me = splitnode;
    n = n2 - n1;                /* total number of observations */
    me->id = nodenum;
    
//#ifdef DEBUG
    
    //fptr=fopen("C:\\Users\\vikasr\\Documents\\debug_text.txt","w");
    //fprintf(fptr,"test print\n");
    //fclose(fptr);
    R_FlushConsole();
    //Rprintf("test print\n");
    //R_ShowMessage("R_show_message\n");
    
//#endif
    
    if (nodenum > 1) {
        twt = 0;
        ttr = 0;
	    k = 0;
	    for (i = n1; i < n2; i++) {
	      j = ct.sorts[0][i]; /* any variable would do, use first */
	      if (j < 0)
		      j = -(1 + j);   /* if missing, value = -(1+ true index) */
	      ct.wtemp[k] = ct.wt[j];
          ct.trtemp[k] = ct.treatment[j];
	      ct.ytemp[k] = ct.ydata[j];
	      twt += ct.wt[j];
          ttr += ct.treatment[j] * ct.wt[j];
	      k++;
	    }
	    if (split_Rule == 1) {
	        // tot
	        (*ct_eval) (n, ct.ytemp, me->response_est, me->controlMean, me->treatMean, 
          &(me->risk), ct.wtemp, ct.trtemp, ct.max_y, ct.propensity);
	    } else if (split_Rule == 2) {
	        // ct
	        (*ct_eval) (n, ct.ytemp, me->response_est, me->controlMean, me->treatMean, 
          &(me->risk), ct.wtemp, ct.trtemp, ct.max_y, alpha, train_to_est_ratio);
	    } else if (split_Rule == 3) {
	        // fit
	        (*ct_eval) (n, ct.ytemp, me->response_est, me->controlMean, me->treatMean,
          &(me->risk), ct.wtemp, ct.trtemp, ct.max_y, alpha, train_to_est_ratio);
	    } else if (split_Rule == 4) {
	        //tstats
	        (*ct_eval) (n, ct.ytemp, me->response_est, me->controlMean, me->treatMean, 
          &(me->risk), ct.wtemp, ct.trtemp, ct.max_y, alpha, train_to_est_ratio);
	    } else if (split_Rule == 5) {
	        // totD
	        (*ct_eval) (n, ct.ytemp, me->response_est, me->controlMean, me->treatMean,
          &(me->risk), ct.wtemp, ct.trtemp, ct.max_y, ct.propensity);
	    } else if (split_Rule == 6) {
	        // CTD
	        (*ct_eval) (n, ct.ytemp, me->response_est, me->controlMean, me->treatMean, 
          &(me->risk), ct.wtemp, ct.trtemp, ct.max_y, alpha, train_to_est_ratio);
	    } else if (split_Rule == 7) {
	        //fitD
	        (*ct_eval) (n, ct.ytemp, me->response_est, me->controlMean, me->treatMean, 
          &(me->risk), ct.wtemp, ct.trtemp, ct.max_y, alpha, train_to_est_ratio);
	    } else if (split_Rule == 8) {
	        //tstatsD
	        (*ct_eval) (n, ct.ytemp, me->response_est, me->controlMean, me->treatMean, 
          &(me->risk), ct.wtemp, ct.trtemp, ct.max_y, alpha, train_to_est_ratio);
	    } else if (split_Rule == 9) {
	        // user (temporarily set as CT)
	        (*ct_eval) (n, ct.ytemp, me->response_est, me->controlMean, me->treatMean, 
          &(me->risk), ct.wtemp, ct.trtemp, ct.max_y, alpha, train_to_est_ratio);
	    } else if (split_Rule == 10) {
	        // userD (temporarily set as CTD)
	        (*ct_eval) (n, ct.ytemp, me->response_est, me->controlMean, me->treatMean, 
          &(me->risk), ct.wtemp, ct.trtemp, ct.max_y, alpha, train_to_est_ratio);
	    }else if (split_Rule == 11) {
	      // change the call to policy only, make a different function table
	      //is this needed? probably, since we do not wanted to change function signs for other options
	      // also look at pNode (me) vars below
	      // policy (temporarily set as CTD)
	      (*ct_eval_multi) (n, ct.ytemp, me->response_est, me->controlMean, me->treatMean, 
        &(me->risk_multi), ct.wtemp, ct.trtemp, ct.max_y, alpha, train_to_est_ratio);
	    }else if (split_Rule == 12) {
	      // policyD (temporarily set as CTD)
	      (*ct_eval) (n, ct.ytemp, me->response_est, me->controlMean, me->treatMean, 
        &(me->risk), ct.wtemp, ct.trtemp, ct.max_y, alpha, train_to_est_ratio);
	    }

	    me->num_obs = n;
	    me->sum_wt = twt;
        me->sum_tr = ttr;
        if(split_Rule==11)
        {
          tempcp_multi = me->risk_multi;
          tempcp = me->risk_multi[0];
        }
        else
        {
	         tempcp = me->risk;
        }
        //how do we compare now for vector complexity_multi? run a loop here
        //if(split_Rule==11)
        if(0)
        {
          int tmp1 = 0;
          for(tmp1=0;tmp1<ct.ntreats;tmp1++)
            if(tempcp_multi[tmp1]>me->complexity_multi[tmp1])
              tempcp_multi[tmp1]=me->complexity_multi[tmp1];
        }
      else{
	    if (tempcp > me->complexity)
	      tempcp = me->complexity;
        }
    } else
    {
	    if(split_Rule==11)
	    {
	     tempcp_multi = me->risk_multi;
	     tempcp = me->risk_multi[0]; 
	    }
	    else
	    {
	    tempcp = me->risk; 
	    }
    }

    /*
     * Can I quit now ?
     */
  
  //if(split_Rule==11)
  if(0)
  {
    
    //check on vector
    //make multi sum for tempcp_multi
    int tsum_var =0, sumc = 0 ;
    for(tsum_var=0;tsum_var<ct.ntreats;tsum_var++)
    {
      if(tempcp_multi[tsum_var]<=ct.alpha) //ct.alpha_multi[tsum_var]) //check this (do we need alpha_multi?): tbd
        sumc++;
    }
    if(me->num_obs < ct.min_split || sumc >= ntreats || nodenum > ct.maxnode)
    {
      //multi version
      me->complexity = ct.alpha;
      *sumrisk = me->risk;
      for(tsum_var=0;tsum_var<ct.ntreats;tsum_var++)
      {
      me->complexity_multi[tsum_var] = ct.alpha; //ct.alpha_multi[tsum_var];
      sumrisk_multi[tsum_var] = me->risk_multi;//need to fix this up in causaltree.c ('temp' vector version)
      }
      me->leftson = (pNode)  NULL;
      me->rightson = (pNode) NULL;
      me->primary = (pSplit) NULL;
      me->surrogate = (pSplit) NULL;
      return 0;
      
    }
  }
  else{
    if (me->num_obs < ct.min_split || tempcp <= ct.alpha || nodenum > ct.maxnode) {
        
        if(split_Rule==11)
        {
          me->complexity = ct.alpha_multi[0];
          *sumrisk = me->risk_multi[0];  
        }
        else
        {
        me->complexity = ct.alpha;
  	    *sumrisk = me->risk;
        }
      
	/*
	 * make sure the split doesn't have random pointers to somewhere
	 * i.e., don't trust that whoever allocated memory set it to zero
	 */
	    me->leftson = (pNode)  NULL;
	    me->rightson = (pNode) NULL;
	    me->primary = (pSplit) NULL;
	    me->surrogate = (pSplit) NULL;
	    return 0;
    }
  }
    /*
     * Guess I have to do the split
     */
    
    bsplit(me, n1, n2, min_node_size, split_Rule, alpha, bucketnum, bucketMax, train_to_est_ratio);
    if(split_Rule!=11)
    {
    if (!me->primary) {
	/*
	 * This is rather rare -- but I couldn't find a split worth doing
	 */
	    me->complexity = ct.alpha;
	    me->leftson = (pNode) NULL;
	    me->rightson = (pNode) NULL;
	    me->primary = (pSplit) NULL;
	    me->surrogate = (pSplit) NULL;
	    *sumrisk = me->risk;
	    return 0;
    }
    }
#ifdef DEBUG
    print_tree(me, 4);
#endif
    if (ct.maxsur > 0)
	surrogate(me, n1, n2);
    else
	me->surrogate = (pSplit) NULL;
    
    if(split_Rule!=11)
    nodesplit(me, nodenum, n1, n2, &nleft, &nright);

    /*
     * split the leftson
     */
    me->leftson = (pNode) CALLOC(1, nodesize);
    (me->leftson)->parent = me;
    //tbd: tempcp_multi
    (me->leftson)->complexity = tempcp - ct.alpha;
    int j111;
    if(split_Rule==11)
    {
    for(j111=0;j111<ct.ntreats;j111++)
    {
      (me->leftson)->complexity_multi[j111] = tempcp_multi[j111] -ct.alpha;
    }
    }
    left_split = partition(2 * nodenum, me->leftson, &left_risk, left_risk_multi,n1, n1 + nleft,
                           min_node_size, split_Rule, alpha, bucketnum, bucketMax,
                           train_to_est_ratio);

    /*
     * Update my estimate of cp, and split the right son.
     */
    
    //for vectors/split_rule 11, what about tempcp and tempcp2 multi versions?
 if(split_Rule==11)
 {
  int j11;
   //if(split_Rule==11)
   //{
     for(j111=0;j111<ct.ntreats;j111++)
     {
       (me->rightson)->complexity_multi[j111] = tempcp_multi[j111] -ct.alpha;
     }
   //}
   
   //vector version: tbd
   for(j11=0;j11<ct.ntreats;j11++)
   {
  tempcp = (me->risk - left_risk) / (left_split + 1);
   tempcp2 = (me->risk - (me->leftson)->risk);
   if (tempcp < tempcp2)
     tempcp = tempcp2;
   if (tempcp > me->complexity)
     tempcp = me->complexity;
   }
   me->rightson = (pNode) CALLOC(1, nodesize);
   (me->rightson)->parent = me;
   //fix this too for the vector version
   for(j11=0;j11<ct.ntreats;j11++)
   (me->rightson)->complexity_multi[j11] = tempcp_multi[j11] - ct.alpha;
   
   //put the right split call outside the ntreats for loop
   right_split = partition(1 + 2 * nodenum, me->rightson, &right_risk,right_risk_multi,
                           n1 + nleft, n1 + nleft + nright, min_node_size, split_Rule, alpha,
                           bucketnum, bucketMax, train_to_est_ratio);
   
   for(j11=0;j11<ct.ntreats;j11++)
   {
   /*
   * Now calculate my actual C.P., which depends on children nodes, and
   *  on grandchildren who do not collapse before the children.
   * The calculation is done assuming that I am the top node of the
   *  whole tree, an assumption to be fixed up later.
   */
   //use sums for vector version
   int j2;
  for(j2=0;j2<ct.ntreats;j2++)
  {
   tempcp+ = (me->risk_multi[j2] - (left_risk[j2] + right_risk[j2])) /
   (left_split[j2] + right_split[j2] + 1);
  }
   /* Who goes first -- minimum of tempcp, leftson, and rightson */
   if ((me->rightson)->complexity > (me->leftson)->complexity) {
     if (tempcp > (me->leftson)->complexity) {
       /* leftson collapses first */
       left_risk = (me->leftson)->risk;
       left_split = 0;
       
       tempcp = (me->risk - (left_risk + right_risk)) /
       (left_split + right_split + 1);
       if (tempcp > (me->rightson)->complexity) {
         /* right one goes too */
         right_risk = (me->rightson)->risk;
         right_split = 0;
       }
     }
   } else if (tempcp > (me->rightson)->complexity) {
     /* right hand child goes first */
     right_split = 0;
     right_risk = (me->rightson)->risk;
     
     tempcp = (me->risk - (left_risk + right_risk)) /
     (left_split + right_split + 1);
     if (tempcp > (me->leftson)->complexity) {
       /* left one goes too */
       left_risk = (me->leftson)->risk;
       left_split = 0;
     }
   }
   
   me->complexity = (me->risk - (left_risk + right_risk)) /
   (left_split + right_split + 1);
   
   
   if (me->complexity <= ct.alpha) {
     /*
     * All was in vain!  This node doesn't split after all.
     */
     free_tree(me, 0);
     *sumrisk = me->risk;
     for (i = n1; i < n2; i++) {
       j = ct.sorts[0][i];
       if (j < 0)
         j = -(1 + j);
       ct.which[j] = nodenum;      /* revert to the old nodenumber */
     }
     return 0;               /* return # of splits */
   } else {
     *sumrisk = left_risk + right_risk;
     return left_split + right_split + 1;
   }
  }
 }
 else
 {
    tempcp = (me->risk - left_risk) / (left_split + 1);
    tempcp2 = (me->risk - (me->leftson)->risk);
    if (tempcp < tempcp2)
	tempcp = tempcp2;
    if (tempcp > me->complexity)
	tempcp = me->complexity;

    me->rightson = (pNode) CALLOC(1, nodesize);
    (me->rightson)->parent = me;
    (me->rightson)->complexity = tempcp - ct.alpha;
    right_split = partition(1 + 2 * nodenum, me->rightson, &right_risk,
  		    n1 + nleft, n1 + nleft + nright, min_node_size, split_Rule, alpha,
  		    bucketnum, bucketMax, train_to_est_ratio);


    /*
     * Now calculate my actual C.P., which depends on children nodes, and
     *  on grandchildren who do not collapse before the children.
     * The calculation is done assuming that I am the top node of the
     *  whole tree, an assumption to be fixed up later.
     */
    tempcp = (me->risk - (left_risk + right_risk)) /
	  (left_split + right_split + 1);
    /* Who goes first -- minimum of tempcp, leftson, and rightson */
    if ((me->rightson)->complexity > (me->leftson)->complexity) {
      if (tempcp > (me->leftson)->complexity) {
	    /* leftson collapses first */
	      left_risk = (me->leftson)->risk;
	      left_split = 0;

	      tempcp = (me->risk - (left_risk + right_risk)) /
        (left_split + right_split + 1);
	      if (tempcp > (me->rightson)->complexity) {
		/* right one goes too */
		      right_risk = (me->rightson)->risk;
		      right_split = 0;
	      }
	    }
    } else if (tempcp > (me->rightson)->complexity) {
	/* right hand child goes first */
	  right_split = 0;
	  right_risk = (me->rightson)->risk;

	  tempcp = (me->risk - (left_risk + right_risk)) /
	    (left_split + right_split + 1);
	  if (tempcp > (me->leftson)->complexity) {
	    /* left one goes too */
	    left_risk = (me->leftson)->risk;
	    left_split = 0;
      }
    }
    
    me->complexity = (me->risk - (left_risk + right_risk)) /
	(left_split + right_split + 1);

    
    if (me->complexity <= ct.alpha) {
	/*
	 * All was in vain!  This node doesn't split after all.
	 */
	free_tree(me, 0);
	*sumrisk = me->risk;
	for (i = n1; i < n2; i++) {
	    j = ct.sorts[0][i];
	    if (j < 0)
		j = -(1 + j);
	    ct.which[j] = nodenum;      /* revert to the old nodenumber */
	}
	return 0;               /* return # of splits */
    } else {
	*sumrisk = left_risk + right_risk;
	return left_split + right_split + 1;
    }
 }
}
