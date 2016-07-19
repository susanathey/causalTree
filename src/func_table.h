/*
 * The table of implimented splitting functions
 *
 */

extern int totinit(int n, double *y[], int maxcat, char **error, int *size,
                   int who, double *wt, double *treatment, int bucketnum, int bucketMax, 
                   double *train_to_est_ratio);
    
extern void tot(int n, double *y[], double *x, int nclass,
                  int edge, double *improve, double *split, int *csplit,
                  double myrisk, double *wt, double *treatment, double propensity,int minsize);
extern void totss(int n, double *y[], double *value, double *tr_mean, double *con_mean, double *risk,
                    double *wt, double *treatment, double max_y, double propensity);
extern double totpred(double *y, double wt, double treatment, double *yhat, double p);

extern int totDinit(int n, double *y[], int maxcat, char **error,
             int *size, int who, double *wt, double *treatment,
             int bucketnum, int bucketMax, double *train_to_est_ratio);
extern void totDss(int n, double *y[], double *value, double *tr_mean, double *con_mean, double *risk, double *wt, 
                   double *treatment, double max_y, double propensity);
extern void totD(int n, double *y[], double *x, int nclass, int edge, double *improve, 
                      double *split, int *csplit, double myrisk, double *wt,
                      double *treatment, double propensity, int minsize, int bucketnum,
                      int bucketMax);
extern double totDpred(double *y, double wt, double treatment, double *yhat, double p);

extern int CTinit(int n, double *y[], int maxcat, char **error,
           int *size, int who, double *wt, double *treatment, int bucketnum, 
           int bucketMax, double *train_to_est_ratio);
extern void CTss(int n, double *y[], double *value, double *tr_mean, double *con_mean, double *risk, 
                 double *wt, double *treatment, double max_y, double alpha, double train_to_est_ratio);
extern void CT(int n, double *y[], double *x, int nclass,int edge, double *improve, double *split, 
               int *csplit, double myrisk, double *wt, double *treatment, int minsize, double alpha,
               double train_to_est_ratio);
extern double CTpred(double *y, double wt, double treatment, double *yhat, double propensity);

extern int CTDinit(int n, double *y[], int maxcat, char **error, int *size,
                   int who, double *wt, double *treatment, int bucketnum, int bucketMax,
                   double *train_to_est_ratio);

extern void CTDss(int n, double *y[], double *value, double *tr_mean, double *con_mean, double *risk, double *wt, 
                  double *treatment, double max_y, double alpha, double train_to_est_ratio);
extern void CTD(int n, double *y[], double *x, int nclass,
        int edge, double *improve, double *split, int *csplit,
        double myrisk, double *wt, double *treatment, int minsize, double alpha, int bucketnum,
        int bucketMax, double train_to_est_ratio);
extern double CTDpred(double *y, double wt, double treatment, double *yhat, 
                      double propensity); 

extern int fitinit(int n, double *y[], int maxcat, char **error,
            int *size, int who, double *wt, double *treatment, int bucketnum, 
            int bucketMax, double *train_to_est_ratio);

extern void fitss(int n, double *y[], double *value, double *tr_mean, double *con_mean, 
                  double *risk, double *wt, double *treatment, double max_y,
                  double alpha, double train_to_est_ratio);

extern void fit(int n, double *y[], double *x, int nclass,
                int edge, double *improve, double *split, int *csplit,
                double myrisk, double *wt, double *treatment, int minsize,
                double alpha, double train_to_est_ratio);
extern double fitpred(double *y, double wt, double treatment, double *yhat, double propensity);

extern int fitDinit(int n, double *y[], int maxcat, char **error,
                   int *size, int who, double *wt, double *treatment, int bucketnum, 
                   int bucketMax, double *train_to_est_ratio);

extern void fitDss(int n, double *y[], double *value, double *tr_mean, double *con_mean, 
                  double *risk, double *wt, double *treatment, double max_y,
                  double alpha, double train_to_est_ratio);

extern void fitD(int n, double *y[], double *x, int nclass,
                int edge, double *improve, double *split, int *csplit,
                double myrisk, double *wt, double *treatment, int minsize,
                int bucketnum, int bucketMax, double alpha, double train_to_est_ratio);
extern double fitDpred(double *y, double wt, double treatment, double *yhat, double propensity);

extern int tstatsinit(int n, double *y[], int maxcat, char **error,
                      int *size, int who, double *wt, double *treatment,
                      int bucketnum, int bucketMax, double *train_to_est_ratio);

extern void tstatsss(int n, double *y[], double *value, double *tr_mean, double *con_mean, double *risk,
                     double *wt, double *treatment, double max_y, double alpha, double train_to_est_ratio);
extern void tstats(int n, double *y[], double *x, int nclass,
                   int edge, double *improve, double *split, int *csplit, 
                   double myrisk, double *wt, double *treatment, int minsize, 
                   double alpha, double train_to_est_ratio);
extern double tstatspred(double *y, double wt, double treatment, double *yhat, double propensity);

extern int tstatsDinit(int n, double *y[], int maxcat, char **error,
                int *size, int who, double *wt, double *treatment, 
                int bucketnum, int bucektMax, double *train_to_est_ratio);
extern void tstatsDss(int n, double *y[], double *value, double *tr_mean, double *con_mean, double *risk,
                      double *wt, double *treatment, double max_y, double alpha, double train_to_est_ratio);
extern void tstatsD(int n, double *y[], double *x, int nclass,
            int edge, double *improve, double *split, int *csplit,
            double myrisk, double *wt, double *treatment, int minsize, double alpha, 
            int bucketnum, int bucketMax, double train_to_est_ratio);
extern double tstatsDpred(double *y, double wt, double treatment, double *yhat, double propensity);



extern int userinit(int n, double *y[], int maxcat, char **error,
                  int *size, int who, double *wt, double *treatment, int bucketnum, 
                  int bucketMax, double *train_to_est_ratio);
extern void userss(int n, double *y[], double *value, double *tr_mean, double *con_mean, double *risk, 
                 double *wt, double *treatment, double max_y, double alpha, double train_to_est_ratio);
extern void user(int n, double *y[], double *x, int nclass,int edge, double *improve, double *split, 
               int *csplit, double myrisk, double *wt, double *treatment, int minsize, double alpha,
               double train_to_est_ratio);
extern double userpred(double *y, double wt, double treatment, double *yhat, double propensity);

extern int userDinit(int n, double *y[], int maxcat, char **error, int *size,
                   int who, double *wt, double *treatment, int bucketnum, int bucketMax,
                   double *train_to_est_ratio);

extern void userDss(int n, double *y[], double *value, double *tr_mean, double *con_mean, double *risk, double *wt, 
                  double *treatment, double max_y, double alpha, double train_to_est_ratio);
extern void userD(int n, double *y[], double *x, int nclass,
                int edge, double *improve, double *split, int *csplit,
                double myrisk, double *wt, double *treatment, int minsize, double alpha, int bucketnum,
                int bucketMax, double train_to_est_ratio);
extern double userDpred(double *y, double wt, double treatment, double *yhat, 
                      double propensity); 

///// adding policy
extern int policyinit(int n, double *y[], int maxcat, char **error,
                    int *size, int who, double *wt, double *treatment, int bucketnum, 
                    int bucketMax, double *train_to_est_ratio);
extern void policyss(int n, double *y[], double *value, double *tr_mean, double *con_mean, double *risk, 
                   double *wt, double *treatment, double max_y, double alpha, double train_to_est_ratio);
extern void policy(int n, double *y[], double *x, int nclass,int edge, double *improve, double *split, 
                 int *csplit, double myrisk, double *wt, double *treatment, int minsize, double alpha,
                 double train_to_est_ratio);
extern double policypred(double *y, double wt, double treatment, double *yhat, double propensity);

extern int policyDinit(int n, double *y[], int maxcat, char **error, int *size,
                     int who, double *wt, double *treatment, int bucketnum, int bucketMax,
                     double *train_to_est_ratio);

extern void policyDss(int n, double *y[], double *value, double *tr_mean, double *con_mean, double *risk, double *wt, 
                    double *treatment, double max_y, double alpha, double train_to_est_ratio);
extern void policyD(int n, double *y[], double *x, int nclass,
                  int edge, double *improve, double *split, int *csplit,
                  double myrisk, double *wt, double *treatment, int minsize, double alpha, int bucketnum,
                  int bucketMax, double train_to_est_ratio);
extern double policyDpred(double *y, double wt, double treatment, double *yhat, 
                        double propensity); 


/* ------------------------------------------------------ --------------------------- */
extern double tot_xpred(double *y, double wt, double treatment, double *yhat, double propensity); 
extern double matching_xpred(double *y1, double *y2, double wt1, double wt2, double treatment1,
                                  double treatment2, double pred1, double pred2);
extern double fitH_xpred(double *y, double wt, double treatment, double tr_mean, double con_mean,
                         double trs, double cons, double alpha, double xtrain_to_est_ratio);
extern double fitA_xpred(double *y, double wt, double treatment, double tr_mean, double con_mean);
extern double CTH_xpred(double *y, double wt, double treatment, double tr_mean, double con_mean,
                        double trs, double cons, double alpha, double xtrain_to_est_ratio,
                        double propensity);
extern double CTA_xpred(double *y, double wt, double treatment, double tr_mean, double con_mean,
                 double tree_tr_mean, double tree_con_mean, double alpha);
extern double userH_xpred(double *y, double wt, double treatment, double tr_mean, double con_mean,
                        double trs, double cons, double alpha, double xtrain_to_est_ratio,
                        double propensity);
extern double userA_xpred(double *y, double wt, double treatment, double tr_mean, double con_mean,
                        double tree_tr_mean, double tree_con_mean, double alpha);

///// adding policy
extern double policyH_xpred(double *y, double wt, double treatment, double tr_mean, double con_mean,
                          double trs, double cons, double alpha, double xtrain_to_est_ratio,
                          double propensity);
extern double policyA_xpred(double *y, double wt, double treatment, double tr_mean, double con_mean,
                          double tree_tr_mean, double tree_con_mean, double alpha, double gamma);


//extern double totxeval(int *unique_leaf, int **val_leaf_mat, int cp_id, int t, int *sorts, double *wt, 
//                      double *treatment,  double *y[], double propensity, int k, int nobs, 
//                       double val_sum_wt, int val_count);


    
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
extern int anovainit(int n, double *y[], int maxcat, char **error,
  	     double *parm, int *size, int who, double *wt, double *treatment);


extern void anovass(int n, double *y[], double *value, double *risk,
  	    double *wt, double *treatment, double max_y);


extern void anova(int n, double *y[], double *x, int nclass,
  	  int edge, double *improve, double *split, int *csplit,
		  double myrisk, double *wt, double *treatment, int minsize);

extern double anovapred(double *y, double wt, double treatment, double *yhat, double p);


static struct {
    int (*init_split) ();
    void (*choose_split) ();
    void (*eval) ();
    double (*error) ();
} func_table[] = {
    {anovainit, anova, anovass, anovapred}, // test for TOT:
};
#define NUM_METHODS 1        /* size of the above structure */  

static struct {
    int (*init_split) ();
    void (*choose_split) ();
    void (*eval) ();
    double (*error) ();
} split_func_table[] = {
    {totinit, tot, totss, totpred},
    {CTinit, CT, CTss, CTpred},
    {fitinit, fit, fitss, fitpred},
    {tstatsinit, tstats, tstatsss, tstatspred},
    {totDinit, totD, totDss, totDpred},
    {CTDinit, CTD, CTDss, CTDpred},
    {fitDinit, fitD, fitDss, fitDpred},
    {tstatsDinit, tstatsD, tstatsDss, tstatsDpred},
    {userinit, user, userss, userpred},
    {userDinit, userD, userDss, userDpred},
    {policyinit, policy, policyss, policypred},
    {policyDinit, policyD, policyDss, policyDpred},
};

static struct {

    double (* xeval)();
} cv_func_table[] = {
    {tot_xpred},
    {matching_xpred},
    {fitH_xpred},
    {fitA_xpred},
    {CTH_xpred},
    {CTA_xpred},
    {userH_xpred},
    {userA_xpred},
    {policyH_xpred},
    {policyA_xpred},
};

#define NUM_SPLIT_RULE 12
#define NUM_CROSSMETH 11
