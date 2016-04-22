/*
 * The error function for anova splitting
 */
 
 // only used temporarily for TOT cv method
// static double p = 0.5;

double
//anovapred(double *y, double wt, double *yhat)
anovapred(double *y, double wt, double treatment, double *yhat, double p) // pass in ct.which
{
    double temp;
    if (treatment == 1)  temp = y[0] / p;
    else temp = - y[0] / (1 - p);
      
    //double temp = y[0] - *yhat;
    temp -= *yhat;
    
    return temp * temp * wt;
}
