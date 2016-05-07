/*
 * The error function for anova splitting
 */
 
double
anovapred(double *y, double wt, double treatment, double *yhat, double p) 
{
    double temp;
    if (treatment == 1)  temp = y[0] / p;
    else temp = - y[0] / (1 - p);

    temp -= *yhat;
    return temp * temp * wt;
}
