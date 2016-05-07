/*
 * some tool functions for cross validation
 */
#include <math.h>
#include "causalTree.h"
#include "node.h"
#include "causalTreeproto.h"


#ifdef NAN
/* NAN is supported */
#endif
static double INFTY = 9999999999; // infinity

int findNeighbor(int obs, int k) { 
  int i, j, temp, neighbor;
  int obs2 = (obs < 0)? -(1 + obs) : obs;
  double dist, min = INFTY;

 neighbor = 0;
  
  for (i = k; i < ct.n; i++) {
    j = ct.sorts[0][i];
    temp = (j < 0)? -(1 + j) : j;
    if (ct.treatment[temp] != ct.treatment[obs2]) {
      dist = measureDistance(obs2, temp);
      if (dist < min) {
        neighbor = j;
        min = dist;       
      } 
    }       
  }
  return neighbor;
}


double measureDistance(int i, int j) {
  int k;
  double distance = 0;
  for (k = 0; k < ct.nvar; k++) {
    distance += (ct.xdata[k][i] - ct.xdata[k][j]) * (ct.xdata[k][i] - ct.xdata[k][j]) / ct.xvar[k];   
  }
  return distance;
}
