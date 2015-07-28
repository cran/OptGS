#include "normaldistribution.h"
#include "Rmath.h"


using namespace std;

double normalpdf(double z)
{
  return dnorm(z,0,1,0);
}


double normalcdf(double z){
 
  return pnorm(z,0,1,1,0);
}


double inversenormalcdf(double p)
{

  return qnorm(p,0,1,1,0);
}


