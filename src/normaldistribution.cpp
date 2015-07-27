#include <iostream>
#include <math.h>
#include <vector>
#include "normaldistribution.h"


using namespace std;

double normalpdf(double z)
{
  return ((1.0)/sqrt(2*3.14159265))*exp(-pow(z,2)/2);
}


double normalcdf(double z){
  if(z>6.0){return 1.0;}
  if(z<-6.0){return 0.0;}

  double b1=0.31938153;
  double b2=-0.356563782;
  double b3=1.781477937;
  double b4=-1.821255978;
  double b5=1.330274429;
  double p=0.2316419;
  double c2=0.3989423;

  double a=fabs(z);
  double t=1.0/(1.0+a*p);

  double b=c2*exp((-z)*(z/2.0));
  double n=((((b5*t+b4)*t+b3)*t+b2)*t+b1)*t;
  n=1.0-b*n;
  if(z<0.0) 
    {
      n=1.0-n;
    }
  return n;
}


double inversenormalcdf(double p)
{

  double  A1 = (-3.969683028665376e+01);
  double  A2 =  2.209460984245205e+02;
  double  A3 = (-2.759285104469687e+02);
  double  A4 =  1.383577518672690e+02;
  double  A5 = (-3.066479806614716e+01);
  double  A6 =  2.506628277459239e+00;
  
  double  B1 = (-5.447609879822406e+01);
  double  B2 =  1.615858368580409e+02;
  double  B3 = (-1.556989798598866e+02);
  double  B4 =  6.680131188771972e+01;
  double  B5 = (-1.328068155288572e+01);
  
  double  C1 = (-7.784894002430293e-03);
  double  C2 = (-3.223964580411365e-01);
  double  C3 = (-2.400758277161838e+00);
  double  C4 = (-2.549732539343734e+00);
  double  C5 =  4.374664141464968e+00;
  double  C6 =  2.938163982698783e+00;

  double  D1  = 7.784695709041462e-03;
  double  D2 =  3.224671290700398e-01;
  double  D3 =  2.445134137142996e+00;
  double  D4 =  3.754408661907416e+00;
  
  double P_LOW = 0.02425;
  double P_HIGH = 0.97575;

  double x;
  double q,r,u,e;

  if((0<p)&&(p<P_LOW)){
    q = sqrt(-2*log(p));
    x = (((((C1*q+C2)*q+C3)*q+C4)*q+C5)*q+C6) / ((((D1*q+D2)*q+D3)*q+D4)*q+1);
  }
  
  else{
    if ((P_LOW <= p) && (p <= P_HIGH)){
      q = p - 0.5;
      r = q*q;
      x = (((((A1*r+A2)*r+A3)*r+A4)*r+A5)*r+A6)*q /(((((B1*r+B2)*r+B3)*r+B4)*r+B5)*r+1);
    }
        else{
	  if ((P_HIGH < p)&&(p < 1)){
	    q = sqrt(-2*log(1-p));
	    x = -(((((C1*q+C2)*q+C3)*q+C4)*q+C5)*q+C6) / ((((D1*q+D2)*q+D3)*q+D4)*q+1);
	  }
        }
  }
   
  if(( 0 < p)&&(p < 1)){
  e = 0.5 * erfc(-x/sqrt(2)) - p;
  u = e * sqrt(2*M_PI) * exp(x*x/2);
  x = x - u/(1 + x*u/2);
  }
  
  return x;
}


