//multistagefunctions_v2 provides functions that allow finding the operating characteristics of a multistage sequential trial using z-statistics. If the number of stages is below a certain threshold, integration of the multivariate normal is used. If above, the method of sequential integration is used.


#include <iostream>
#include <math.h>
#include <vector>
#include <time.h>
#include "multistagefunctions_v3.h"
#include "normaldistribution.h"


using namespace std;


double onestagesamplesize(double difference,double sigma,double typeIerror,double typeIIerror,double R)
{

return((pow(inversenormalcdf(1-(typeIerror))+inversenormalcdf(1-typeIIerror),2)*(pow(sigma*sqrt((1+R)/R),2))/(difference*difference)));
}


//information calculates the information given number of individuals, delta, and sigma:

double information(double numberindividuals,double sigma)
{

  return(numberindividuals/(2*sigma*sigma));

}

void converthtophi(vector<vector<double> >& h,vector<vector<double> >& z,vector<double>& phi,vector<double>& parameters,double delta0,double newdelta,double sigma)
{//converthtophi calculates probability of stopping for efficacy at each stage
  int i,j;

 phi.clear();
 
  phi.push_back(1-normalcdf(parameters.at(2)-newdelta*sqrt(parameters.at(0))/sqrt(2*sigma*sigma)));
 
  for(i=1;i<=z.size();i++)
    {
      phi.push_back(0);
     
      for(j=0;j<z.at(i-1).size();j++)
	{
	  
	    phi.at(i)+=(exp((newdelta-delta0)*z.at(i-1).at(j)*sqrt(information(i*parameters.at(0),sigma))-(pow(newdelta,2)-pow(delta0,2))*information(i*parameters.at(0),sigma)/2)*(h.at(i-1).at(j)*(normalcdf((z.at(i-1).at(j)*sqrt(information(i*parameters.at(0),sigma))-parameters.at((i+1)*2)*sqrt(information((i+1)*parameters.at(0),sigma))+newdelta*(information((i+1)*parameters.at(0),sigma)-information(i*parameters.at(0),sigma)))/(sqrt(information((i+1)*parameters.at(0),sigma)-information(i*parameters.at(0),sigma)))))));

 
	}
    
    
    }

}


double expectedsamplesize(vector<double> phi,vector<double> psi,vector<double> parameters)
{

  int i;
  double expectedsamplesize=0;

  for(i=0;i<phi.size();i++)
    {
      expectedsamplesize+=(i+1)*parameters.at(0)*(phi.at(i)+psi.at(i));
    }

  return expectedsamplesize;
}

void converthtopsi(vector<vector<double> >& h,vector<vector<double> >& z,vector<double>& psi,vector<double>& parameters,double delta0,double newdelta,double sigma)
{
  int i,j;

 psi.clear();
 
  psi.push_back(normalcdf(parameters.at(1)-newdelta*sqrt(parameters.at(0))/sqrt(2*sigma*sigma)));
 
  for(i=1;i<=z.size();i++)
    {
      psi.push_back(0);
     
      for(j=0;j<z.at(i-1).size();j++)
	{
	  
	    psi.at(i)+=(exp((newdelta-delta0)*z.at(i-1).at(j)*sqrt(information(i*parameters.at(0),sigma))-(pow(newdelta,2)-pow(delta0,2))*information(i*parameters.at(0),sigma)/2)*(h.at(i-1).at(j)*(1-normalcdf((z.at(i-1).at(j)*sqrt(information(i*parameters.at(0),sigma))-parameters.at((i+1)*2-1)*sqrt(information((i+1)*parameters.at(0),sigma))+newdelta*(information((i+1)*parameters.at(0),sigma)-information(i*parameters.at(0),sigma)))/(sqrt(information((i+1)*parameters.at(0),sigma)-information(i*parameters.at(0),sigma)))))));

	   
	}
    
         }

}





void finddeltaminimax_seq(vector<vector<double> >& h,vector<vector<double> >& z,vector<double>& phi,vector<double>& parameters,double delta0,double newdelta,double sigma,double *deltaminimax,double *maxen)
{

  //finds the delta which gives highest expected sample size

  double lowerdelta=delta0,middelta=newdelta,upperdelta=delta0+2*(newdelta-delta0),tempdelta,loweren,miden,upperen,tempen;
  vector<double> lowerphi,lowerpsi,midphi,midpsi,upperphi,upperpsi,tempphi,temppsi;

converthtophi(h,z,lowerphi,parameters,delta0,lowerdelta,sigma);
converthtophi(h,z,midphi,parameters,delta0,middelta,sigma);
converthtophi(h,z,upperphi,parameters,delta0,upperdelta,sigma);
converthtopsi(h,z,lowerpsi,parameters,delta0,lowerdelta,sigma);
converthtopsi(h,z,midpsi,parameters,delta0,middelta,sigma);
converthtopsi(h,z,upperpsi,parameters,delta0,upperdelta,sigma);
 loweren=expectedsamplesize(lowerphi,lowerpsi,parameters);
 miden=expectedsamplesize(midphi,midpsi,parameters); 
upperen=expectedsamplesize(upperphi,upperpsi,parameters);

//check whether miden is higher than both loweren and upperen

 if(miden<loweren || miden<upperen)
   {
     do
       {
	 lowerdelta-=(newdelta-delta0);
	 upperdelta+=(newdelta-delta0);
	 converthtophi(h,z,lowerphi,parameters,delta0,lowerdelta,sigma);
	 converthtophi(h,z,upperphi,parameters,delta0,upperdelta,sigma);
	 converthtopsi(h,z,lowerpsi,parameters,delta0,lowerdelta,sigma);
	 converthtopsi(h,z,upperpsi,parameters,delta0,upperdelta,sigma);
	 loweren=expectedsamplesize(lowerphi,lowerpsi,parameters);
	 upperen=expectedsamplesize(upperphi,upperpsi,parameters);
	 
       }
     while(miden<loweren || miden<upperen);
   }


//find delta which gives maximum expected sample size:

do
  {
  
    if((upperdelta-middelta)>(middelta-lowerdelta))
      {
	tempdelta=(upperdelta+middelta)/2;
	converthtophi(h,z,tempphi,parameters,delta0,tempdelta,sigma);
	converthtopsi(h,z,temppsi,parameters,delta0,tempdelta,sigma);
	tempen=expectedsamplesize(tempphi,temppsi,parameters);
	if(tempen<miden)
	  {
	    upperen=tempen;
	    upperdelta=tempdelta;
	    upperphi=tempphi;
	    upperpsi=temppsi;
	  }
	else
	  {
	    loweren=miden;
	    lowerdelta=middelta;
	    lowerphi=midphi;
	    lowerpsi=midpsi;
	    miden=tempen;
	    middelta=tempdelta;
	    midphi=tempphi;
	    midpsi=temppsi;
	  }
      }

    else
      {

	tempdelta=(lowerdelta+middelta)/2;
	converthtophi(h,z,tempphi,parameters,delta0,tempdelta,sigma);
	converthtopsi(h,z,temppsi,parameters,delta0,tempdelta,sigma);
	tempen=expectedsamplesize(tempphi,temppsi,parameters);

	if(tempen<miden)
	  {
	    loweren=tempen;
	    lowerdelta=tempdelta;
	    lowerphi=tempphi;
	    lowerpsi=temppsi;
	  }
	else
	  {
	    upperen=miden;
	    upperdelta=middelta;
	    upperphi=midphi;
	    upperpsi=midpsi;
	    miden=tempen;
	    middelta=tempdelta;
	    midphi=tempphi;
	    midpsi=temppsi;
	  }

      }

  }
 while((middelta-lowerdelta)>1e-3 || (upperdelta-middelta)>1e-3);


 *deltaminimax=middelta;
 *maxen=miden;

}




//trialproperties_seq uses the method given in Section 19.2 of Jennison and Turnbull (2000) to find the probability of stopping at each stage in a sequential trial using Z-tests

void trialproperties_seq(vector<double>& parameters,double delta0,double delta1,double sigma,double *typeIerror,double *power,double *expectedsamplesize_null,double *expectedsamplesize_crd,double *worstcasedelta,double *expectedsamplesize_dm,int checkdm)
{
  //Function will find typeIerror and power for trial parameters. If checkdm==1, the worst-case scenario delta will be found together with its expected sample size. Else, both will be returned as 0

  *worstcasedelta=0;
  *expectedsamplesize_dm=0;

  int i,j;
 
  //get grid of points to use
  vector<vector<double> > x;

  vector<double> tempvector;

  for(i=0;i<(parameters.size()-1.0)/2-1;i++)
    {
      
      
      tempvector.clear();
      tempvector.push_back(parameters.at(i*2+1));
      for(j=1;j<=15;j++)
	{
	  if(delta0*sqrt(information((i+1)*parameters.at(0),sigma))-(3+4*log(16.0/j))<parameters.at(i*2+2) && delta0*sqrt(information((i+1)*parameters.at(0),sigma))-(3.0+4*log(16.0/j))>parameters.at(i*2+1))
	    {
	      tempvector.push_back(delta0*sqrt(information((i+1)*parameters.at(0),sigma))-(3.0+4*log(16.0/j)));
	    }
	  
	}
     
      for(j=16;j<=5*16;j++)
	{
	  
	  if(delta0*sqrt(information((i+1)*parameters.at(0),sigma))-(3.0-3*(j-16.0)/(2*16.0))<parameters.at(i*2+2) && delta0*sqrt(information((i+1)*parameters.at(0),sigma))-(3.0-3*(j-16.0)/(2*16.0))>parameters.at(i*2+1))
	    {
	      tempvector.push_back(delta0*sqrt(information((i+1)*parameters.at(0),sigma))-(3.0-3*(j-16.0)/(2*16.0)));
	      
	    }
	  
	}
     
      for(j=5*16+1;j<=6*16-1;j++)
	{
	  
	  if(delta0*sqrt(information((i+1)*parameters.at(0),sigma))+(3.0+4*log(16.0/(6*16-j)))<parameters.at(i*2+2) && delta0*sqrt(information((i+1)*parameters.at(0),sigma))+(3.0+4*log(16.0/(6*16-j)))>parameters.at(i*2+1))
	    {
	      
	      tempvector.push_back(delta0*sqrt(information((i+1)*parameters.at(0),sigma))+(3.0+4*log(16.0/(6*16-j))));
	    }
	  
	}

    tempvector.push_back(parameters.at(i*2+2));
     

           x.push_back(tempvector);
	  

	}
 
  //get z's - odd numbered points are the x's, even numbered points are their midpoints 

  vector<vector<double> > z;

  for(i=0;i<x.size();i++)
    {
      tempvector.clear();
      for(j=0;j<x.at(i).size()-1;j++)
	{
	  tempvector.push_back(x.at(i).at(j));
	  tempvector.push_back((x.at(i).at(j)+x.at(i).at(j+1))/2);
	}
      tempvector.push_back(x.at(i).at(x.at(i).size()-1));
      z.push_back(tempvector);

    }

  //z's define the weights used in the integration:

  vector<vector<double> > weights;

  for(i=0;i<z.size();i++)
    {
      tempvector.clear();
      tempvector.push_back((z.at(i).at(2)-z.at(i).at(0))/6);
      for(j=2;j<=z.at(i).size()-1;j++)
	{
	  if(j%2==0)
	    {
	      tempvector.push_back(4.0*(z.at(i).at(j)-z.at(i).at(j-2))/6);
	    }
	  else if(j%2==1)
	    {
tempvector.push_back((z.at(i).at(j+1)-z.at(i).at(j-3))/6);
	    }
	}
      tempvector.push_back((z.at(i).at(z.at(i).size()-1)-z.at(i).at(z.at(i).size()-3))/6);
      weights.push_back(tempvector);
      
     
     
    }

  
  //h is a matrix which has values of h at each point

  vector<vector<double> > h;
  int k;

  for(i=0;i<z.size();i++)
    {
      tempvector.clear();
      if(i==0)
	{
	  for(j=0;j<z.at(i).size();j++)
	    {
	      tempvector.push_back(weights.at(i).at(j)*normalpdf(z.at(i).at(j)-delta0*sqrt(information((i+1)*parameters.at(0),sigma))));
	     
	    }
	}
      else
	{
	  for(j=0;j<z.at(i).size();j++)
	    {
	      tempvector.push_back(0);
	      for(k=0;k<z.at(i-1).size();k++)
		{
		  tempvector.at(j)+=h.at(i-1).at(k)*weights.at(i).at(j)*(sqrt(information((i+1)*parameters.at(0),sigma))/sqrt(information((i+1)*parameters.at(0),sigma)-information(i*parameters.at(0),sigma)))*normalpdf((z.at(i).at(j)*sqrt(information((i+1)*parameters.at(0),sigma))-z.at(i-1).at(k)*sqrt(information(i*parameters.at(0),sigma))-delta0*(information((i+1)*parameters.at(0),sigma)-information(i*parameters.at(0),sigma)))/(sqrt(information((i+1)*parameters.at(0),sigma)-information(i*parameters.at(0),sigma))));
	 
	     
		}
	    }
	}
      h.push_back(tempvector);
 
    }


  //phi gives the probability of stopping for efficacy at each stage
  vector<double> phi,psi;

 
  converthtophi(h,z,phi,parameters,delta0,delta0,sigma);
 
 
  //calculate typeIerror by summing up phi:

  *typeIerror=0;
  for(i=0;i<phi.size();i++)
    {
      *typeIerror+=phi.at(i);
    }
 

converthtopsi(h,z,psi,parameters,delta0,delta0,sigma);

 *expectedsamplesize_null=expectedsamplesize(phi,psi,parameters);

  converthtophi(h,z,phi,parameters,delta0,delta1,sigma);
  
 *power=0;
  for(i=0;i<phi.size();i++)
    {
          *power+=phi.at(i);
    }
 
converthtopsi(h,z,psi,parameters,delta0,delta1,sigma);
    

 *expectedsamplesize_crd=expectedsamplesize(phi,psi,parameters);
 double deltaminimax;

 if(checkdm==1)
   {

 finddeltaminimax_seq(h,z,phi,parameters,delta0,delta1,sigma,&deltaminimax,expectedsamplesize_dm);
*worstcasedelta=deltaminimax;
   }

     }


void getexpectedsamplesizes(vector<double>& parameters,double delta0,double delta1,double sigma,double *typeIerror,double *power,double *expectedsamplesize_null,double *expectedsamplesize_crd,double *worstcasedelta,double *expectedsamplesize_dm,int checkdm)
{
  //Function will find typeIerror and power for trial parameters. If checkdm==1, the worst-case scenario delta will be found together with its expected sample size. Else, both will be returned as 0
  *worstcasedelta=0;
  *expectedsamplesize_dm=0;

  int i,j;
 
  //get grid of points to use
  vector<vector<double> > x;

  vector<double> tempvector;

  for(i=0;i<(parameters.size()-1.0)/2-1;i++)
    {
      
      
      tempvector.clear();
      tempvector.push_back(parameters.at(i*2+1));
      for(j=1;j<=15;j++)
	{
	  if(delta0*sqrt(information((i+1)*parameters.at(0),sigma))-(3+4*log(16.0/j))<parameters.at(i*2+2) && delta0*sqrt(information((i+1)*parameters.at(0),sigma))-(3.0+4*log(16.0/j))>parameters.at(i*2+1))
	    {
	      tempvector.push_back(delta0*sqrt(information((i+1)*parameters.at(0),sigma))-(3.0+4*log(16.0/j)));
	    }
	  
	}
     
      for(j=16;j<=5*16;j++)
	{
	  
	  if(delta0*sqrt(information((i+1)*parameters.at(0),sigma))-(3.0-3*(j-16.0)/(2*16.0))<parameters.at(i*2+2) && delta0*sqrt(information((i+1)*parameters.at(0),sigma))-(3.0-3*(j-16.0)/(2*16.0))>parameters.at(i*2+1))
	    {
	      tempvector.push_back(delta0*sqrt(information((i+1)*parameters.at(0),sigma))-(3.0-3*(j-16.0)/(2*16.0)));
	      
	    }
	  
	}
     
      for(j=5*16+1;j<=6*16-1;j++)
	{
	  
	  if(delta0*sqrt(information((i+1)*parameters.at(0),sigma))+(3.0+4*log(16.0/(6*16-j)))<parameters.at(i*2+2) && delta0*sqrt(information((i+1)*parameters.at(0),sigma))+(3.0+4*log(16.0/(6*16-j)))>parameters.at(i*2+1))
	    {
	      
	      tempvector.push_back(delta0*sqrt(information((i+1)*parameters.at(0),sigma))+(3.0+4*log(16.0/(6*16-j))));
	    }
	  
	}

    tempvector.push_back(parameters.at(i*2+2));
     

           x.push_back(tempvector);
	  

	}
 
  //get z's - odd numbered points are the x's, even numbered points are their midpoints 

  vector<vector<double> > z;

  for(i=0;i<x.size();i++)
    {
      tempvector.clear();
      for(j=0;j<x.at(i).size()-1;j++)
	{
	  tempvector.push_back(x.at(i).at(j));
	  tempvector.push_back((x.at(i).at(j)+x.at(i).at(j+1))/2);
	}
      tempvector.push_back(x.at(i).at(x.at(i).size()-1));
      z.push_back(tempvector);

    }

  //z's define the weights used in the integration:

  vector<vector<double> > weights;

  for(i=0;i<z.size();i++)
    {
      tempvector.clear();
      tempvector.push_back((z.at(i).at(2)-z.at(i).at(0))/6);
      for(j=2;j<=z.at(i).size()-1;j++)
	{
	  if(j%2==0)
	    {
	      tempvector.push_back(4.0*(z.at(i).at(j)-z.at(i).at(j-2))/6);
	    }
	  else if(j%2==1)
	    {
 tempvector.push_back((z.at(i).at(j+1)-z.at(i).at(j-3))/6);
	    }
	}
      tempvector.push_back((z.at(i).at(z.at(i).size()-1)-z.at(i).at(z.at(i).size()-3))/6);
      weights.push_back(tempvector);
          
     
    }

  
  //h is a matrix which has values of h at each point

  vector<vector<double> > h;
  int k;

  for(i=0;i<z.size();i++)
    {
      tempvector.clear();
      if(i==0)
	{
	  for(j=0;j<z.at(i).size();j++)
	    {
	      tempvector.push_back(weights.at(i).at(j)*normalpdf(z.at(i).at(j)-delta0*sqrt(information((i+1)*parameters.at(0),sigma))));
	     
	    }
	}
      else
	{
	  for(j=0;j<z.at(i).size();j++)
	    {
	      tempvector.push_back(0);
	      for(k=0;k<z.at(i-1).size();k++)
		{
		  tempvector.at(j)+=h.at(i-1).at(k)*weights.at(i).at(j)*(sqrt(information((i+1)*parameters.at(0),sigma))/sqrt(information((i+1)*parameters.at(0),sigma)-information(i*parameters.at(0),sigma)))*normalpdf((z.at(i).at(j)*sqrt(information((i+1)*parameters.at(0),sigma))-z.at(i-1).at(k)*sqrt(information(i*parameters.at(0),sigma))-delta0*(information((i+1)*parameters.at(0),sigma)-information(i*parameters.at(0),sigma)))/(sqrt(information((i+1)*parameters.at(0),sigma)-information(i*parameters.at(0),sigma))));
	 
	     
		}
	    }
	}
      h.push_back(tempvector);
 
    }


  //phi gives the probability of stopping for efficacy at each stage
  vector<double> phi,psi;
   
  converthtophi(h,z,phi,parameters,delta0,delta0,sigma);
 converthtopsi(h,z,psi,parameters,delta0,delta0,sigma);

 *expectedsamplesize_null=expectedsamplesize(phi,psi,parameters);

 
}



void getallexpectedsamplesizes(vector<double>& parameters,double delta0,double delta1,double sigma,vector<double>& expectedsamplesizes)
{
  //This function returns the expected sample size at 50000 values of delta between delta0 and 2*delta1; this is used to plot the expected sample size of a design

  int i,j;
 
  //get grid of points to use
  vector<vector<double> > x;

  vector<double> tempvector;

  for(i=0;i<(parameters.size()-1.0)/2-1;i++)
    {
      
      
      tempvector.clear();
      tempvector.push_back(parameters.at(i*2+1));
      for(j=1;j<=15;j++)
	{
	  if(delta0*sqrt(information((i+1)*parameters.at(0),sigma))-(3+4*log(16.0/j))<parameters.at(i*2+2) && delta0*sqrt(information((i+1)*parameters.at(0),sigma))-(3.0+4*log(16.0/j))>parameters.at(i*2+1))
	    {
	      tempvector.push_back(delta0*sqrt(information((i+1)*parameters.at(0),sigma))-(3.0+4*log(16.0/j)));
	    }
	  
	}
     
      for(j=16;j<=5*16;j++)
	{
	  
	  if(delta0*sqrt(information((i+1)*parameters.at(0),sigma))-(3.0-3*(j-16.0)/(2*16.0))<parameters.at(i*2+2) && delta0*sqrt(information((i+1)*parameters.at(0),sigma))-(3.0-3*(j-16.0)/(2*16.0))>parameters.at(i*2+1))
	    {
	      tempvector.push_back(delta0*sqrt(information((i+1)*parameters.at(0),sigma))-(3.0-3*(j-16.0)/(2*16.0)));
	      
	    }
	  
	}
     
      for(j=5*16+1;j<=6*16-1;j++)
	{
	  
	  if(delta0*sqrt(information((i+1)*parameters.at(0),sigma))+(3.0+4*log(16.0/(6*16-j)))<parameters.at(i*2+2) && delta0*sqrt(information((i+1)*parameters.at(0),sigma))+(3.0+4*log(16.0/(6*16-j)))>parameters.at(i*2+1))
	    {
	      
	      tempvector.push_back(delta0*sqrt(information((i+1)*parameters.at(0),sigma))+(3.0+4*log(16.0/(6*16-j))));
	    }
	  
	}

    tempvector.push_back(parameters.at(i*2+2));
     

           x.push_back(tempvector);
	  

	}
 
  //get z's - odd numbered points are the x's, even numbered points are their midpoints 

  vector<vector<double> > z;

  for(i=0;i<x.size();i++)
    {
      tempvector.clear();
      for(j=0;j<x.at(i).size()-1;j++)
	{
	  tempvector.push_back(x.at(i).at(j));
	  tempvector.push_back((x.at(i).at(j)+x.at(i).at(j+1))/2);
	}
      tempvector.push_back(x.at(i).at(x.at(i).size()-1));
      z.push_back(tempvector);

    }

  //z's define the weights used in the integration:

  vector<vector<double> > weights;

  for(i=0;i<z.size();i++)
    {
      tempvector.clear();
      tempvector.push_back((z.at(i).at(2)-z.at(i).at(0))/6);
      for(j=2;j<=z.at(i).size()-1;j++)
	{
	  if(j%2==0)
	    {
	      tempvector.push_back(4.0*(z.at(i).at(j)-z.at(i).at(j-2))/6);
	    }
	  else if(j%2==1)
	    {
tempvector.push_back((z.at(i).at(j+1)-z.at(i).at(j-3))/6);
	    }
	}
      tempvector.push_back((z.at(i).at(z.at(i).size()-1)-z.at(i).at(z.at(i).size()-3))/6);
      weights.push_back(tempvector);
      
     
     
    }

  
  //h is a matrix which has values of h at each point

  vector<vector<double> > h;
  int k;

  for(i=0;i<z.size();i++)
    {
      tempvector.clear();
      if(i==0)
	{
	  for(j=0;j<z.at(i).size();j++)
	    {
	      tempvector.push_back(weights.at(i).at(j)*normalpdf(z.at(i).at(j)-delta0*sqrt(information((i+1)*parameters.at(0),sigma))));
	     
	    }
	}
      else
	{
	  for(j=0;j<z.at(i).size();j++)
	    {
	      tempvector.push_back(0);
	      for(k=0;k<z.at(i-1).size();k++)
		{
		  tempvector.at(j)+=h.at(i-1).at(k)*weights.at(i).at(j)*(sqrt(information((i+1)*parameters.at(0),sigma))/sqrt(information((i+1)*parameters.at(0),sigma)-information(i*parameters.at(0),sigma)))*normalpdf((z.at(i).at(j)*sqrt(information((i+1)*parameters.at(0),sigma))-z.at(i-1).at(k)*sqrt(information(i*parameters.at(0),sigma))-delta0*(information((i+1)*parameters.at(0),sigma)-information(i*parameters.at(0),sigma)))/(sqrt(information((i+1)*parameters.at(0),sigma)-information(i*parameters.at(0),sigma))));
	 
	     
		}
	    }
	}
      h.push_back(tempvector);
 
    }
  double delta;

  //phi gives the probability of stopping for efficacy at each stage
  vector<double> phi,psi;
  for(i=0;i<50000;i++)
    {
      phi.clear();
      psi.clear();
      delta=delta0+i*(2*delta1-delta0)/50000;
      converthtophi(h,z,phi,parameters,delta0,delta,sigma);
      converthtopsi(h,z,psi,parameters,delta0,delta,sigma);

 expectedsamplesizes.push_back(expectedsamplesize(phi,psi,parameters));
    }
 
}



