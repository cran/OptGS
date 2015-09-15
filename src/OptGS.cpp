#include <vector>
#include "multistagefunctions_v3.h"
#include "normaldistribution.h"
#include "neldermead.h"
#include <Rmath.h>
#include <R.h>



using namespace std;


double functionvalue_twoparameter(double logC1,double logC2,double Deltaf,double Deltae,double requiredtypeIerror,double requiredpower,double delta,double sigma,int J,double *typeIerror,double *power)
{
 
 
  double C1=exp(logC1),C2=exp(logC2);
double expectedsamplesize_null,expectedsamplesize_crd,worstcasedelta,expectedsamplesize_dm;


int j;

  vector<double> parameters;
  double n=(2*sigma*sigma*pow(C1+C2,2)/(delta*delta))/J;
  parameters.push_back(n);
double information=pow(C1+C2,2)/(delta*delta);


              parameters.clear();
 n=(2*sigma*sigma*pow(C1+C2,2)/(delta*delta))/J;

 parameters.push_back(n);
 information=pow(C1+C2,2)/(delta*delta);
 
 for(j=0;j<J;j++)
   {
     parameters.push_back(delta*sqrt(information*(j+1.0)/J)-C2*(pow((j+1.0)/J,Deltaf-0.5)));
     parameters.push_back(C1*(pow((j+1.0)/J,Deltae-0.5)));
   }
 
 trialproperties_seq(parameters,0,delta,sigma,typeIerror,power,&expectedsamplesize_null,&expectedsamplesize_crd,&worstcasedelta,&expectedsamplesize_dm,0);

 return sqrt(pow((requiredtypeIerror-*typeIerror)/requiredtypeIerror,2)+pow((*power-requiredpower)/(1-requiredpower),2));

}


void findc1andc2_twoparameter(double Deltaf,double Deltae,double requiredtypeIerror,double requiredpower, double delta, double sigma, int J,vector<double>& finallogc1andc2,double *error)
{
  int i;

  GetRNGstate();
 
 
 finallogc1andc2.clear();
 //find c1 and c2 such that typeI error and power are correct

 vector<double> start;
 double ynewlo,reqmin=1e-06,step[2],minx,typeIerror,power;
 int konvge=10,kcount=500,icount,numres,ifault;
 start.clear();
 start.push_back(0.5);
 start.push_back(0.5);
 do
   {


     icount=0;
     numres=0;
    
 
 step[0]=1;
 step[1]=1;


 nelmin (start,finallogc1andc2,&ynewlo,reqmin,step,konvge,kcount,&icount,&numres,&ifault,Deltaf,Deltae,requiredtypeIerror,requiredpower,delta,sigma,J,&typeIerror,&power);

 finallogc1andc2.at(0)+=0.0001;
finallogc1andc2.at(1)+=0.0001;

 

 double temp=functionvalue_twoparameter(finallogc1andc2.at(0),finallogc1andc2.at(1),Deltaf,Deltae,requiredtypeIerror,requiredpower,delta,sigma,J,&typeIerror,&power);
 //if initial starting values don't result in convergence, try different ones:
 start.clear();
 start.push_back(runif(0,1));
 start.push_back(runif(0,1));
}
 while(typeIerror>requiredtypeIerror+0.001 || typeIerror<requiredtypeIerror-0.001 || power<requiredpower-0.001 || power>requiredpower+0.001);

 PutRNGstate();
}




extern "C" {

void powerfamily_twoparameter_nonintegern(double *Deltaf,double *Deltae,double *requiredtypeIerror,double *requiredpower,double *delta,double *sigma,int *J, double *finalparameters,double *expectedsamplesize,double *error,int *methodnumber,double *balanceparameter)
{
  *error=0;
 double typeIerror,power,expectedsamplesize_null,expectedsamplesize_crd,worstcasedelta,expectedsamplesize_dm=0;

 int checkdm;
 if(*methodnumber==1 || *methodnumber==2)
   {
     checkdm=0;
   }
 else
   {
     checkdm=1;
   }
 //if method number is 4, create balanced design weightings
 vector<double> balancevector;
 double balanceparameter_double;


 if(*methodnumber==4)
   {
     balancevector.push_back(balanceparameter[0]);
     balancevector.push_back(balanceparameter[1]);
     balancevector.push_back(balanceparameter[2]);
     balancevector.push_back(balanceparameter[3]);

     if(balancevector.at(2)==0){checkdm=0;}
   }
 else
   {
     balanceparameter_double=*balanceparameter;
   }
 
 
 //find C1 and C2 such that required type I error and power are correct

 vector<double> finallogc1andc2;
 findc1andc2_twoparameter(*Deltaf,*Deltae,*requiredtypeIerror,*requiredpower,*delta,*sigma,*J,finallogc1andc2,error);

 double logC1=finallogc1andc2.at(0),logC2=finallogc1andc2.at(1);
 
  int j;
  vector<double> parameters;
  double n=(2*pow(*sigma,2)*pow(exp(logC1)+exp(logC2),2)/(pow(*delta,2)))/(*J);
 double information=pow(exp(logC1)+exp(logC2),2)/(pow(*delta,2));
 
 parameters.clear();
 parameters.push_back(n);
 finalparameters[0]=n;
 finalparameters[1]=logC1;
 finalparameters[2]=logC2;
 //finalparameters 1st entry sample size per arm per stage; 2nd entry logc1, 3rd entry logc2
int correctdesign=1;

 for(j=0;j<(*J);j++)
   {
     parameters.push_back(*delta*sqrt(information*(j+1.0)/(*J))-exp(logC2)*(pow((j+1.0)/(*J),*Deltaf-0.5)));
     parameters.push_back(exp(logC1)*(pow((j+1.0)/(*J),(*Deltae)-0.5)));
  
if(parameters.at(parameters.size()-2)>parameters.at(parameters.size()-1))
       {
	 correctdesign=0;
       }

   }

 if(correctdesign==0)
   {
    
    *expectedsamplesize=9e50;
    return;
   }


 trialproperties_seq(parameters,0,*delta,*sigma,&typeIerror,&power,&expectedsamplesize_null,&expectedsamplesize_crd,&worstcasedelta,&expectedsamplesize_dm,checkdm);
 
 if(*methodnumber==1)
   {
 *expectedsamplesize=(balanceparameter_double)*expectedsamplesize_null+(1-balanceparameter_double)*(*J*parameters.at(0));
   }

 else if(*methodnumber==2)
   {
 *expectedsamplesize=balanceparameter_double*expectedsamplesize_crd+(1-balanceparameter_double)*(*J*parameters.at(0));
   }

 else if(*methodnumber==3)
   {
 *expectedsamplesize=(balanceparameter_double)*expectedsamplesize_dm+(1-balanceparameter_double)*(*J*parameters.at(0));
   }

 else
   {
  
     *expectedsamplesize=balancevector.at(0)*expectedsamplesize_null+balancevector.at(1)*expectedsamplesize_crd+balancevector.at(2)*expectedsamplesize_dm+balancevector.at(3)*(*J*parameters.at(0));
   }

}



}





extern "C" {

void powerfamily_fixedn(double *Deltaf,double *Deltae,double *logC1,double *requiredtypeIerror,double *requiredpower,double *requiredn,double *delta,double *sigma,int *J, double *finalparameters,double *functionvalue,double *penaltyfactor,int *methodnumber,double *balanceparameter)
{

 double typeIerror,power,expectedsamplesize_null,expectedsamplesize_crd,worstcasedelta,expectedsamplesize_dm=0;

  
  int j;

  vector<double> parameters;
  double n=*requiredn;


 int checkdm;
 if(*methodnumber==1 || *methodnumber==2)
   {
     checkdm=0;
   }
 else
   {
     checkdm=1;
   }



 //if method number is 4, create balanced design weightings
 vector<double> balancevector;
 double balanceparameter_double;

 if(*methodnumber==4)
   {
     balancevector.push_back(balanceparameter[0]);
 balancevector.push_back(balanceparameter[1]);
 balancevector.push_back(balanceparameter[2]);
 balancevector.push_back(balanceparameter[3]);
 if(balancevector.at(2)==0){checkdm=0;}
   }
 else
   {
     balanceparameter_double=*balanceparameter;
   }
 
 
  double logC2=log(sqrt(((*J)*pow(*delta,2)*n)/(2*pow(*sigma,2)))-exp(*logC1));

 parameters.clear();
 parameters.push_back(n);
 finalparameters[0]=n;
 finalparameters[1]=logC2;

double information=pow(exp(*logC1)+exp(logC2),2)/(pow(*delta,2));
 int correctdesign=1;

for(j=0;j<(*J);j++)
   {
     parameters.push_back(*delta*sqrt(information*(j+1.0)/(*J))-exp(logC2)*(pow((j+1.0)/(*J),*Deltaf-0.5)));
     parameters.push_back(exp(*logC1)*(pow((j+1.0)/(*J),(*Deltae)-0.5)));
     if(parameters.at(parameters.size()-2)>parameters.at(parameters.size()-1)+1e-10)
       {
	 correctdesign=0;
       }
   }
//check if any futility boundaries are above the efficacy boundares; if so, reutrn 9e50:

 if(correctdesign==0)
   {
    
    *functionvalue=9e50;
    return;
   }

  trialproperties_seq(parameters,0,*delta,*sigma,&typeIerror,&power,&expectedsamplesize_null,&expectedsamplesize_crd,&worstcasedelta,&expectedsamplesize_dm,checkdm);


 if(*methodnumber==1)
   {
 *functionvalue=(balanceparameter_double)*expectedsamplesize_null+(1-balanceparameter_double)*(*J*parameters.at(0));
   }

 else if(*methodnumber==2)
   {
 *functionvalue=balanceparameter_double*expectedsamplesize_crd+(1-balanceparameter_double)*(*J*parameters.at(0));
   }

 else if(*methodnumber==3)
   {
 *functionvalue=(balanceparameter_double)*expectedsamplesize_dm+(1-balanceparameter_double)*(*J*parameters.at(0));
   }

 else
   {
     *functionvalue=balancevector.at(0)*expectedsamplesize_null+balancevector.at(1)*expectedsamplesize_crd+balancevector.at(2)*expectedsamplesize_dm+balancevector.at(3)*(*J*parameters.at(0));
   }


*functionvalue+=(*penaltyfactor)*sqrt(pow(typeIerror-*requiredtypeIerror,2));
if(typeIerror>*requiredtypeIerror)
  {
*functionvalue+=(*penaltyfactor)*(typeIerror-*requiredtypeIerror);
  }


 *functionvalue+=(*penaltyfactor)*sqrt(pow(*requiredpower-power,2));
 
 if(power<*requiredpower)
  {
*functionvalue+=(*penaltyfactor)*(*requiredpower-power);
  }


}

}




extern "C" {

void operatingcharacteristics(double *Deltaf,double *Deltae,double *logC1,double *requiredn,double *logC2,double *delta,double *sigma,int *J,double *typeIerror,double *power, double *finalparameters,double *expectedsamplesize_null,double *expectedsamplesize_crd,double *expectedsamplesize_dm)
{

 double worstcasedelta;
 int j;

  vector<double> parameters;
  double n=*requiredn;

 *logC2=log(sqrt(((*J)*pow(*delta,2)*n)/(2*pow(*sigma,2)))-exp(*logC1));

  

 parameters.clear();
 parameters.push_back(n);
 finalparameters[0]=n;
 

double information=pow(exp(*logC1)+exp(*logC2),2)/(pow(*delta,2));

for(j=0;j<(*J);j++)
   {
     parameters.push_back(*delta*sqrt(information*(j+1.0)/(*J))-exp(*logC2)*(pow((j+1.0)/(*J),*Deltaf-0.5)));
     parameters.push_back(exp(*logC1)*(pow((j+1.0)/(*J),(*Deltae)-0.5)));
     finalparameters[j*2+1]=*delta*sqrt(information*(j+1.0)/(*J))-exp(*logC2)*(pow((j+1.0)/(*J),*Deltaf-0.5));
     finalparameters[j*2+2]=exp(*logC1)*(pow((j+1.0)/(*J),(*Deltae)-0.5));

    
   }

  trialproperties_seq(parameters,0,*delta,*sigma,typeIerror,power,expectedsamplesize_null,expectedsamplesize_crd,&worstcasedelta,expectedsamplesize_dm,1);




}




}



extern "C" {
 void extendedpowerfamily(double *Deltaf,double *Deltae,double *requiredtypeIerror,double *requiredpower,double *delta,double *sigma,int *J, double *finalparameters)
  {

 vector<double> finallogc1andc2;
 double error;
 findc1andc2_twoparameter(*Deltaf,*Deltae,*requiredtypeIerror,*requiredpower,*delta,*sigma,*J,finallogc1andc2,&error);

 double logC1=finallogc1andc2.at(0),logC2=finallogc1andc2.at(1);

 int j;
  vector<double> parameters;
  double n=(2*pow(*sigma,2)*pow(exp(logC1)+exp(logC2),2)/(pow(*delta,2)))/(*J);
 

 finalparameters[0]=n;
 finalparameters[1]=logC1;
 finalparameters[2]=logC2;


 return;
  }
}



extern "C" {

void powerfamily_fourparameter_nonintegern(double *Deltaf,double *Deltae,double *logC1,double *logC2,double *requiredtypeIerror,double *requiredpower,double *delta,double *sigma,int *J, double *finalparameters,double *functionvalue,double *penaltyfactor)
{

 double typeIerror,power,expectedsamplesize_null,expectedsamplesize_crd,worstcasedelta,expectedsamplesize_dm;

  
  int j;

  vector<double> parameters;
 
double n=(2*pow(*sigma,2)*pow(exp(*logC1)+exp(*logC2),2)/(pow(*delta,2)))/(*J);

 parameters.clear();
 parameters.push_back(n);
 finalparameters[0]=n;

double information=pow(exp(*logC1)+exp(*logC2),2)/(pow(*delta,2));
 
for(j=0;j<(*J);j++)
   {
     parameters.push_back(*delta*sqrt(information*(j+1.0)/(*J))-exp(*logC2)*(pow((j+1.0)/(*J),*Deltaf-0.5)));
     parameters.push_back(exp(*logC1)*(pow((j+1.0)/(*J),(*Deltae)-0.5)));
     finalparameters[j*2+1]=(*delta)*sqrt(information*(j+1.0)/(*J))-exp(*logC2)*(pow((j+1.0)/(*J),*Deltaf-0.5));
finalparameters[j*2+2]=(*delta)*sqrt(information*(j+1.0)/(*J))-exp(*logC2)*(pow((j+1.0)/(*J),*Deltaf-0.5));
   }
 trialproperties_seq(parameters,0,*delta,*sigma,&typeIerror,&power,&expectedsamplesize_null,&expectedsamplesize_crd,&worstcasedelta,&expectedsamplesize_dm,0);

 *functionvalue=expectedsamplesize_null;

 if(typeIerror>*requiredtypeIerror)
   {
     *functionvalue+=(*penaltyfactor)*(typeIerror-*requiredtypeIerror);
   }
 if(power<*requiredpower)
   {
     *functionvalue+=(*penaltyfactor)*(*requiredpower-power);
   }
 


}
}





extern "C" {

void getexpectedsamplesizes(double *parameters,double *delta,double *sigma,int *J, double *expectedsamplesizes)
{

 vector<double> parameters_vec;
 int i;
 parameters_vec.push_back(parameters[0]);
 for(i=1;i<(*J*2+1);i++)
   {
     parameters_vec.push_back(parameters[i]);
   }
 

 vector<double> expectedsamplesizes_vec;

 getallexpectedsamplesizes(parameters_vec,0,*delta,*sigma,expectedsamplesizes_vec);

 for(i=0;i<expectedsamplesizes_vec.size();i++)
   {
     expectedsamplesizes[i]=expectedsamplesizes_vec.at(i);
   }




}


}

