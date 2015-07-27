//multistagefunctions_v2 provides functions that allow finding the operating characteristics of a multistage sequential trial using z-statistics. If the number of stages is below a certain threshold, integration of the multivariate normal is used. If above, the method of sequential integration is used.


#include <iostream>
#include <math.h>
#include <vector>
#include <time.h>

using namespace std;



double onestagesamplesize(double difference,double sigma,double typeIerror,double typeIIerror,double R);

double information(double numberindividuals,double sigma);


void converthtophi(vector<vector<double> >& h,vector<vector<double> >& z,vector<double>& phi,vector<double>& parameters,double delta0,double newdelta,double sigma);


double expectedsamplesize(vector<double> phi,vector<double> psi,vector<double> parameters);


void converthtopsi(vector<vector<double> >& h,vector<vector<double> >& z,vector<double>& psi,vector<double>& parameters,double delta0,double newdelta,double sigma);


void finddeltaminimax_seq(vector<vector<double> >& h,vector<vector<double> >& z,vector<double>& phi,vector<double>& parameters,double delta0,double newdelta,double sigma,double *deltaminimax,double *maxen);




void trialproperties_seq(vector<double>& parameters,double delta0,double delta1,double sigma,double *typeIerror,double *power,double *expectedsamplesize_null,double *expectedsamplesize_crd,double *worstcasedelta,double *expectedsamplesize_dm,int checkdm);

void getallexpectedsamplesizes(vector<double>& parameters,double delta0,double delta1,double sigma,vector<double>& expectedsamplesizes);
