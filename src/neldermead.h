# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <ctime>
# include <cmath>
# include <vector>

using namespace std;

void nelmin (vector<double> start,vector<double>& xmin, 
  double *ynewlo, double reqmin, double step[], int konvge, int kcount, 
	     int *icount, int *numres, int *ifault,double Deltaf,double Deltae,double requiredtypeIerror,double requiredpower,double delta,double sigma,int J,double *typeIerror,double *power);

