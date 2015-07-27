//This program takes in a covariance matrix, and two limits of integration, then uses Genz' method to find the probability of the multivariance normal with mean 0 and covariance matrix as specified being within the two limits

#include <iostream>
#include <math.h>
#include <vector>

using namespace std;

double normalpdf(double z);

double normalcdf(double z);

double inversenormalcdf(double p);
