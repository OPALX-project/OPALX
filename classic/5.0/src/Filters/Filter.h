#ifndef CLASSIC_FILTER_HH
#define CLASSIC_FILTER_HH

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

class Filter
{

 public:
  Filter(){ ;}
  virtual void apply(vector<double> &histogram) = 0;
  virtual void calc_derivative(vector<double> &histogram, const double &h) = 0;

};

#endif //CLASSIC_FILTER_HH
