#ifndef OPAL_FILTERS_HH
#define OPAL_FILTERS_HH

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

class SavitzkyGolayFilter
{
 public:
  SavitzkyGolayFilter(int np, int nl, int nr, int m);
  ~SavitzkyGolayFilter()
    { ;}

  void apply(vector<double> &LineDensity);
 
 private:
  int NumberPoints;
  int NumberPointsLeft;
  int NumberPointsRight;
  int PolynomialOrder;
  vector<double> Coefs;

};

class IlyaPogorelovFilter
{
 public:
  static void apply(vector<double> &LineDensity, vector<double> &firstDerivative, const double &hz);
};

class S_G_FilterOptions
{
 public:
  S_G_FilterOptions(const int &np = 33, const int &nl = 16, const int &nr = 16, const int &m = 2):
    np_m(np),
    nl_m(nl),
    nr_m(nr),
    m_m(m)
    {  
      check();
    }

  void check()
    {
      if (nl_m < 0)
        {
          cerr << "FilterOptions: the number of points to the left has to be greater than or \n"
               << "               equal to zero; resetting nl to 16" << endl;
          nl_m = 16;
        }
      if (nr_m < 0)
        {
          cerr << "FilterOptions: the number of points to the right has to be greater than or \n"
               << "               equal to zero; resetting nr to 16" << endl;
          nr_m = 16;
        }
      if ( nl_m + nr_m < m_m)
        {
          cerr << "FilterOptions: the sum of the number of points to the left and to the  right \n"
               << "has to be greater than or equal to polynomial order; resetting nr and nl to (m+1)/2" << endl;
          nl_m = nr_m = (int)((m_m+1)/2);
          np_m = nl_m + nr_m + 1;
        }
        
      if (np_m < nl_m + nr_m + 1)
        {
          cerr << "FilterOptions: the sum of the number of points to the left and to the  right \n"
               << "has to be less than or equal to the total number of points minus 1; resetting np = nl + nr" << endl;
          np_m = nl_m + nr_m + 1;
        }
    }
  int np_m; //number of points
  int nl_m; //number of points to the left
  int nr_m; //number of points to the right
  int m_m;  //polynomial order 

};

void savgol(vector<double> &c, const int &np, const int &nl, const int &nr, const int &ld, const int &m);
void convlv(const vector<double> &data, const vector<double> &respns, const int &isign, vector<double> &ans);
void ludcmp(vector<double> &a, vector<int> &indx, double &d);
void lubksb(vector<double> &a, vector<int> &indx, vector<double> &b);


#endif //OPAL_FILTERS_HH
