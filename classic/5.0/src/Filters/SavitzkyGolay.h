#ifndef CLASSIC_SAVITZKY_GOLAY_FILTER_HH
#define CLASSIC_SAVITZKY_GOLAY_FILTER_HH

#include "Filters/Filter.h"

class SavitzkyGolayFilter: public Filter {
public:
    SavitzkyGolayFilter(int np, int nl, int nr, int m);
    ~SavitzkyGolayFilter()
    { ;}

    void apply(vector<double> &histogram);
    void calc_derivative(vector<double> &histogram, const double &h);

private:
    int NumberPoints_m;
    int NumberPointsLeft_m;
    int NumberPointsRight_m;
    int PolynomialOrder_m;
    vector<double> Coefs_m;
    vector<double> CoefsDeriv_m;

};

void savgol(vector<double> &c, const int &np, const int &nl, const int &nr, const int &ld, const int &m);
void convlv(const vector<double> &data, const vector<double> &respns, const int &isign, vector<double> &ans);
void ludcmp(vector<double> &a, vector<int> &indx, double &d);
void lubksb(vector<double> &a, vector<int> &indx, vector<double> &b);

#endif // CLASSIC_SAVITZKY_GOLAY_FILTER_HH
