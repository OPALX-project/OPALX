#ifndef CLASSIC_field_HH
#define CLASSIC_field_HH
#include<string>
#include"Point.hh"

using namespace std;
//define field map class
class Field 
{
public:

  Field(
        const double Rmin,
        const double Rmax,
        const double dR,
        const double Thetamin,
        const double Thetamax,
        const double dTheta
        );
   virtual ~Field();

  bool setField(Point* V, int NumPoints, const double Bfield);
  bool saveField(const string filename);
  

private:
  
  double* bfld_m;   //Bz       unit: kGauss
  double* dbt_m;    //dBz/dtheta
  double* dbtt_m;   //d2Bz/dtheta2
  double* dbttt_m;  //d3Bz/dtheta3
  Point*  R_m;      // (x,y)   unit: mm

  double Rmin_m;    // unit: mm
  double Rmax_m;    // unit: mm
  double dR_m;      // unit: mm

  double Thetamin_m;// unit: degree
  double Thetamax_m;// unit: degree
  double dTheta_m;  // unit: degree

  int NpTh_m;  // don't include the points of the arc
  int NpR_m; //  include the outer point of the Radius 
  long Ntotal_m; // total points 
};

#endif
