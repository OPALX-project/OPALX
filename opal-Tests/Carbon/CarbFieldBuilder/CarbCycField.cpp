#include"CarbCycField.hh"
#include"CheckPoint.hh"
#include<iostream>
#include<cmath>
#include<fstream>
using namespace checkPoint;
using namespace std;


const double pi         = 3.14159265358979323846;
const double two_pi     = 2.0 * pi;
const double u_two_pi   = 1.0 / (2.0 * pi);

Field::Field(const double Rmin,const double Rmax,const double dR,const double Thetamin,const double Thetamax,const double dTheta):
  Rmin_m(Rmin),Rmax_m(Rmax),dR_m(dR), Thetamin_m(Thetamin), Thetamax_m(Thetamax) 
{
  if (dTheta < 0.0)
    dTheta_m = -1.0/dTheta;
  else dTheta_m = dTheta;

  if (dR < 0.0)
    dR_m = -1.0/dR;
  else dR_m = dR;
 
  NpTh_m =(int)((Thetamax_m-Thetamin_m)/dTheta_m); 
  NpR_m  = (int)((Rmax_m-Rmin_m)/dR_m + 1);
  Ntotal_m = (long)(NpR_m * NpTh_m);

  bfld_m =  new double[Ntotal_m];
  dbt_m =   new double[Ntotal_m];
  dbtt_m =  new double[Ntotal_m];
  dbttt_m = new double[Ntotal_m];
  R_m =     new Point [Ntotal_m];

  for (int ii = 0; ii < Ntotal_m; ii++){
    bfld_m[ii] = 0.0;
    dbt_m[ii] = 0.0;
    dbtt_m[ii] = 0.0;
    dbttt_m[ii] = 0.0;
  }
  
    
}

bool Field::setField(Point* V, int NumPoints, const double Bfield)
{
  ofstream outf;
  outf.setf(ios::scientific, ios::floatfield);
  outf.precision(8);
  outf.open("Field.data", ios::app);

  for(int ii = 0; ii < NpR_m; ii++){

    double r = Rmin_m + ii*dR_m;
  
    for (int jj = 0; jj < NpTh_m; jj++)
    {
      int index = ii*NpTh_m + jj;
          
      double theta = (Thetamin_m + jj*dTheta_m);
     
      R_m[index].x = r*cos(theta / 180.0 * pi);
      R_m[index].y = r*sin(theta / 180.0 * pi);
      
      int flag = cn_PnPoly( R_m[index], V, NumPoints );

      if (flag != 0)
      {
        //cout<<"R = ("<< r <<" , "<<theta<<" ) Inside !"<<endl;
        outf<< R_m[index].x<<"  "<<R_m[index].y<<endl;

        bfld_m[index] = Bfield;  // kGauss
        dbt_m[index]  = 0.0;
        dbtt_m[index] = 0.0;
        dbttt_m[index]= 0.0;
      }
    }
  }
  outf.close();
  return true;
}

bool Field::saveField(const string filename)
{
 
  ofstream outf;
  outf.setf(ios::scientific, ios::floatfield );
  outf.precision(8);
  outf.open(filename.c_str());

  if (dTheta_m < 1.0) dTheta_m = -abs(1/dTheta_m);
  if (dR_m < 1.0) dR_m = -abs(1/dR_m);
  outf<<Rmin_m<<endl;
  outf<<dR_m<<endl;
  outf<<Thetamin_m<<endl;
  outf<<dTheta_m<<endl;
  outf<<NpTh_m<<endl;
  outf<<NpR_m<<endl; 
  for (int ii=0; ii<Ntotal_m; ii++)
  {
    outf<<" "<<bfld_m[ii]<<" ";
    if ((ii+1)%5 == 0) outf<<endl;
  }

  outf.close();
  
  return true;
}

Field::~Field()
{

  delete[] bfld_m;
  delete[] dbt_m;
  delete[] dbtt_m;
  delete[] dbttt_m;
  delete[] R_m;

  
}
