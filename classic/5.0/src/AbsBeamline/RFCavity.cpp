// ------------------------------------------------------------------------
// $RCSfile: RFCavity.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: RFCavity
//   Defines the abstract interface for an accelerating structure.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/RFCavity.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "Physics/Physics.h"
#include <iostream>
#include <fstream>
#include "fftw3.h"


// Class RFCavity
// ------------------------------------------------------------------------

RFCavity::RFCavity():
  Component(),
  dx_m(0.0),
  dy_m(0.0),
  dz_m(0.0),
  fast_m(false),
  lengthUnit_m(1.0),
  type_m(SW),
  NumCells_m(0)
{ }


RFCavity::RFCavity(const RFCavity &right):
  Component(right),
  filename_m(right.filename_m),
  scale_m(right.scale_m),
  frequency_m(right.frequency_m),
  phase_m(right.phase_m),
  type_m(right.type_m),
  NumCells_m(right.NumCells_m),
  rmin_m(right.rmin_m),
  rmax_m(right.rmax_m),
  angle_m(right.angle_m),
  pdis_m(right.pdis_m),
  gapwidth_m(right.gapwidth_m),
  phi0_m(right.phi0_m),
  fast_m(right.fast_m),
  dx_m(right.dx_m),
  dy_m(right.dy_m),
  dz_m(right.dz_m),
  lengthUnit_m(right.lengthUnit_m),
  myFieldmap(right.myFieldmap),
  RNormal_m(NULL),
  VrNormal_m(NULL),
  DvDr_m(NULL)
{}


RFCavity::RFCavity(const string &name):
  Component(name)
{}


RFCavity::~RFCavity()
{
    Fieldmap::deleteFieldmap(filename_m);
    delete[] RNormal_m;
    delete[] VrNormal_m;
    delete[] DvDr_m;
}


void RFCavity::accept(BeamlineVisitor &visitor) const
{
  visitor.visitRFCavity(*this);
}

void RFCavity::setFieldMapFN(string fn)
{
  filename_m = fn;
}

string RFCavity::getFieldMapFN() const
{
  return filename_m;
}

void RFCavity::setAmplitudem(double vPeak)
{
  scale_m = vPeak;
}

void RFCavity::setFrequencym(double freq)
{
  frequency_m = freq;
}

void RFCavity::setPhasem(double phase)
{
  phase_m = phase;
}

void RFCavity::setCavityType(string type)
{

}

string RFCavity::getCavityType() const
{

}

void RFCavity::setFast(bool fast)
{
  fast_m = fast;
}


bool RFCavity::getFast() const
{
  return fast_m;
}

void RFCavity::setMisalignment(double x, double y, double z)
{
  dx_m = x;
  dy_m = y;
  dz_m = z;
}

void RFCavity::getMisalignment(double &x, double &y, double &z) const
{
  x = dx_m;
  y = dy_m;
  z = dz_m;
}


bool RFCavity::apply(const int &i, const double &t, double E[], double B[])
{
  Vector_t Ev(0,0,0), Bv(0,0,0);
  if (apply(RefPartBunch_m->R[i],t,Ev,Bv)) return true;
      
  E[0] = Ev(0); E[1] = Ev(1); E[2] = Ev(2);
  B[0] = Bv(0); B[1] = Bv(1); B[2] = Bv(2);

  return false;

}
bool RFCavity::apply(const int &i, const double &t, Vector_t &E, Vector_t &B)
{
  const double phase = frequency_m * t + phase_m;
    
  const Vector_t tmpR(RefPartBunch_m->R[i](0),RefPartBunch_m->R[i](1),RefPartBunch_m->R[i](2) - startField_m);
  Vector_t tmpE(0.0,0.0,0.0), tmpB(0.0,0.0,0.0);

  const bool out_of_bounds = myFieldmap->getFieldstrength(tmpR,tmpE,tmpB);
  E +=  scale_m * cos( phase ) * tmpE;
  B += -scale_m * sin( phase ) * tmpB;

  return out_of_bounds;
}
  
bool RFCavity::apply( const Vector_t &R, const double &t, Vector_t &E, Vector_t &B)
{
  const double phase = frequency_m * t + phase_m;
    
  const Vector_t tmpR(R(0),R(1),R(2) - startField_m);
  Vector_t  tmpE(0.0,0.0,0.0), tmpB(0.0,0.0,0.0);

  bool out_of_bounds = myFieldmap->getFieldstrength(tmpR,tmpE,tmpB);
  E +=  scale_m * cos( phase ) * tmpE;
  B += -scale_m * sin( phase ) * tmpB;

  return out_of_bounds;
}

void RFCavity::initialise(const PartBunch *bunch, double &startField, double &endField, const double &scaleFactor)
{
  using Physics::two_pi;
 
  Inform msg("RFCavity ");

  RefPartBunch_m = bunch;

  msg << getName() << " using file " ;
  myFieldmap = Fieldmap::getFieldmap(filename_m, fast_m);
  myFieldmap->getInfo(&msg);
  if (frequency_m != myFieldmap->getFrequency())
    {
      msg << "************ WARNING ********************************************************" << endl;
      msg << " FREQUENCY IN INPUT FILE DIFFERENT THAN IN FIELD MAP; " << endl;
      msg << frequency_m << " <> " << myFieldmap->getFrequency() << "; TAKE ON THE LATTER" << endl;
      msg << "*****************************************************************************" << endl;
      frequency_m = myFieldmap->getFrequency();
    }
  double zBegin = 0.0, zEnd = 0.0, rBegin = 0.0, rEnd = 0.0;
  myFieldmap->getFieldDimensions(zBegin, zEnd, rBegin, rEnd);
  
  endField = startField + zEnd;
  startField += zBegin;
  startField_m = startField;  
}
  
// In current version ,this function read in the cavity voltage profile data from file.  
void RFCavity::initialise(const PartBunch *bunch, const double &scaleFactor)
{
  using Physics::pi;
  
  Inform msg("visitRFCavity read voltage");

  RefPartBunch_m = bunch;

  ifstream in(filename_m.c_str());
  if (!in.good()){
    msg<< "Error in Cyclotron::readFieldMap() !"<<endl;
    msg<<" Cannot open file "<< filename_m<<", please check if it really exists."<<endl;
    exit(1);
  }
  msg<<"----------------------------------------------"<<endl;
  msg<<"            READ IN VOLTAGE DATA              "<<endl;
  msg<<"        (data format: s/L, v, dV/dr)          "<<endl;
  msg<<"----------------------------------------------"<<endl;

  in>>points_num;

  if (RNormal_m!=NULL) delete[] RNormal_m;
  if (VrNormal_m!=NULL)delete[] VrNormal_m;
  if (DvDr_m!=NULL)    delete[] DvDr_m;    

  RNormal_m  = new double[points_num];
  VrNormal_m = new double[points_num];
  DvDr_m     = new double[points_num];

  for (int i=0;i<points_num;i++)
  {
    if(in.eof()){
      msg<< "Error in Cyclotron::readFieldMap() !"<<endl;
      msg<<" Not enough data in"<< filename_m<<", please check data format."<<endl;
      exit(1);
    }
    in>>RNormal_m[i]>>VrNormal_m[i]>>DvDr_m[i];

  }
  sinAngle_m = sin(angle_m/180.0*pi);
  cosAngle_m = cos(angle_m/180.0*pi);
  
  msg<<"Cavity voltage data read successfully!"<<endl;
}

void RFCavity::finalise()
{}

void RFCavity::rescaleFieldMap(const double &scaleFactor)
{
  startField_m *= scaleFactor/lengthUnit_m;
  dx_m *= scaleFactor/lengthUnit_m;
  dy_m *= scaleFactor/lengthUnit_m;
  dz_m *= scaleFactor/lengthUnit_m;
  
  myFieldmap->rescale(scaleFactor);
  lengthUnit_m = scaleFactor;
}

bool RFCavity::bends() const
{
  return false;
}

void RFCavity::goOnline()
{
  Fieldmap::readMap(filename_m);
  online_m = true;
}

void RFCavity::goOffline()
{
  Fieldmap::freeMap(filename_m);
  online_m = false;
}

void  RFCavity::setRmin(double rmin){
  rmin_m = rmin;
}
  
void  RFCavity::setRmax(double rmax){
  rmax_m = rmax;
}
  
void  RFCavity::setAzimuth(double angle){
  angle_m = angle;
}

void  RFCavity::setPerpenDistance(double pdis){
  pdis_m = pdis;
}

void  RFCavity::setGapWidth(double gapwidth){
  gapwidth_m = gapwidth;
}

void RFCavity::setPhi0(double phi0){
  phi0_m = phi0;

}

double  RFCavity::getRmin() const{
  return rmin_m;

}

double  RFCavity::getRmax() const{
  return rmax_m;

}
  
double  RFCavity::getAzimuth() const{
  return angle_m;
  
}

double  RFCavity::getSinAzimuth() const{
  return sinAngle_m;
  
}

double  RFCavity::getCosAzimuth() const{
  return cosAngle_m;
  
}

double  RFCavity::getPerpenDistance() const{
  return pdis_m;
  
}
  
double  RFCavity::getGapWidth() const{
  return gapwidth_m;
  
}

double RFCavity::getPhi0() const{

  return phi0_m;

}

void RFCavity::setComponentType(string name){
  
  if (name == "STANDING")
    type_m = SW;
  else if (name =="SINGLEGAP")
    type_m = SGSW;
  else
    {
      if (name != "")
        {
          Inform msg("RFCavity ");
          msg << "************ WARNING ********************************************************" << endl;
          msg << " CAVITY TYPE " << name << " DOES NOT EXIST; CHANGING TO REGULAR STANDING WAVE" << endl;
          msg << "*****************************************************************************" << endl;
        }
      type_m = SW;
    }
  
}

string RFCavity::getComponentType()const{

  if (type_m == SGSW)
    return string("SINGLEGAP");
  else
    return string("STANDING");
}

double RFCavity::getCycFrequency()const{

  return  frequency_m;

}
void RFCavity::getMomentaKick(double nomalRadius,double momentum[], double t,double dtCorrt)
{
  
  using Physics::m_p; //  0.93827231GeV
  using Physics::two_pi;
  using Physics::pi;
  using Physics::c;
  double derivate;
  double Voltage,dVdR;
  
  double momentum2  = momentum[0]*momentum[0]+momentum[1]*momentum[1]+momentum[2]*momentum[2];
  double betgam =sqrt(momentum2);
  
  double gamma = sqrt(1.0 +momentum2);
  double beta = betgam/gamma;

  Voltage = spline(nomalRadius, &derivate )*scale_m; // V
  // *gmsg<<" Voltage = " << Voltage/1000.0 <<"[kV]"<<endl;
  dVdR = derivate;

  if( fabs(Voltage)==0.0 ) dVdR = 0.0;

  //correct for transit time effect
  // U = 1/2*omega*deltT Vnew = V*sinU/U;
  //frequency_m = 2pi*Frf ,  [ rad/s]

  double transit_factor = 0.0;
  double U= 1.0;
  
  if (gapwidth_m > 0.0)
  {
    transit_factor= 0.5 * frequency_m * gapwidth_m*1.0e-3 / (c*beta);
    U = sin(transit_factor)/transit_factor;
  }
  
  Voltage *= U;

  double dgam=0.0, defl = 0.0;
  double nphase = ( frequency_m *(t+dtCorrt)*1.0e-9) - phi0_m/180.0*pi ; // rad/s, ns --> rad

  dgam = Voltage*cos(nphase)/(m_p*1.0e9);

  double tempdegree =fmod (nphase*360.0/two_pi,360.0);
  if (tempdegree > 270.0) tempdegree -=360.0;
  
  gamma += dgam;
  double newmomentum2 = pow(gamma,2 ) - 1.0;
  
  double pr = momentum[0]* cosAngle_m + momentum[1]* sinAngle_m;
  double ptheta = sqrt(newmomentum2 - pow(pr,2));
  momentum[0] = pr*cosAngle_m - ptheta*sinAngle_m ; // x
  momentum[1] = pr*sinAngle_m + ptheta*cosAngle_m; // y

  // *gmsg<<"After Cavity crossing, pr = "<< pr << "          [CU] "<<"ptheta = " <<ptheta<<"           [CU]"<<endl;
  // *gmsg<<"Cavity phase = " << tempdegree <<"       [deg.] , t = " <<t<<"  [ns], transit factor =  " <<U <<endl;
  // *gmsg<<"Energy Gain in Gap-crossing dE = "<< (dgam)*m_p*1000.0<<"    [MeV]"<<endl;
  // *gmsg<<"Current Energy = "<< (gamma-1.0)*m_p*1000.0<<"    [MeV]"<<endl;
  
  /*
    if (phi0_m > 0.0){
    //TODO: need to check units.
    double Bcav = -2.5e6*scale_m*Voltage*derivate*sin(nphase)*U/(gapwidth_m* frequency_m*(Rmax_m-Rmin_m)*two_pi); // T 
    defl = -Bcav*gapwidth_m*1.0e-3/betgam/(m_p*1.0e9);
  }

  pr += DEFL*pl;
  ptheta = sqrt(newmomentum2 - pow(pr,2));
  momentum[0] = pr*cosAngle_m - ptheta*sinAngle_m ; // x
  momentum[1] = pr*sinAngle_m + ptheta*cosAngle_m; // y
  *gmsg<<"After correction for cavity magnetic field: "<<endl;
  *gmsg<<" pr = " << pr <<""<<endl<<", ptheta = " << ptheta <<endl;
  */
}

/* cubic spline subrutine */ 
double RFCavity::spline( double z, double *za){
  double splint;

  // domain-test and handling of case "1-support-point"
  if(points_num < 1){
    printf("Error in RFCavity::SPLINT(): No Support-Points ! \n");
    exit(1);
  }
  if(points_num == 1){
    splint=RNormal_m[0];    
    *za = 0.0;
    return splint;
  }
  
  // search the two support-points
  int il,ih;
  il = 0;
  ih = points_num-1;
  while(( ih -il ) > 1 )
  {
    int i= (int)(( il+ih ) / 2.0);
    if     (z < RNormal_m[i]){
       ih=i;
    }
    else if(z > RNormal_m[i]){
       il=i;
    }
    else if(z == RNormal_m[i]){
       il=i;
       ih=i+1;
       break; 
    }
  }

  // Inform msg("visitRFCavity read voltage");

  //*gmsg <<"points_num = "<<points_num<<", ih = "<<ih<<", il = "<<il<<endl;
  
  double x1 =  RNormal_m[il];
  double x2 =  RNormal_m[ih];
  double y1 =  VrNormal_m[il];
  double y2 =  VrNormal_m[ih];
  double y1a = DvDr_m[il];
  double y2a = DvDr_m[ih];
  //
  // determination of the requested function-values and its derivatives
  //
  double dx =x2-x1;
  double dy = y2-y1;
  double u  = (z-x1)/dx;
  double u2 = pow(u,2);
  double u3 = pow(u,3);
  double dx2= pow(dx,2);
  double dx3= pow(dx,3);
  double dy2= -2.0*dy;
  double ya2= y2a+2.0*y1a;
  double dy3= 3.0*dy;
  double ya3=y2a+y1a;
  double yb2=dy2+dx*ya3;
  double yb4=dy3-dx*ya2;
  splint=y1+u*dx*y1a+u2*yb4+u3*yb2;
  *za=y1a+2.0*u/dx*yb4+3.0*u2/dx*yb2;
  // if(m>=1) za=y1a+2.0*u/dx*yb4+3.0*u2/dx*yb2;
  // if(m>=2) za[1]=2.0/dx2*yb4+6.0*u/dx2*yb2;
  // if(m>=3) za[2]=6.0/dx3*yb2;
 
  return splint;
}
