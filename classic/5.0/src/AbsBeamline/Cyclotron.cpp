// ------------------------------------------------------------------------
// $RCSfile: Cyclotron.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Definitions for class: Cyclotron
//   Defines the abstract interface for a sector bend magnet.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Cyclotron.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "Physics/Physics.h"
#include "assert.h"

using Physics::pi;

// Class Cyclotron
// ------------------------------------------------------------------------

Cyclotron::Cyclotron():
  Component()
{}


Cyclotron::Cyclotron(const Cyclotron &right):
  Component(right),
  rinit_m(right.rinit_m),
  prinit_m(right.prinit_m),
  phiinit_m(right.phiinit_m),
  fmapfn_m(right.fmapfn_m),
  rffrequ_m(right.rffrequ_m),
  symmetry_m(right.symmetry_m),
  type_m(right.type_m),
  harm_m(right.harm_m)
{}



Cyclotron::Cyclotron(const string &name):
  Component(name)
{}


Cyclotron::~Cyclotron()
{
  delete[] Bfield.bfld;
  delete[] Bfield.dbt;
  delete[] Bfield.dbtt;
  delete[] Bfield.dbttt;

  delete[] BP.rarr;

  delete[] Bfield.dbr;
  delete[] Bfield.dbrr;
  delete[] Bfield.dbrrr;

  delete[] Bfield.dbrt;
  delete[] Bfield.dbrrt;
  delete[] Bfield.dbrtt;

  delete[] Bfield.f2;
  delete[] Bfield.f3;
  delete[] Bfield.g3;

}


void Cyclotron::accept(BeamlineVisitor &visitor) const
{
  visitor.visitCyclotron(*this);
}

void Cyclotron::setRinit(double rinit)
{
  rinit_m=rinit;
}

double Cyclotron::getRinit() 
{  
  return rinit_m;
}


void Cyclotron::setPRinit(double prinit)
{
  prinit_m=prinit;
}

double Cyclotron::getPRinit() 
{  
  return prinit_m;
}

void Cyclotron::setPHIinit(double phiinit)
{
  phiinit_m=phiinit;
}

double Cyclotron::getPHIinit() 
{  
  return phiinit_m;
}

void Cyclotron::setFieldMapFN(string f)
{
  fmapfn_m=f;
}

string Cyclotron::getFieldMapFN()
{  
  return fmapfn_m;
}

void Cyclotron::setRfFrequ(double f)
{
  rffrequ_m=f;
}

double Cyclotron::getRfFrequ() 
{  
  return rffrequ_m;
}

void Cyclotron::setSymmetry(double s)
{
  symmetry_m=s;
}

double Cyclotron::getSymmetry() 
{  
  return symmetry_m;
}


void Cyclotron::setType(string t)
{
  type_m=t;
}

string Cyclotron::getType()
{  
  return type_m;
}

void Cyclotron::setCyclHarm(double h)
{
  harm_m=h;
}

double Cyclotron::getCyclHarm() 
{  
  return harm_m;
}

double Cyclotron::getRmin() 
{  
  return BP.rmin;
}


double Cyclotron::getRmax() 
{  
  return BP.rmin + (Bfield.nrad-1)*BP.delr;
}

// This function aims at obtaining magentic field at any given location R by interpolation.
// arguments t is useless here.
// arguments E is set to zero. 

bool Cyclotron::apply(const int &i, const double &t, double E[], double B[])
{
  Vector_t Ev(0,0,0), Bv(0,0,0);
  if (apply(RefPartBunch_m->R[i],t,Ev,Bv)) return true;
      
  E[0] = Ev(0); E[1] = Ev(1); E[2] = Ev(2);
  B[0] = Bv(0); B[1] = Bv(1); B[2] = Bv(2);
      
  return false;
}

bool Cyclotron::apply(const int &i, const double &t, Vector_t &E, Vector_t &B)
{
  return apply(RefPartBunch_m->R[i],t,E,B);
}
  
bool Cyclotron::apply(const Vector_t &R, const double &t, Vector_t &E, Vector_t &B)
{
  const double rad = sqrt(R[0]*R[0] + R[1]*R[1]);
  const double xir = (rad-BP.rmin)/BP.delr;

  // ir : the mumber of path whoes radius is less then the 4 points of cell which surrond particle.
  // that is Rmin = 1900, dR = 20 , r = 1911, then ir = 0     
  const int    ir = (int)xir;

  // wr1 : the relative distance to the inner path radius 
  const double wr1 = xir - (double)ir;
  // wr2 : the relative distance to the outer path radius
  const double wr2 = 1.0 - wr1;

  double tet,tet_map, xit;
  int it;

  if      ((R[0]>0) && (R[1]>=0)) tet = atan(R[1]/R[0]);
  else if ((R[0]<0) && (R[1]>=0)) tet = pi+atan(R[1]/R[0]);
  else if ((R[0]<0) && (R[1]<=0)) tet = pi+atan(R[1]/R[0]); 
  else if ((R[0]>0) && (R[1]<=0)) tet = 2.0*pi+atan(R[1]/R[0]);
  else if ((R[0]==0) && (R[1]> 0)) tet = pi/2.0;
  else if ((R[0]==0) && (R[1]< 0)) tet = 3.0/2.0*pi;

  double tet_rad = tet;
  
  // the actual angle of particle
  tet = tet/pi*180.0;
 
  // the corresponding angle on the field map
  // Note: this does not work if the start point of field map does not equal zero.
  tet_map = fmod(tet, 360.0/symmetry_m);

  xit = tet_map/BP.dtet;

  it = (int) xit;
  const double wt1 = xit - (double)it;
  const double wt2 = 1.0 - wt1;

  // it : the number of point on the inner path whoes angle is less then the particle' corresponding angle. 
  // include zero degree point   
  it = it + 1; 
  
  int r1t1, r2t1, r1t2, r2t2;
  int ntetS = Bfield.ntet+1;

  // r1t1 : the index of the "min angle, min radius" point in the 2D field array. 
  // considering  the array start with index of zero, minus 1.  
  r1t1 = it + ntetS*ir - 1; 
  // the index of other three points
  r1t2 = r1t1 + 1;
  r2t1 = r1t1 + ntetS;
  r2t2 = r2t1 + 1 ;

  /* debug
  *gmsg << "x,y,z  = ( "<< R[0] << " , "<< R[1] << " , "<< R[2] << ")[mm] **********, "<<endl;
  *gmsg << "r2t2, r2t1, r1t2, r1t1 = (" << r2t2 <<" , "<<r2t1 <<" , "<< r1t2<<" , "<<r1t1 <<" ) " <<endl;
  *gmsg << "wr1, wt1= (" << wr1 <<" , "<<wt1 <<" ) " <<endl;
  */

  double bzcub=0.0, bzf=0.0, bz=0.0;
  double brcub=0.0, brf=0.0, br=0.0;
  double btcub=0.0, btf=0.0, bt=0.0;

  if((it >= 0) && (ir >= 0) && (it < Bfield.ntetS) && (ir < Bfield.nrad)){

    /* Bz */
    bzf =  (Bfield.bfld[r1t1] * wr2*wt2 + Bfield.bfld[r2t1] * wr1*wt2 +
           Bfield.bfld[r1t2] * wr2*wt1 + Bfield.bfld[r2t2] * wr1*wt1 );
    
    bzcub =(Bfield.f2[r1t1] * wr2*wt2 +
            Bfield.f2[r2t1] * wr1*wt2 +
            Bfield.f2[r1t2] * wr2*wt1 +
            Bfield.f2[r2t2] * wr1*wt1 ) * pow(R[2],2.0);
    
    // bz = -( bzf - bzcub );
    bz = - bzf ;
    

     /* Br */
    brf =  (Bfield.dbr[r1t1] * wr2*wt2 +
            Bfield.dbr[r2t1] * wr1*wt2 +
            Bfield.dbr[r1t2] * wr2*wt1 +
            Bfield.dbr[r2t2] * wr1*wt1 ) * R[2];


    brcub= (Bfield.f3[r1t1] * wr2*wt2 +
            Bfield.f3[r2t1] * wr1*wt2 +
            Bfield.f3[r1t2] * wr2*wt1 +
            Bfield.f3[r2t2] * wr1*wt1 ) * pow(R[2],3.0);

    // br = -( brf - brcub );
    br = - brf;
    

     /* Btheta */
    btf = ( Bfield.dbt[r1t1] * wr2*wt2 +
            Bfield.dbt[r2t1] * wr1*wt2 +
            Bfield.dbt[r1t2] * wr2*wt1 +
            Bfield.dbt[r2t2] * wr1*wt1 )/rad * R[2];


    btcub = (Bfield.g3[r1t1] * wr2*wt2 +
             Bfield.g3[r2t1] * wr1*wt2 +
             Bfield.g3[r1t2] * wr2*wt1 +
             Bfield.g3[r2t2] * wr1*wt1 )/rad * pow(R[2],3.0);

    // bt = -( btf - btcub );
    bt = - btf;
    
    /* Br Btheta -> Bx By */
    
    B[0] = br*cos(tet_rad) - bt*sin(tet_rad);
    B[1] = br*sin(tet_rad) + bt*cos(tet_rad);
    B[2] = bz;

    /* Test for homo field.
      B(0) =0.0;
      B(1) =0.0;
      B(2) = -5.00000;// 5000Gauss
    */
   
  }
  /* error output, out of field! */
  else{
    *gmsg <<"Error!" << getName() << ".getFieldstrength(): out of boudaries (z,r) = (" << R[0] << "," << R[1] << "," << R[2] << ")" << endl
          <<" rad="<<rad<<" [mm], theta="<<tet<<" [deg], it="<<it<<", ir="<<ir<<endl;
   return true;

  }
  return false;
}

void Cyclotron::initialise(const PartBunch *bunch, double &startField, double &endField, const double &scaleFactor)
{
  RefPartBunch_m = bunch;
  online_m = true;
}

void Cyclotron::finalise()
{
  online_m = false;
}

void Cyclotron::rescaleFieldMap(const double &scaleFactor)
{}

bool Cyclotron::bends() const
{
  return true;
}





// calculate derivatives with 5-point lagrange's formula.
double Cyclotron::gutdf5d(double *f, double dx, const int kor, const int krl, const int lpr)

{
  double C[5][5][3], FAC[3];
  double result;
  int j;
  /* CALCULATE DERIVATIVES WITH 5-POINT LAGRANGE FORMULA
   * PARAMETERS:
   * F  STARTADDRESS FOR THE 5 SUPPORT POINTS
   * DX STEPWIDTH FOR ARGUMENT
   * KOR        ORDER OF DERIVATIVE (KOR=1,2,3).
   * KRL        NUMBER OF SUPPORT POINT, WHERE THE DERIVATIVE IS TO BE CALCULATED
   *  (USUALLY 3, USE FOR BOUNDARY 1 ,2, RESP. 4, 5)
   * LPR        DISTANCE OF THE 5 STORAGE POSITIONS (=1 IF THEY ARE NEIGHBORS OR LENGTH
   * OF COLUMNLENGTH OF A MATRIX, IF THE SUPPORT POINTS ARE ON A LINE).
   * ATTENTION! THE INDICES ARE NOW IN C-FORMAT AND NOT IN FORTRAN-FORMAT.*/      

  /* COEFFICIENTS FOR THE 1ST DERIVATIVE: */
  C[0][0][0]=-50.0; C[1][0][0]= 96.0; C[2][0][0]=-72.0; C[3][0][0]= 32.0; C[4][0][0]= -6.0;
  C[0][1][0]= -6.0; C[1][1][0]=-20.0; C[2][1][0]= 36.0; C[3][1][0]=-12.0; C[4][1][0]=  2.0;
  C[0][2][0]=  2.0; C[1][2][0]=-16.0; C[2][2][0]=  0.0; C[3][2][0]= 16.0; C[4][2][0]= -2.0;
  C[0][3][0]= -2.0; C[1][3][0]= 12.0; C[2][3][0]=-36.0; C[3][3][0]= 20.0; C[4][3][0]=  6.0;
  C[0][4][0]=  6.0; C[1][4][0]=-32.0; C[2][4][0]= 72.0; C[3][4][0]=-96.0; C[4][4][0]= 50.0;

  /* COEFFICIENTS FOR THE 2ND DERIVATIVE: */
  C[0][0][1]= 35.0; C[1][0][1]=-104;  C[2][0][1]=114.0; C[3][0][1]=-56.0; C[4][0][1]= 11.0;
  C[0][1][1]= 11.0; C[1][1][1]=-20.0; C[2][1][1]=  6.0; C[3][1][1]=  4.0; C[4][1][1]= -1.0;
  C[0][2][1]= -1.0; C[1][2][1]= 16.0; C[2][2][1]=-30.0; C[3][2][1]= 16.0; C[4][2][1]= -1.0;
  C[0][3][1]= -1.0; C[1][3][1]=  4.0; C[2][3][1]=  6.0; C[3][3][1]=-20.0; C[4][3][1]= 11.0;
  C[0][4][1]= 11.0; C[1][4][1]=-56.0; C[2][4][1]=114.0; C[3][4][1]=-104;  C[4][4][1]= 35.0;


  /* COEFFICIENTS FOR THE 3RD DERIVATIVE: */ 
  C[0][0][2]=-10.0; C[1][0][2]= 36.0; C[2][0][2]=-48.0; C[3][0][2]= 28.0; C[4][0][2]= -6.0;
  C[0][1][2]= -6.0; C[1][1][2]= 20.0; C[2][1][2]=-24.0; C[3][1][2]= 12.0; C[4][1][2]= -2.0;
  C[0][2][2]= -2.0; C[1][2][2]=  4.0; C[2][2][2]=  0.0; C[3][2][2]= -4.0; C[4][2][2]=  2.0;
  C[0][3][2]=  2.0; C[1][3][2]=-12.0; C[2][3][2]= 24.0; C[3][3][2]=-20.0; C[4][3][2]=  6.0;
  C[0][4][2]=  6.0; C[1][4][2]=-28.0; C[2][4][2]= 48.0; C[3][4][2]=-36.0; C[4][4][2]= 10.0;

  /* FACTOR: */
  FAC[0] = 24.0;    FAC[1] = 12.0;    FAC[2] = 4.0;
  
  result = 0.0;
  for(j=0;j<5;j++){
    result += C[j][krl][kor] * *(f +j*lpr);
  }

  return result/(FAC[kor]*pow(dx,(kor+1)));
}


// evaulate other derivative of magnetic field.
void Cyclotron::getdiffs(){

  assert(Bfield.dbr   = new double[Bfield.ntot]);
  assert(Bfield.dbrr  = new double[Bfield.ntot]);
  assert(Bfield.dbrrr = new double[Bfield.ntot]);

  assert(Bfield.dbrt  = new double[Bfield.ntot]);
  assert(Bfield.dbrrt = new double[Bfield.ntot]);
  assert(Bfield.dbrtt = new double[Bfield.ntot]);
  
  assert(Bfield.f2    = new double[Bfield.ntot]);
  assert(Bfield.f3    = new double[Bfield.ntot]);
  assert(Bfield.g3    = new double[Bfield.ntot]);

  for (int i =0; i< Bfield.nrad; i++){
    
    for(int k=0; k < Bfield.ntet; k++){

      double dtheta = pi/180.0*BP.dtet;

      int kEdge;
      
      kEdge = max(k-2,0);
      kEdge = min(kEdge, Bfield.ntet-5);

      int dkFromEdge = k-kEdge;
      int index = idx(i,k);
      int indexkEdge =idx(i,kEdge);

      
      Bfield.dbt[index]    = gutdf5d(&Bfield.bfld[indexkEdge], dtheta, 0, dkFromEdge, 1 );
      Bfield.dbtt[index]   = gutdf5d(&Bfield.bfld[indexkEdge], dtheta, 1, dkFromEdge, 1 );
      Bfield.dbttt[index]  = gutdf5d(&Bfield.bfld[indexkEdge], dtheta, 2, dkFromEdge, 1 );
    }
  }
  
  

  for(int k=0; k < Bfield.ntet; k++){
     // inner loop varies R
     for(int i=0; i < Bfield.nrad; i++){
       double rac = BP.rarr[i];
       // define iredg, the reference index for radial interpolation
       // standard: i-2 minimal: 0 (not negative!)  maximal: nrad-4
       int iredg = max(i-2,0);
       iredg = min(iredg, Bfield.nrad-5);
       int irtak = i-iredg;
       int index = idx(i,k);
       int indexredg = idx(iredg,k);


       Bfield.dbr[index]    = gutdf5d(&Bfield.bfld[indexredg], BP.delr, 0, irtak, Bfield.ntetS );
       Bfield.dbrr[index]   = gutdf5d(&Bfield.bfld[indexredg], BP.delr, 1, irtak, Bfield.ntetS );
       Bfield.dbrrr[index]  = gutdf5d(&Bfield.bfld[indexredg], BP.delr, 2, irtak, Bfield.ntetS );

       Bfield.dbrt[index]   = gutdf5d(&Bfield.dbt[indexredg], BP.delr, 0, irtak, Bfield.ntetS );
       Bfield.dbrrt[index]  = gutdf5d(&Bfield.dbt[indexredg], BP.delr, 1, irtak, Bfield.ntetS );
       Bfield.dbrtt[index]  = gutdf5d(&Bfield.dbtt[indexredg], BP.delr, 0, irtak, Bfield.ntetS );

       // fehlt noch!! f2,f3,g3, 
       Bfield.f2[index] = ( Bfield.dbrr[index]
                            + Bfield.dbr[index]/rac
                            + Bfield.dbtt[index]/rac/rac )/2.0;
       
       Bfield.f3[index] = ( Bfield.dbrrr[index]
                            + Bfield.dbrr[index]/rac
                            + (Bfield.dbrtt[index] - Bfield.dbr[index])/rac/rac
                            - 2.0*Bfield.dbtt[index]/rac/rac/rac )/6.0;

       Bfield.g3[index] = (  Bfield.dbrrt[index]
                           + Bfield.dbrt[index]/rac 
                           + Bfield.dbttt[index]/rac/rac )/6.0;
     } // Radius Loop
  } // Azimuth loop

  // copy 1st azimuth to last + 1 to always yield an interval
  for(int i = 0; i< Bfield.nrad; i++){
    int iend = idx(i,Bfield.ntet);
    int istart = idx(i,0);

    Bfield.bfld[iend]   = Bfield.bfld[istart];
    Bfield.dbt[iend]    = Bfield.dbt[istart];
    Bfield.dbtt[iend]   = Bfield.dbtt[istart];
    Bfield.dbttt[iend]  = Bfield.dbttt[istart];

    Bfield.dbr[iend]    = Bfield.dbr[istart];
    Bfield.dbrr[iend]   = Bfield.dbrr[istart];
    Bfield.dbrrr[iend]  = Bfield.dbrrr[istart];

    Bfield.dbrt[iend]   = Bfield.dbrt[istart];
    Bfield.dbrtt[iend]  = Bfield.dbrtt[istart];
    Bfield.dbrrt[iend]  = Bfield.dbrrt[istart];

    Bfield.f2[iend]     = Bfield.f2[istart];
    Bfield.f3[iend]     = Bfield.f3[istart];
    Bfield.g3[iend]     = Bfield.g3[istart];
    
  }
}

// read field map from external file.
void Cyclotron::getFieldFromFile( const double &scaleFactor)
{	

  FILE *f=NULL;
  int lpar;
  char fout[100];
  double dtmp;
  
  Inform msg("visitCyclotron read field ");
  msg<<"----------------------------------------------"<<endl;
  msg<<"            READ IN RING FIELD MAP            "<<endl;
  msg<< "     (The first data block is useless)       "<<endl;
  msg<<"----------------------------------------------"<<endl;

  //  Bfield.filename = "s03av.nar"
  BP.Bfact = scaleFactor;

 if( (f=fopen(fmapfn_m.c_str(),"r")) == NULL )
 {
   msg<< "Error in Cyclotron::getFieldFromFile()!"<<endl;
   msg<<" Cannot open file, please check if it really exists."<<endl;
   exit(1);
}
  
  assert(fscanf(f, "%lf", &BP.rmin));
  msg << "Minimal radius of measured field map: "<<BP.rmin<<" [mm]"<<endl;

  assert(fscanf(f, "%lf", &BP.delr));
  msg << "Stepsize in radial direction: "<<BP.delr<<" [mm]"<<endl;
  
  assert(fscanf(f, "%lf", &BP.tetmin));
  msg << "Minimal angle of measured field map: "<<BP.tetmin<<" [deg.]"<<endl;

  assert(fscanf(f, "%lf", &BP.dtet));
  //if the value is nagtive, the actual value is its reciprocal. 
  if (BP.dtet<0.0) BP.dtet = 1.0/(-BP.dtet);
  msg << "Stepsize in azimuth direction: "<<BP.dtet<<" [deg.]"<<endl;

  for(int i=0; i<13; i++)assert(fscanf(f, "%s",fout));

  assert(fscanf(f, "%d", &Bfield.nrad));
  msg << "Index in radial direction: "<<Bfield.nrad<<endl;

  assert(fscanf(f, "%d", &Bfield.ntet));
  msg << "Index in azimuthal direction: "<<Bfield.ntet<<endl;

  Bfield.ntetS = Bfield.ntet+1;
  msg << "Accordingly, total grid point along azimuth:  "<<Bfield.ntetS<<endl;

  for(int i=0; i<5; i++){
    assert(fscanf(f, "%s",fout));
  }
  assert(fscanf(f, "%d", &lpar));
  // msg<< "READ"<<lpar<<" DATA ENTRIES"<<endl;

  for(int i=0; i<4; i++){
    assert(fscanf(f, "%s",fout));
  }

  for(int i=0; i<lpar; i++){
    assert(fscanf(f, "%16lE", &dtmp));
  }
  for(int i=0; i<6; i++){
    assert(fscanf(f, "%s",fout));
  }
  //msg << "READ FILE DESCRIPTION..." <<endl;
  for(int i=0; i<10000; i++){
    assert(fscanf(f, "%s",fout));
    if(strcmp(fout,"LREC=")==0)break;
  }

  for(int i=0; i<5; i++){
    assert(fscanf(f, "%s",fout));
  }
  Bfield.ntot = idx(Bfield.nrad-1,Bfield.ntet)+1;
  //jjyang
  msg << "Total stored grid point number ( ntetS * nrad ) : "<<Bfield.ntot<<endl;
  assert(Bfield.bfld  = new double[Bfield.ntot]);
  assert(Bfield.dbt   = new double[Bfield.ntot]);
  assert(Bfield.dbtt  = new double[Bfield.ntot]);
  assert(Bfield.dbttt = new double[Bfield.ntot]);

  msg << "read -in loop one block per radius"<< endl; 
  msg << "rescaling of the fields with factor: "<< BP.Bfact <<endl;
  for(int i=0; i < Bfield.nrad; i++){

    if(i>0){
      for(int dummy=0; dummy<6; dummy++){
        assert(fscanf(f, "%s",fout)); // INFO-LINE
      }}
    for(int k=0; k < Bfield.ntet; k++){
      assert(fscanf(f, "%16lE", &(Bfield.bfld[idx(i,k)])));
      Bfield.bfld[idx(i,k)] *= BP.Bfact;
    }
    for(int k=0; k < Bfield.ntet; k++){
      assert(fscanf(f, "%16lE", &(Bfield.dbt[idx(i,k)])));
      Bfield.dbt[idx(i,k)] *= BP.Bfact;
    }
    for(int k=0; k < Bfield.ntet; k++){
      assert(fscanf(f, "%16lE", &(Bfield.dbtt[idx(i,k)])));
      Bfield.dbtt[idx(i,k)] *= BP.Bfact;
    }
    for(int k=0; k < Bfield.ntet; k++){
      assert(fscanf(f, "%16lE", &(Bfield.dbttt[idx(i,k)])));
      Bfield.dbttt[idx(i,k)] *= BP.Bfact;
    }
  }
  fclose(f);


  msg<<"Field Map read successfully!"<<endl;
}

// Calculates Radiae of initial grid. 
// dimensions in [mm]!
void Cyclotron::initR(double rmin, double dr, int nrad)
{
  assert(BP.rarr = new double[nrad]);
  for(int i=0; i< nrad; i++){ 
    BP.rarr[i] = rmin + i*dr;
  }
  BP.delr = dr;
}


void Cyclotron::initialise(const PartBunch *bunch, const int &fieldflag, const double &scaleFactor)
{
  RefPartBunch_m = bunch;

  // for your own format field, you should add your own getFieldFromFile() function by yourself. 
  
  if(fieldflag == 1){
    //*gmsg<<"Read field data from PSI format field map file."<<endl;
    getFieldFromFile(scaleFactor);
  }else if (fieldflag == 2)
  {
    // *gmsg<<"Read data from 450MeV Carbon cyclotron field file"<<endl;
    getFieldFromFile_Carbon(scaleFactor);
  }else
    ERRORMSG("The field reading function of this TYPE of CYCLOTRON has not implemented yet.!" << endl);

  // calculate the radii of initial grid.
  initR(BP.rmin, BP.delr, Bfield.nrad); 

  // calculate the remaining derivatives
  getdiffs();

  //*gmsg<<"----------------------------------------------"<<endl;
  
}

// read field map from external file.
void Cyclotron::getFieldFromFile_Carbon( const double &scaleFactor)
{	

  FILE *f=NULL;
  int lpar;
  char fout[100];
  double dtmp;
  
  Inform msg("visitCyclotron read field ");
  msg<<"----------------------------------------------"<<endl;
  msg<<"     READ IN CARBON CYCLOTRON FIELD MAP       "<<endl;
  msg<<"----------------------------------------------"<<endl;

  BP.Bfact = scaleFactor;

 if( (f=fopen(fmapfn_m.c_str(),"r")) == NULL )
 {
   msg<< "Error in Cyclotron::getFieldFromFile_Carbon()!"<<endl;
   msg<<" Cannot open file, please check if it really exists."<<endl;
   exit(1);
}
  
  assert(fscanf(f, "%lf", &BP.rmin));
  msg << "Minimal radius of measured field map: "<<BP.rmin<<" [mm]"<<endl;

  assert(fscanf(f, "%lf", &BP.delr));
  msg << "Stepsize in radial direction: "<<BP.delr<<" [mm]"<<endl;
  
  assert(fscanf(f, "%lf", &BP.tetmin));
  msg << "Minimal angle of measured field map: "<<BP.tetmin<<" [deg.]"<<endl;

  assert(fscanf(f, "%lf", &BP.dtet));
  //if the value is nagtive, the actual value is its reciprocal. 
  if (BP.dtet<0.0) BP.dtet = 1.0/(-BP.dtet);
  msg << "Stepsize in azimuth direction: "<<BP.dtet<<" [deg.]"<<endl;

  assert(fscanf(f, "%d", &Bfield.ntet));
  msg << "Index in azimuthal direction: "<<Bfield.ntet<<endl;
  
  assert(fscanf(f, "%d", &Bfield.nrad));
  msg << "Index in radial direction: "<<Bfield.nrad<<endl;
  
  Bfield.ntetS = Bfield.ntet+1;
  msg << "Accordingly, total grid point along azimuth:  "<<Bfield.ntetS<<endl;

  Bfield.ntot = idx(Bfield.nrad-1,Bfield.ntet)+1;

  msg << "Total stored grid point number ( ntetS * nrad ) : "<<Bfield.ntot<<endl;
  assert(Bfield.bfld  = new double[Bfield.ntot]);
  assert(Bfield.dbt   = new double[Bfield.ntot]);
  assert(Bfield.dbtt  = new double[Bfield.ntot]);
  assert(Bfield.dbttt = new double[Bfield.ntot]);

  msg << "rescaling of the fields with factor: "<< BP.Bfact <<endl;
  for(int i=0; i < Bfield.nrad; i++){
    for(int k=0; k < Bfield.ntet; k++){
      assert(fscanf(f, "%16lE", &(Bfield.bfld[idx(i,k)])));
      Bfield.bfld[idx(i,k)] *= BP.Bfact;
      //debug
      /*if( Bfield.bfld[idx(i,k)] !=0.0){
        
        double tempR = BP.rmin + i*BP.delr;
        double temptheta = BP.tetmin + k*BP.dtet;
        msg<<"NonZeroP :"<<tempR*cos(temptheta / 180.0 * pi)<<" "<<tempR*sin(temptheta / 180.0 * pi)<<endl;
        }*/
      //end debug
    }
  }
  fclose(f);

  msg<<"Field Map read successfully!"<<endl;
}


void Cyclotron::getDimensions(double &zBegin, double &zEnd) const
{

}

