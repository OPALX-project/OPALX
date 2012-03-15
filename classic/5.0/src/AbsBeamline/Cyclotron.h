#ifndef CLASSIC_Cyclotron_HH
#define CLASSIC_Cyclotron_HH

// ------------------------------------------------------------------------
// $RCSfile: Cyclotron.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1.2.1 $
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
// $Date: 2004/11/12 18:57:53 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Component.h"
#include "BeamlineGeometry/PlanarArcGeometry.h"
#include "Fields/BMultipoleField.h"

struct BfieldData
{
  std::string filename;
  // known from file: field and three theta derivatives
  double* bfld;   //Bz
  double* dbt;    //dBz/dtheta
  double* dbtt;   //d2Bz/dtheta2
  double* dbttt;  //d3Bz/dtheta3

  // to be calculated in getdiffs: all other derivatives:
  double* dbr;    // dBz/dr
  double* dbrr;   // ...
  double* dbrrr;

  double* dbrt;
  double* dbrrt;
  double* dbrtt;

  // used to get (Br,Btheta,Bz) at any off-plane point
  double* f2;  // for Bz
  double* f3;  // for Br
  double* g3;  // for Btheta

  // Grid-Size
  //need to be read from inputfile.
  int nrad, ntet;

  // one more grid line is stored in azimuthal direction:
  int ntetS;

  // total grid points number.
  int ntot;

  // Mean and Maximas
  double bacc, dbtmx, dbttmx, dbtttmx; 

  
};

struct BPositions
{
  // this 4 parameters are need to be read from field file.
  double  rmin, delr;
  double  tetmin, dtet;
  
  // Radii and step width of initial Grid
  double* rarr;
 
  //  int     ThetaPeriodicity; // Periodicity of Magnetic field
  double  Bfact;	    // MULTIPLICATION FACTOR FOR MAGNETIC FIELD
};


// Class Cyclotron
// ------------------------------------------------------------------------
/// Interface for a Cyclotron.
//  This class defines the abstract interface for a Cyclotron.

class Cyclotron: public Component {

public:

  /// Constructor with given name.
  explicit Cyclotron(const string &name);

  Cyclotron();
  Cyclotron(const Cyclotron &);
  virtual ~Cyclotron();

  /// Apply visitor to Cyclotron.
  virtual void accept(BeamlineVisitor &) const;

  /// Get number of slices.
  //  Slices and stepsize used to determine integration step.
  virtual double getSlices() const = 0;

  /// Get stepsize.
  //  Slices and stepsize used to determine integration step.
  virtual double getStepsize() const = 0;
 
  void setFieldMapFN(string fmapfn);
  virtual string getFieldMapFN();

  void setType(string t);
  virtual string getType();

  void setCyclHarm(double h );
  virtual double getCyclHarm();

  void setRfFrequ(double f );
  virtual double getRfFrequ();

  void setSymmetry(double symmetry);
  virtual double getSymmetry();

  void   setRinit(double rinit);
  virtual double getRinit();

  void   setPRinit(double prinit);
  virtual double getPRinit();

  void   setPHIinit(double phiinit);
  virtual double getPHIinit();

  virtual bool apply(const int &i, const double &t, double E[], double B[]);

  virtual bool apply(const int &i, const double &t, Vector_t &E, Vector_t &B);
  
  virtual bool apply(const Vector_t &R, const double &t, Vector_t &E, Vector_t &B);

  virtual void initialise(const PartBunch *bunch, double &startField, double &endField, const double &scaleFactor);

  virtual void initialise(const PartBunch *bunch, const double &scaleFactor);

  virtual void finalise();

  virtual void rescaleFieldMap(const double &scaleFactor);

  virtual bool bends() const;

  virtual double getRmax();

  virtual double getRmin();

 private:

  string fmapfn_m; /* stores the filename of the fieldmap */
  double rffrequ_m; 
  double symmetry_m;

  double rinit_m;
  double prinit_m;
  double phiinit_m;

  string type_m; /* what type of field we use */
  double harm_m;

  
  // Not implemented.
  void operator=(const Cyclotron &);

  // object of Matrics including magnetic field map and its derivates   
  BfieldData Bfield;

  // object of parameters about the map grid
  BPositions BP;

  void   getdiffs();

  double gutdf5d(double *f, double dx, const int kor, const int krl, const int lpr);

  void   initR(double rmin, double dr, int nrad);

  void   getFieldFromFile(const double &scaleFactor);
  void   getFieldFromFile_Carbon(const double &scaleFactor);

  inline int idx(int irad, int ktet) {return (ktet+Bfield.ntetS*irad);}

};

#endif // CLASSIC_Cyclotron_HH
