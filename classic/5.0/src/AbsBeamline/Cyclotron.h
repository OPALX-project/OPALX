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
// $Date: 2007/08/01 $
// $Author: Yang, Adelmann $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Component.h"
#include "BeamlineGeometry/PlanarArcGeometry.h"
#include "Fields/BMultipoleField.h"

class Fieldmap;
class LossDataSink;

enum BFieldType {PSIBF, CARBONBF,ANSYSBF,AVFEQBF, FFAGBF,BANDRF};

struct BfieldData {
    std::string filename;
    // known from file: field and three theta derivatives
    //~ double *bfld;   //Bz
    //~ double *dbt;    //dBz/dtheta
    //~ double *dbtt;   //d2Bz/dtheta2
    //~ double *dbttt;  //d3Bz/dtheta3
//~ 
    //~ // to be calculated in getdiffs: all other derivatives:
    //~ double *dbr;    // dBz/dr
    //~ double *dbrr;   // ...
    //~ double *dbrrr;
//~ 
    //~ double *dbrt;
    //~ double *dbrrt;
    //~ double *dbrtt;
//~ 
    //~ // used to get (Br,Btheta,Bz) at any off-plane point
    //~ double *f2;  // for Bz
    //~ double *f3;  // for Br
    //~ double *g3;  // for Btheta
//~ 
    std::vector<double> bfld;   //Bz
    std::vector<double> dbt;    //dBz/dtheta
    std::vector<double> dbtt;   //d2Bz/dtheta2
    std::vector<double> dbttt;  //d3Bz/dtheta3

    // to be calculated in getdiffs: all other derivatives:
    std::vector<double> dbr;    // dBz/dr
    std::vector<double> dbrr;   // ...
    std::vector<double> dbrrr;

    std::vector<double> dbrt;
    std::vector<double> dbrrt;
    std::vector<double> dbrtt;

    // used to get (Br,Btheta,Bz) at any off-plane point
    std::vector<double> f2;  // for Bz
    std::vector<double> f3;  // for Br
    std::vector<double> g3;  // for Btheta

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

struct BPositions {
    // this 4 parameters are need to be read from field file.
    double  rmin, delr;
    double  tetmin, dtet;

    // Radii and step width of initial Grid
    std::vector<double> rarr;

    //  int     ThetaPeriodicity; // Periodicity of Magnetic field
    double  Bfact;      // MULTIPLICATION FACTOR FOR MAGNETIC FIELD
};


// Class Cyclotron
// ------------------------------------------------------------------------
/// Interface for a Cyclotron.
//  This class defines the abstract interface for a Cyclotron.

class Cyclotron: public Component {
public:

    /// Constructor with given name.
    explicit Cyclotron(const std::string &name);

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

    void setFieldMapFN(std::string fmapfn);
    virtual std::string getFieldMapFN() const;

    void setRfFieldMapFN(std::vector<std::string> rffmapfn);

    void setType(std::string t);
    virtual const std::string &getType() const;

    virtual void getDimensions(double &zBegin, double &zEnd) const;

    void setCyclHarm(double h);
    virtual double getCyclHarm() const;

    void setRfPhi(std::vector<double> f);

    void setRfFrequ(std::vector<double> f);
    double getRfFrequ() const;

    void setSymmetry(double symmetry);
    virtual double getSymmetry() const;

    void   setRinit(double rinit);
    virtual double getRinit() const;

    void   setPRinit(double prinit);
    virtual double getPRinit() const;

    void   setPHIinit(double phiinit);
    virtual double getPHIinit() const;

    void   setZinit(double zinit);
    virtual double getZinit() const;
    
    void   setPZinit(double zinit);
    virtual double getPZinit() const;

    void   setBScale(double bs);
    virtual double getBScale() const;

    void   setEScale(std::vector<double> bs);

    void   setTCr1(double tcr1);
    virtual  double getTCr1() const;

    void   setTCr2(double tcr2);
    virtual  double getTCr2() const;

    void   setMBtc(double mbtc);
    virtual  double getMBtc() const;

    void   setSLPtc(double slptc);
    virtual  double getSLPtc() const;

    void   setSuperpose(std::vector<bool> flag);
    //    virtual bool getSuperpose() const;

    void setMinR(double r);
    virtual double getMinR() const;
    void setMaxR(double r);
    virtual double getMaxR() const;

    void setMinZ(double z);
    virtual double getMinZ() const;
    void setMaxZ(double z);
    virtual double getMaxZ() const;

    void setFMLowE(double e);
    virtual double getFMLowE() const;
    void setFMHighE(double e);
    virtual double getFMHighE() const;

    virtual bool apply(const size_t &id, const double &t, double E[], double B[]);

    virtual bool apply(const size_t &id, const double &t, Vector_t &E, Vector_t &B);

    virtual bool apply(const Vector_t &R, const Vector_t &centroid, const double &t, Vector_t &E, Vector_t &B);

    virtual void initialise(PartBunch *bunch, double &startField, double &endField, const double &scaleFactor);

    virtual void initialise(PartBunch *bunch, const int &fieldflag, const double &scaleFactor);

    virtual void finalise();

    virtual bool bends() const;

    virtual double getRmax() const;

    virtual double getRmin() const;

private:

    std::string fmapfn_m; /* stores the filename of the fieldmap */
    std::vector<double> rffrequ_m;
    std::vector<double> rfphi_m;
    std::vector<double> escale_m; // a scale factor for the E-field
    std::vector<bool> superpose_m; // electric fields are superposed or not 

    double symmetry_m;

    double rinit_m;
    double prinit_m;
    double phiinit_m;
    double zinit_m;
    double pzinit_m;

    std::string type_m; /* what type of field we use */
    double harm_m;

    double bscale_m; // a scale factor for the B-field

    double tcr1_m;
    double tcr2_m;
    double mbtc_m;
    double slptc_m;

    double minr_m;
    double maxr_m;

    double minz_m;
    double maxz_m;

    double fmLowE_m;
    double fmHighE_m;

    // Not implemented.
    void operator=(const Cyclotron &);

    // object of Matrics including magnetic field map and its derivates
    BfieldData Bfield;

    // object of parameters about the map grid
    BPositions BP;

    BFieldType myBFieldType_m;
    
    // RF field map handler
    //    Fieldmap *RFfield;
    std::vector<Fieldmap *> RFfields_m; 
    std::vector<std::string> RFfilename_m;

    // handling for store the particle out of region 
    std::unique_ptr<LossDataSink> lossDs_m;

    void   getdiffs();

    double gutdf5d(double *f, double dx, const int kor, const int krl, const int lpr);

    void   initR(double rmin, double dr, int nrad);

    void   getFieldFromFile(const double &scaleFactor);
    void   getFieldFromFile_Carbon(const double &scaleFactor);
    void   getFieldFromFile_CYCIAE(const double &scaleFactor);
    void   getFieldFromFile_AVFEQ(const double &scaleFactor);
    void   getFieldFromFile_FFAG(const double &scaleFactor);
    void   getFieldFromFile_BandRF(const double &scaleFactor);

    inline int idx(int irad, int ktet) {return (ktet + Bfield.ntetS * irad);}

};

#endif // CLASSIC_Cyclotron_HH
