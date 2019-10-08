#ifndef CLASSIC_BeamStripping_HH
#define CLASSIC_BeamStripping_HH

// Class category: AbsBeamline
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: BeamStripping
//   Defines the abstract interface for a beam BeamStripping
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
// $Date: 2018/11 $
// $Author: PedroCalvo$
// ------------------------------------------------------------------------

#include "AbsBeamline/Component.h"
#include "AbsBeamline/Cyclotron.h"

#include "gsl/gsl_spline.h"
#include "gsl/gsl_interp.h"

#include <string>
#include <vector>

class BeamlineVisitor;
class LossDataSink;

struct PFieldData {
    std::string filename;
    // known from file: field and three theta derivatives
    std::vector<double> pfld;   //Pz

    // Grid-Size
    //need to be read from inputfile.
    int nrad, ntet;

    // one more grid line is stored in azimuthal direction:
    int ntetS;

    // total grid points number.
    int ntot;
};

struct PPositions {
    // these 4 parameters are need to be read from field file.
    double  rmin, delr;
    double  tetmin, dtet;

    // Radii and step width of initial Grid
    std::vector<double> rarr;

    double  Pfact; // MULTIPLICATION FACTOR FOR PRESSURE MAP
};



// Class BeamStripping
// ------------------------------------------------------------------------

class BeamStripping: public Component {

public:

    /// Constructor with given name.
    explicit BeamStripping(const std::string &name);

    BeamStripping();
    BeamStripping(const BeamStripping &rhs);
    virtual ~BeamStripping();

    /// Apply visitor to BeamStripping.
    virtual void accept(BeamlineVisitor &) const;

    virtual bool apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B);

    virtual bool apply(const Vector_t &R, const Vector_t &P, const double &t, Vector_t &E, Vector_t &B);

    virtual bool checkBeamStripping(PartBunchBase<double, 3> *bunch, Cyclotron* cycl, const int turnnumber, const double t, const double tstep);

    virtual bool checkBeamStripping(Vector_t r, Vector_t rmin, Vector_t rmax);

    virtual void initialise(PartBunchBase<double, 3> *bunch, double &startField, double &endField);

    virtual void initialise(PartBunchBase<double, 3> *bunch, const double &scaleFactor);

    virtual void finalise();

    virtual bool bends() const;

    virtual void goOnline(const double &kineticEnergy);

    virtual void goOffline();

    virtual ElementBase::ElementType getType() const;

    virtual void getDimensions(double &zBegin, double &zEnd) const;

    void print();

    std::string  getBeamStrippingShape();
    void setOutputFN(std::string fn);
    std::string getOutputFN();

    unsigned int getLosses() const;

    int checkPoint(const double &x, const double &y, const double &z);
    
    double checkPressure(const double &x, const double &y);

    void setPressure(double pressure) ;
    double getPressure() const;

    void setTemperature(double temperature) ;
    double getTemperature() const;

    void setPressureMapFN(std::string pmapfn);
    virtual std::string getPressureMapFN() const;
    
    void setPScale(double ps);
    virtual double getPScale() const;

    void setResidualGas(std::string gas);
    virtual std::string getResidualGas() const;

    void setStop(bool stopflag);
    virtual bool getStop() const;
    
protected:

    void   initR(double rmin, double dr, int nrad);

    void   getPressureFromFile(const double &scaleFactor);

    inline int idx(int irad, int ktet) {return (ktet + PField.ntetS * irad);}

private:

    std::string filename_m;               /**< The name of the outputfile*/
    bool informed_m;

    ///@{ parameters for BeamStripping
    std::string gas_m;
    double pressure_m;
    std::string pmapfn_m; /// stores the filename of the pressure map
    double pscale_m;      /// a scale factor for the P-field
    double temperature_m; ///< K
    double stop_m;
    ///@}

    ///@{ size limits took from cyclotron
    double minr_m; 
    double maxr_m;
    double minz_m;
    double maxz_m;
    ///@}

    unsigned int losses_m;
    std::unique_ptr<LossDataSink> lossDs_m;

    ParticleMatterInteractionHandler *parmatint_m;
    
protected:
    // object of Matrices including pressure field map and its derivates
    PFieldData PField;

    // object of parameters about the map grid
    PPositions PP;
};

inline
unsigned int BeamStripping::getLosses() const {
    return losses_m;
}

#endif // CLASSIC_BeamStripping_HH
