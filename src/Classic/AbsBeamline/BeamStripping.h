#ifndef CLASSIC_BeamStripping_HH
#define CLASSIC_BeamStripping_HH

// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: BeamStripping
//   Defines the abstract interface for a beam BeamStripping.
//   *** MISSING *** BeamStripping interface is still incomplete.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Component.h"

#include "gsl/gsl_spline.h"
#include "gsl/gsl_interp.h"

#include <string>
#include <vector>

using namespace std;

class BeamlineVisitor;
class LossDataSink;

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

    virtual bool applyToReferenceParticle(const Vector_t &R, const Vector_t &P, const double &t, Vector_t &E, Vector_t &B);

    virtual bool checkBeamStripping(PartBunchBase<double, 3> *bunch, const int turnnumber, const double t, const double tstep);

    virtual bool checkBeamStripping(Vector_t r, Vector_t rmin, Vector_t rmax);

    virtual void initialise(PartBunchBase<double, 3> *bunch, double &startField, double &endField);

    virtual void initialise(PartBunchBase<double, 3> *bunch);

    virtual void finalise();

    virtual bool bends() const;

    virtual void goOnline(const double &kineticEnergy);

    virtual void goOffline();

    virtual ElementBase::ElementType getType() const;

    virtual void getDimensions(double &zBegin, double &zEnd) const;

    void print();

    string  getBeamStrippingShape();
    void setOutputFN(string fn);
    string getOutputFN();

    unsigned int getLosses() const;

    int checkPoint(const double &x, const double &y, const double &z);

    // --------Cyclotron beam stripping

    void setPressure(double pressure) ;
    double getPressure() const;

    void setTemperature(double temperature) ;
    double getTemperature() const;

    void setCrossSection(vector<double> sigma);
    vector<double> getCrossSection() const;

    void setEnergyCS(vector<double> energycs);
    vector<double> getEnergyCS() const;

    void setMinR(double r);
    double getMinR() const;
    void setMaxR(double r);
    double getMaxR() const;

    void setMinZ(double z);
    double getMinZ() const;
    void setMaxZ(double z);
    double getMaxZ() const;

//    double CrossSection();


private:

    // Not implemented.
    void operator=(const BeamStripping &);

    string filename_m;               /**< The name of the outputfile*/
    bool informed_m;

    //parameters for BeamStripping
    double pressure_m;
    double temperature_m;
    vector<double> sigma_m;
	vector<double> energycs_m;

    double minr_m;
    double maxr_m;
    double minz_m;
    double maxz_m;
    double rpos;
    double zpos;

//	double CS;
//	double Eng;

    unsigned int losses_m;
    unique_ptr<LossDataSink> lossDs_m;

    ParticleMatterInteractionHandler *parmatintbst_m;
};

inline
unsigned int BeamStripping::getLosses() const {
    return losses_m;
}

#endif // CLASSIC_BeamStripping_HH
