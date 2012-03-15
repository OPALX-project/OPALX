#ifndef CLASSIC_RFCavity_HH
#define CLASSIC_RFCavity_HH

// ------------------------------------------------------------------------
// $RCSfile: RFCavity.h,v $
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


#include "AbsBeamline/Component.h"

class Fieldmap;

// Class RFCavity
// ------------------------------------------------------------------------
/// Interface for RF cavity.
//  Class RFCavity defines the abstract interface for RF cavities.


class RFCavity: public Component {

public:

    enum CavityType { SW, SGSW };
    /// Constructor with given name.
    explicit RFCavity(const string &name);

    RFCavity();
    RFCavity(const RFCavity &);
    virtual ~RFCavity();

    /// Apply visitor to RFCavity.
    virtual void accept(BeamlineVisitor &) const;

    /// Get RF amplitude.
    virtual double getAmplitude() const = 0;

    /// Get RF frequencey.
    virtual double getFrequency() const = 0;

    /// Get RF phase.
    virtual double getPhase() const = 0;

    /// Set the name of the field map
    void setFieldMapFN(string fmapfn);

    string getFieldMapFN() const;

    void setAmplitudem(double vPeak);

    void setFrequencym(double freq);

    void setPhasem(double phase);

    double getPhasem() const;

    void setCavityType(string type);

    string getCavityType() const;

    void setFast(bool fast);

    bool getFast() const;
  
    void TestField();
  
    virtual void addKR(int i, double t, Vector_t &K);

    virtual void addKT(int i, double t, Vector_t &K);

    virtual bool apply(const int &i, const double &t, double E[], double B[]);

    virtual bool apply(const int &i, const double &t, Vector_t &E, Vector_t &B);

    virtual bool apply(const Vector_t &R, const Vector_t &centroid, const double &t, Vector_t &E, Vector_t &B);

    virtual void initialise(PartBunch *bunch, double &startField, double &endField, const double &scaleFactor);

    virtual void initialise(PartBunch *bunch, const double &scaleFactor);

    virtual void finalise();

    virtual bool bends() const;

    virtual void goOnline();

    virtual void goOffline();

    void setRmin(double rmin);

    void setRmax(double rmax);

    void setAzimuth(double angle);

    void setPerpenDistance(double pdis);

    void setGapWidth(double gapwidth);

    void setPhi0(double phi0);

    virtual double getRmin() const;

    virtual double getRmax() const;

    virtual double getAzimuth() const;

    virtual double getCosAzimuth() const;

    virtual double getSinAzimuth() const;

    virtual double getPerpenDistance() const;

    virtual double getGapWidth() const;

    virtual double getPhi0() const;

    virtual void setComponentType(string name);

    virtual string getComponentType()const;

    virtual double getCycFrequency()const;
  
    double calcDeriv(double z, double t);

    double calcKR(double z, double t, double beta);

    void getMomentaKick(const double normalRadius,double momentum[], const double t, const double dtCorrt, const int PID, const double restMass );

    double spline(double z, double *za);

    virtual const string& getType() const;

    virtual void getDimensions(double &zBegin, double &zEnd) const;

private:
    string filename_m;             /**< The name of the inputfile*/
    Fieldmap *fieldmap_m;
    double scale_m;              /**< scale multiplier*/
    double phase_m;              /**< phase shift of time varying field(degrees)*/
    double frequency_m;          /**< Read in frequency of time varying field(MHz)*/
    double ElementEdge_m;
    double startField_m;         /**< starting point of field(m)*/
    double endField_m;

    CavityType type_m;

    bool fast_m;

    double rmin_m;
    double rmax_m;
    double angle_m;
    double sinAngle_m;
    double cosAngle_m;
    double pdis_m;
    double gapwidth_m;
    double phi0_m;

    double* RNormal_m;
    double* VrNormal_m;
    double* DvDr_m;
    int num_points_m;


    // Not implemented.
    void operator=(const RFCavity &);
};

#endif // CLASSIC_RFCavity_HH
