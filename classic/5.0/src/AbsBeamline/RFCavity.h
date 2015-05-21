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
#include "Physics/Physics.h"

class Fieldmap;

// Class RFCavity
// ------------------------------------------------------------------------
/// Interface for RF cavity.
//  Class RFCavity defines the abstract interface for RF cavities.


class RFCavity: public Component {

public:

    enum CavityType { SW, SGSW };
    /// Constructor with given name.
    explicit RFCavity(const std::string &name);

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

    void dropFieldmaps();

    /// Set the name of the field map
    void setFieldMapFN(std::string fmapfn);

    std::string getFieldMapFN() const;

    void setAmplitudem(double vPeak);

    void setFrequencym(double freq);

    void setFrequency(double freq);

    double getFrequencym() const ;

    void setPhasem(double phase);

    void updatePhasem(double phase);

    double getPhasem() const;

    void setCavityType(std::string type);

    std::string getCavityType() const;

    void setFast(bool fast);

    bool getFast() const;

    void setAutophaseVeto(bool veto);

    bool getAutophaseVeto() const;

    double getAutoPhaseEstimate(const double & E0, const double & t0, const double & q, const double & m);

    std::pair<double, double> trackOnAxisParticle(const double & p0, 
                                                  const double & t0, 
                                                  const double & dt, 
                                                  const double & q, 
                                                  const double & mass);

    virtual void addKR(int i, double t, Vector_t &K);

    virtual void addKT(int i, double t, Vector_t &K);

    virtual bool apply(const size_t &i, const double &t, double E[], double B[]);

    virtual bool apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B);

    virtual bool apply(const Vector_t &R, const Vector_t &centroid, const double &t, Vector_t &E, Vector_t &B);

    virtual void initialise(PartBunch *bunch, double &startField, double &endField, const double &scaleFactor);

    virtual void initialise(PartBunch *bunch, const double &scaleFactor);

    virtual void finalise();

    virtual bool bends() const;

    virtual void goOnline(const double &kineticEnergy);

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

    virtual void setComponentType(std::string name);

    virtual std::string getComponentType()const;

    virtual double getCycFrequency()const;

    void getMomentaKick(const double normalRadius, double momentum[], const double t, const double dtCorrt, const int PID, const double restMass,const int chargenumber);

    double spline(double z, double *za);

    virtual ElementBase::ElementType getType() const;

    virtual void getDimensions(double &zBegin, double &zEnd) const;

private:
    std::string filename_m;             /**< The name of the inputfile*/
    std::vector<std::string> multiFilenames_m;
    std::vector<Fieldmap*> multiFieldmaps_m;
    double scale_m;              /**< scale multiplier*/
    std::vector<double> multiScales_m;
    double phase_m;              /**< phase shift of time varying field(degrees)*/
    std::vector<double> multiPhases_m;
    double frequency_m;          /**< Read in frequency of time varying field(MHz)*/
    std::vector<double> multiFrequencies_m;
    double ElementEdge_m;
    double startField_m;         /**< starting point of field(m)*/
    double endField_m;
    std::vector<std::pair<double, double> > multi_start_end_field_m;

    CavityType type_m;

    bool fast_m;
    bool autophaseVeto_m;

    double rmin_m;
    double rmax_m;
    double angle_m;
    double sinAngle_m;
    double cosAngle_m;
    double pdis_m;
    double gapwidth_m;
    double phi0_m;

    std::unique_ptr<double[]> RNormal_m;
    std::unique_ptr<double[]> VrNormal_m;
    std::unique_ptr<double[]> DvDr_m;
    int num_points_m;

    inline size_t numFieldmaps() const {
        return multiFilenames_m.size();
    }

    
    double getdE(const int & i, 
                 const std::vector<double> & t, 
                 const double & dz,
                 const double & phi,
                 const double & frequency,
                 const std::vector<double> & F) const;

    double getdT(const int & i, 
                 const std::vector<double> & E,
                 const double & dz,
                 const double mass) const;
    
    double getdA(const int & i,
                 const std::vector<double> & t, 
                 const double & dz,
                 const double & phi,
                 const double & frequency,
                 const std::vector<double> & F) const;

    double getdB(const int & i,
                 const std::vector<double> & t, 
                 const double & dz,
                 const double & phi,
                 const double & frequency,
                 const std::vector<double> & F) const;

    // Not implemented.
    void operator=(const RFCavity &);
};

inline
double RFCavity::getdE(const int & i,
                       const std::vector<double> & t, 
                       const double & dz,
                       const double & phi,
                       const double & frequency,
                       const std::vector<double> & F) const {
    return dz / (frequency * frequency * (t[i] - t[i-1]) * (t[i] - t[i-1])) * 
        (frequency * (t[i] - t[i-1]) * (F[i] * sin(frequency * t[i] + phi) - F[i-1] * sin(frequency * t[i-1] + phi)) +
         (F[i] - F[i-1]) * (cos(frequency * t[i] + phi) - cos(frequency * t[i-1] + phi)));
}

inline
double RFCavity::getdT(const int & i, 
                       const std::vector<double> & E,
                       const double & dz,
                       const double mass) const {
    double gamma1  = 1. + (19. * E[i-1] + 1. * E[i]) / (20. * mass);
    double gamma2  = 1. + (17. * E[i-1] + 3. * E[i]) / (20. * mass);
    double gamma3  = 1. + (15. * E[i-1] + 5. * E[i]) / (20. * mass);
    double gamma4  = 1. + (13. * E[i-1] + 7. * E[i]) / (20. * mass);
    double gamma5  = 1. + (11. * E[i-1] + 9. * E[i]) / (20. * mass);
    double gamma6  = 1. + (9. * E[i-1] + 11. * E[i]) / (20. * mass);
    double gamma7  = 1. + (7. * E[i-1] + 13. * E[i]) / (20. * mass);
    double gamma8  = 1. + (5. * E[i-1] + 15. * E[i]) / (20. * mass);
    double gamma9  = 1. + (3. * E[i-1] + 17. * E[i]) / (20. * mass);
    double gamma10 = 1. + (1. * E[i-1] + 19. * E[i]) / (20. * mass);
    return dz *
        (1. / sqrt(1. - 1. / (gamma1 * gamma1)) +
         1. / sqrt(1. - 1. / (gamma2 * gamma2)) +
         1. / sqrt(1. - 1. / (gamma3 * gamma3)) +
         1. / sqrt(1. - 1. / (gamma4 * gamma4)) +
         1. / sqrt(1. - 1. / (gamma5 * gamma5)) +
         1. / sqrt(1. - 1. / (gamma6 * gamma6)) +
         1. / sqrt(1. - 1. / (gamma7 * gamma7)) +
         1. / sqrt(1. - 1. / (gamma8 * gamma8)) +
         1. / sqrt(1. - 1. / (gamma9 * gamma9)) +
         1. / sqrt(1. - 1. / (gamma10 * gamma10))) / (10. * Physics::c);
}

inline
double RFCavity::getdA(const int & i,
                       const std::vector<double> & t, 
                       const double & dz,
                       const double & phi,
                       const double & frequency,
                       const std::vector<double> & F) const {
    double dt = t[i] - t[i-1];
    return dz / (frequency * frequency * dt * dt) *
        (frequency * dt * (F[i] * cos(frequency * t[i] + phi) - F[i-1] * cos(frequency * t[i-1] + phi)) -
         (F[i] - F[i-1]) * (sin(frequency * t[i] + phi) - sin(frequency * t[i-1] + phi)));
}

inline
double RFCavity::getdB(const int & i,
                       const std::vector<double> & t, 
                       const double & dz,
                       const double & phi,
                       const double & frequency,
                       const std::vector<double> & F) const {
    double dt = t[i] - t[i-1];
    return dz / (frequency * frequency * dt * dt) * 
        (frequency * dt * (F[i] * sin(frequency * t[i] + phi) - F[i-1] * sin(frequency * t[i-1] + phi)) +
         (F[i] - F[i-1]) * (cos(frequency * t[i] + phi) - cos(frequency * t[i-1] + phi)));
}

#endif // CLASSIC_RFCavity_HH
