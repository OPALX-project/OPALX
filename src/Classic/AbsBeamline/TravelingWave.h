#ifndef CLASSIC_TravelingWave_HH
#define CLASSIC_TravelingWave_HH

// ------------------------------------------------------------------------
// $RCSfile: TravelingWave.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: TravelingWave
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


// Class TravelingWave
// ------------------------------------------------------------------------
/// Interface for RF cavity.
//  Class TravelingWave defines the abstract interface for RF cavities.


class TravelingWave: public Component {

public:

    enum CavityType { SW, TW };
    /// Constructor with given name.
    explicit TravelingWave(const std::string &name);

    TravelingWave();
    TravelingWave(const TravelingWave &);
    virtual ~TravelingWave();

    /// Apply visitor to TravelingWave.
    virtual void accept(BeamlineVisitor &) const;

    /// Get RF amplitude.
    virtual double getAmplitude() const = 0;

    /// Get RF frequencey.
    virtual double getFrequency() const = 0;

    /// Get RF phase.
    virtual double getPhase() const = 0;

    /// Set the name of the field map
    void setFieldMapFN(std::string fmapfn);

    /*   void setExitFieldMapFN(std::string fn); */

    /*   void setExitFieldMapFN(std::string fn); */

    std::string getFieldMapFN() const;

    void setFast(bool fast);

    bool getFast() const;

    void setAutophaseVeto(bool veto = true);

    bool getAutophaseVeto() const;

    void setAmplitudem(double vPeak);
    void setAmplitudeError(double voltError);
    double getAmplitudem() const;
    double getAmplitudeError() const;

    void setFrequencym(double freq);
    double getFrequencym() const;

    void setPhasem(double phase);
    void setPhaseError(double phaseError);
    double getPhasem() const;
    double getPhaseError() const;

    void updatePhasem(double phase);

    void setNumCells(int NumCells);

    double getAutoPhaseEstimate(const double & E0, const double & t0, const double & q, const double & m);

    std::pair<double, double> trackOnAxisParticle(const double & p0,
                                                  const double & t0,
                                                  const double & dt,
                                                  const double & q,
                                                  const double & mass);

    virtual void addKR(int i, double t, Vector_t &K);

    virtual void addKT(int i, double t, Vector_t &K);

    virtual bool apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B);

    virtual bool apply(const Vector_t &R, const Vector_t &P, const double &t, Vector_t &E, Vector_t &B);

    virtual bool applyToReferenceParticle(const Vector_t &R, const Vector_t &P, const double &t, Vector_t &E, Vector_t &B);

    virtual void initialise(PartBunch *bunch, double &startField, double &endField);

    virtual void finalise();

    virtual bool bends() const;

    virtual void goOnline(const double &kineticEnergy);

    virtual void goOffline();

    virtual ElementBase::ElementType getType() const;

    virtual void getDimensions(double &zBegin, double &zEnd) const;

    virtual bool isInside(const Vector_t &r) const;

    virtual void setDesignEnergy(double ekin);
    virtual double getDesignEnergy() const;
private:
    std::string CoreFilename_m;             /**< The name of the inputfile*/
    /*   std::string EntryFilename_m; */
    /*   std::string ExitFilename_m; */

    Fieldmap *CoreFieldmap_m;
    /*   Fieldmap *EntryFringeField_m; */
    /*   Fieldmap *ExitFringeField_m; */

    double scale_m;              /**< scale multiplier*/
    double scaleCore_m;
    double scaleError_m;
    double scaleCoreError_m;

    double phase_m;              /**< phase shift of time varying field(degrees)*/
    double phaseCore1_m;
    double phaseCore2_m;
    double phaseExit_m;
    double phaseError_m;

    double frequency_m;          /**< Read in frequency of time varying field(MHz)*/

    double length_m;
    double startCoreField_m;         /**< starting point of field(m)*/
    double startExitField_m;
    double mappedStartExitField_m;

    double PeriodLength_m;
    int NumCells_m;
    double CellLength_m;
    double Mode_m;

    bool fast_m;

    bool autophaseVeto_m;
    double designEnergy_m;

    inline double getdE(const int & i,
                        const int & I,
                        const std::vector<double> & t,
                        const double & phi,
                        const std::vector<std::pair<double, double> > & F) const;

    inline double getdT(const int & i,
                        const int & I,
                        const std::vector<double> & E,
                        const std::vector<std::pair<double, double> > & F,
                        const double mass) const;

    inline double getdA(const int & i,
                        const int & I,
                        const std::vector<double> & t,
                        const double & phi,
                        const std::vector<std::pair<double, double> > & F) const;

    inline double getdB(const int & i,
                        const int & I,
                        const std::vector<double> & t,
                        const double & phi,
                        const std::vector<std::pair<double, double> > & F) const;
    // Not implemented.
    void operator=(const TravelingWave &);
};

inline
void TravelingWave::updatePhasem(double phase) {
    setPhasem(phase);
}

double TravelingWave::getdE(const int & i,
                            const int & I,
                            const std::vector<double> & t,
                            const double & phi,
                            const std::vector<std::pair<double, double> > & F) const {
    return (F[I].first - F[I-1].first) / (frequency_m * frequency_m * (t[i] - t[i-1]) * (t[i] - t[i-1])) *
        (frequency_m * (t[i] - t[i-1]) * (F[I].second * sin(frequency_m * t[i] + phi) - F[I-1].second * sin(frequency_m * t[i-1] + phi)) +
         (F[I].second - F[I-1].second) * (cos(frequency_m * t[i] + phi) - cos(frequency_m * t[i-1] + phi)));
}

double TravelingWave::getdT(const int & i,
                            const int & I,
                            const std::vector<double> & E,
                            const std::vector<std::pair<double, double> > & F,
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
    return (F[I].first - F[I-1].first) *
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

double TravelingWave::getdA(const int & i,
                            const int & I,
                            const std::vector<double> & t,
                            const double & phi,
                            const std::vector<std::pair<double, double> > & F) const {
    double dt = t[i] - t[i-1];
    return (F[I].first - F[I-1].first) / (frequency_m * frequency_m * dt * dt) *
        (frequency_m * dt * (F[I].second * cos(frequency_m * t[i] + phi) - F[I-1].second * cos(frequency_m * t[i-1] + phi)) -
         (F[I].second - F[I-1].second) * (sin(frequency_m * t[i] + phi) - sin(frequency_m * t[i-1] + phi)));
}

double TravelingWave::getdB(const int & i,
                            const int & I,
                            const std::vector<double> & t,
                            const double & phi,
                            const std::vector<std::pair<double, double> > & F) const {
    double dt = t[i] - t[i-1];
    return (F[I].first - F[I-1].first) / (frequency_m * frequency_m * dt * dt) *
        (frequency_m * dt * (F[I].second * sin(frequency_m * t[i] + phi) - F[I-1].second * sin(frequency_m * t[i-1] + phi)) +
         (F[I].second - F[I-1].second) * (cos(frequency_m * t[i] + phi) - cos(frequency_m * t[i-1] + phi)));
}


inline
void TravelingWave::setDesignEnergy(double ekin)
{
    designEnergy_m = ekin;
}

inline
double TravelingWave::getDesignEnergy() const
{
    return designEnergy_m;
}

inline
void TravelingWave::setFieldMapFN(std::string fn) {
    CoreFilename_m = fn;
}

inline
std::string TravelingWave::getFieldMapFN() const {
    return CoreFilename_m;
}

inline
void TravelingWave::setAmplitudem(double vPeak) {
    scale_m = vPeak;
}

inline
void TravelingWave::setAmplitudeError(double voltError) {
    scaleError_m = voltError;
}

inline
double TravelingWave::getAmplitudem() const {
    return scale_m;
}

inline
double TravelingWave::getAmplitudeError() const {
    return scaleError_m;
}

inline
void TravelingWave::setFrequencym(double freq) {
    frequency_m = freq;
}

inline
double TravelingWave::getFrequencym() const {
    return frequency_m;
}

inline
void TravelingWave::setPhasem(double phase) {
    using Physics::pi;

    phase_m = phase;
    phaseCore1_m = phase_m + pi * Mode_m / 2.0;
    phaseCore2_m = phase_m + pi * Mode_m * 1.5;
    phaseExit_m = phase_m - 2.0 * pi * ((NumCells_m - 1) * Mode_m - std::floor((NumCells_m - 1) * Mode_m));
}

inline
void TravelingWave::setPhaseError(double phaseError) {
    phaseError_m = phaseError;
}

inline
double TravelingWave::getPhasem() const {
    return phase_m;
}

inline
double TravelingWave::getPhaseError() const {
    return phaseError_m;
}

inline
void TravelingWave::setNumCells(int NumCells) {
    NumCells_m = NumCells;
}

inline
void TravelingWave::setFast(bool fast) {
    fast_m = fast;
}

inline
bool TravelingWave::getFast() const {
    return fast_m;
}

inline
void TravelingWave::setAutophaseVeto(bool veto) {
    autophaseVeto_m = veto;
}

inline
bool TravelingWave::getAutophaseVeto() const {
    return autophaseVeto_m;
}

#endif // CLASSIC_TravelingWave_HH
