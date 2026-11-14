#ifndef OPAL_EMISSIONSOURCE_HH
#define OPAL_EMISSIONSOURCE_HH

#include "AbstractObjects/Definition.h"
#include "Ippl.h"
#include "Track/ObjectValues.h"

#include <string>

class EmissionSource : public Definition {
public:
    /// The attributes of class EmissionSource.
    enum {
        DISTRIBUTION,
        R0X,
        R0Y,
        R0Z,
        P0X,
        P0Y,
        P0Z,
        T0,
        ZEROFACE_R0Z,
        SHIFTED_GREENS_FUNCTION,
        ZEROFACEPLANEDUMP,
        ZEROFACE_MAXSTEPS,
        EMISSIONMODEL,
        EKIN,
        SIZE
    };

    /// Exemplar constructor.
    EmissionSource();

    virtual ~EmissionSource();

    virtual bool canReplaceBy(Object* object);
    virtual void execute();

    static EmissionSource* find(const std::string& name);

    static std::string getDistributionName(const ObjectValues& values);
    static ippl::Vector<double, 3> getR0(const ObjectValues& values);
    static ippl::Vector<double, 3> getP0(const ObjectValues& values);
    static double getT0(const ObjectValues& values);
    static bool getZeroFaceR0Z(const ObjectValues& values);
    static bool getShiftedGreensFunction(const ObjectValues& values);
    static int getZeroFacePlaneDumpFrequency(const ObjectValues& values);
    static int getZerofaceMaxSteps(const ObjectValues& values);
    static std::string getEmissionModel(const ObjectValues& values);
    static double getKineticEnergy(const ObjectValues& values);

private:
    EmissionSource(const EmissionSource&);
    void operator=(const EmissionSource&);
};

#endif  // OPAL_EMISSIONSOURCE_HH
