#ifndef OPAL_EMISSIONSOURCE_HH
#define OPAL_EMISSIONSOURCE_HH

#include "AbstractObjects/Definition.h"
#include "Ippl.h"

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
        EMISSIONMODEL,
        SIZE
    };

    /// Exemplar constructor.
    EmissionSource();

    virtual ~EmissionSource();

    virtual bool canReplaceBy(Object* object);
    virtual EmissionSource* clone(const std::string& name);
    virtual void execute();

    static EmissionSource* find(const std::string& name);

    std::string getDistributionName() const;
    ippl::Vector<double, 3> getR0() const;
    ippl::Vector<double, 3> getP0() const;
    double getT0() const;
    std::string getEmissionModel() const;

private:
    EmissionSource(const EmissionSource&);
    void operator=(const EmissionSource&);

    EmissionSource(const std::string& name, EmissionSource* parent);
};

#endif  // OPAL_EMISSIONSOURCE_HH
