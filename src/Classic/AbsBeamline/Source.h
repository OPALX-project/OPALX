#ifndef CLASSIC_SOURCE_HH
#define CLASSIC_SOURCE_HH

#include "AbsBeamline/Component.h"
#include "Distribution/Distribution.h"
class OpalBeamline;

template <class T, unsigned Dim>
class PartBunchBase;

class Source: public Component {

public:

    /// Constructor with given name.
    explicit Source(const std::string &name);

    Source();
    Source(const Source &);
    virtual ~Source();

    /// Apply visitor to Source.
    virtual void accept(BeamlineVisitor &) const;

    virtual void addKR(int i, double t, Vector_t &K);

    virtual void addKT(int i, double t, Vector_t &K);

    virtual bool apply(const double &t);

    virtual void initialise(PartBunchBase<double, 3> *bunch, double &startField, double &endField);

    virtual void finalise();

    virtual bool bends() const;

    virtual void goOnline(const double &kineticEnergy);

    virtual void goOffline();

    virtual ElementBase::ElementType getType() const;

    virtual void getDimensions(double &zBegin, double &zEnd) const;

    void setDistribution(std::vector<std::string> distNames);

    void setBeamline(OpalBeamline *beamline);
private:

    double ElementEdge_m;
    double startField_m;           /**< startingpoint of field, m*/
    double endField_m;

    std::vector<Distribution*> distrs_m;

    OpalBeamline *beamline_m;
    // Not implemented.
    void operator=(const Source &);

    static const std::string defaultDistribution;
};

inline
void Source::setBeamline(OpalBeamline *beamline) {
    beamline_m = beamline;
}
#endif // CLASSIC_SOURCE_HH