#ifndef CLASSIC_Solenoid_HH
#define CLASSIC_Solenoid_HH

// ------------------------------------------------------------------------
// $RCSfile: Solenoid.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Solenoid
//   Defines the abstract interface for a solenoid magnet.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:32 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Component.h"
#include "Algorithms/CoordinateSystemTrafo.h"

class Fieldmap;

// Class Solenoid
// ------------------------------------------------------------------------
/// Interface for solenoids.
//  Class Solenoid defines the abstract interface for solenoid magnets.

class Solenoid : public Component {
public:
    /// Constructor with given name.
    explicit Solenoid(const std::string& name);

    Solenoid();
    Solenoid(const Solenoid&);
    virtual ~Solenoid();

    /// Apply visitor to Solenoid.
    virtual void accept(BeamlineVisitor&) const override;

    /// Get solenoid field Bz in Teslas.
    virtual double getBz() const = 0;

    void setKS(double ks);
    void setDKS(double ks);

    virtual bool apply(
        const size_t& i, const double& t, Vector_t<double, 3>& E, Vector_t<double, 3>& B) override;

    virtual bool apply(
        const Vector_t<double, 3>& R, const Vector_t<double, 3>& P, const double& t,
        Vector_t<double, 3>& E, Vector_t<double, 3>& B) override;

    virtual bool applyToReferenceParticle(
        const Vector_t<double, 3>& R, const Vector_t<double, 3>& P, const double& t,
        Vector_t<double, 3>& E, Vector_t<double, 3>& B) override;

    virtual void initialise(PartBunch_t* bunch, double& startField, double& endField) override;

    virtual void finalise() override;

    virtual bool bends() const override;

    virtual void goOnline(const double& kineticEnergy) override;

    virtual void goOffline() override;

    //  Assign the field filename.
    void setFieldMapFN(std::string fn);

    void setFast(bool fast);

    bool getFast() const;

    virtual ElementType getType() const override;

    virtual void getDimensions(double& zBegin, double& zEnd) const override;

    virtual bool isInside(const Vector_t<double, 3>& r) const override;

    virtual void getElementDimensions(double& zBegin, double& zEnd) const override;

    virtual CoordinateSystemTrafo getEdgeToBegin() const override;
    virtual CoordinateSystemTrafo getEdgeToEnd() const override;

private:
    //  std::string name;                   /**< The name of the object*/
    std::string filename_m; /**< The name of the inputfile*/

    Fieldmap* fieldmap_m;

    double scale_m;      /**< scale multiplier*/
    double scaleError_m; /**< scale multiplier error*/

    double startField_m; /**< startingpoint of field, m*/

    bool fast_m;
    // Not implemented.
    void operator=(const Solenoid&);
};

inline CoordinateSystemTrafo Solenoid::getEdgeToBegin() const {
    CoordinateSystemTrafo ret(Vector_t<double, 3>(0, 0, startField_m), Quaternion(1, 0, 0, 0));

    return ret;
}

inline CoordinateSystemTrafo Solenoid::getEdgeToEnd() const {
    CoordinateSystemTrafo ret(
        Vector_t<double, 3>(0, 0, startField_m + getElementLength()), Quaternion(1, 0, 0, 0));
    return ret;
}
#endif  // CLASSIC_Solenoid_HH
