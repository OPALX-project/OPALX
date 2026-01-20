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
// $Date: 20.1.2026 $
// $Author: PSI $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Component.h"
#include "Algorithms/CoordinateSystemTrafo.h"

class Fieldmap;

/**
 * @class Solenoid
 * @brief Abstract class for a solenoid magnet.
 */
class Solenoid : public Component {
public:
/* ============================== Constructors ============================== */
    explicit Solenoid(const std::string& name);
    Solenoid();
    Solenoid(const Solenoid&);
    virtual ~Solenoid();
/* ========================================================================== */
/* ============================== Apply Functions =========================== */
    /**
     * @brief apply the solenoid field to all particles in the bunch
     * 
     * @returns true if at least one particle is lost, false otherwise
     */
    virtual bool apply() override;

    /**
     * @brief apply the solenoid field to particle i
     * 
     * @param i Particle index
     * @param t Time
     * @param E Electric Field
     * @param B Magnetic Field
     * 
     * @returns true if particle is lost, false otherwise
     */
    virtual bool apply(
        const size_t& i, 
        const double& t, 
        Vector_t<double, 3>& E, 
        Vector_t<double, 3>& B) override;

    /**
     * @brief Apply to particle with position R and momentum P
     * 
     * @param R Position
     * @param P Momentum
     * @param t Time
     * @param E Electric Field
     * @param B Magnetic Field
     * 
     * @returns true if particle is lost, false otherwise
     */
    virtual bool apply(
        const Vector_t<double, 3>& R, 
        const Vector_t<double, 3>& P, 
        const double& t,
        Vector_t<double, 3>& E, 
        Vector_t<double, 3>& B) override;

    /**
     * @brief Apply to reference particle with position R and momemtum P
     * 
     * @param R Position
     * @param P Momentum
     * @param t Time
     * @param E Electric Field
     * @param B Magnetic Field
     * 
     * @returns true if particle is lost, false otherwise
     */
    virtual bool applyToReferenceParticle(
        const Vector_t<double, 3>& R, 
        const Vector_t<double, 3>& P, 
        const double& t,
        Vector_t<double, 3>& E, 
        Vector_t<double, 3>& B) override;
/* ========================================================================== */
/* ============================== Functions ================================= */
    /// @brief Apply visitor to Solenoid
    virtual void accept(BeamlineVisitor&) const override;

    /// @brief Get solenoid field Bz
    virtual double getBz() const = 0;

    /// @brief Set the strength scaling factor ks
    void setKS(double ks);

    /// @brief Set the strength scaling error dks
    void setDKS(double ks);

    /**
     * @brief initialise the solenoid element
     * 
     * @param bunch Particle bunch
     * @param startField Starting position of the field 
     * @param endField Ending position of the field
     */
    virtual void initialise(
        PartBunch_t* bunch, 
        double& startField, 
        double& endField) override;

    /// @note not implemented
    virtual void finalise() override;

    virtual bool bends() const override;

    /// @brief Load field map and go online
    virtual void goOnline(const double& kineticEnergy) override;

    /// @brief Free field map and go offline
    virtual void goOffline() override;

    /// @brief Assign the field filename
    void setFieldMapFN(std::string fn);

    void setFast(bool fast);

    bool getFast() const;

    virtual ElementType getType() const override;

    /// @brief Get the begin and end positions of the element 
    virtual void getDimensions(double& zBegin, double& zEnd) const override;
    
    /// @brief Get the begin and end positions of the element
    virtual void getElementDimensions(double& zBegin, double& zEnd) const override;

    /// @brief Check if position r is inside the field map
    virtual bool isInside(const Vector_t<double, 3>& r) const override;

    /// @brief Get the coordinate transformation to the begin of the element
    virtual CoordinateSystemTrafo getEdgeToBegin() const override;

    /// @brief Get the coordinate transformation to the end of the element
    virtual CoordinateSystemTrafo getEdgeToEnd() const override;

private:
/* ========================================================================== */
/* ============================== Variables ================================= */
    
    /// Name of the field map file
    std::string filename_m; 

    /// Fieldmap pointer
    Fieldmap* fieldmap_m;

    /// Scale multiplier
    double scale_m; 
    
    /// Scale error multiplier
    double scaleError_m;

    /// Starting point of the field
    double startField_m; 

    /// Fast tracking flag @note currently not implemented
    bool fast_m;

    /// @note not implemente 
    void operator=(const Solenoid&);
};

/**
 * @brief Get the coordinate transformation to the start of the element
 * 
 * @returns CoordinateSystemTrafo to the begin of the element
 */
inline CoordinateSystemTrafo Solenoid::getEdgeToBegin() const {
    CoordinateSystemTrafo ret(
        Vector_t<double, 3>(0, 0, startField_m), 
        Quaternion(1, 0, 0, 0));
    return ret;
}

/**
 * @brief Get the coordinate transformation to the end of the element
 * 
 * @returns CoordinateSystemTrafo to the end of the element
 */
inline CoordinateSystemTrafo Solenoid::getEdgeToEnd() const {
    CoordinateSystemTrafo ret(
        Vector_t<double, 3>(0, 0, startField_m + getElementLength()), 
        Quaternion(1, 0, 0, 0));
    return ret;
}
#endif  // CLASSIC_Solenoid_HH
