/**
 * @file Component.cpp
 * @brief Defines the abstract interface for a single beamline component in the 
 * accelerator model.
 *
 * @details
 * The Component class provides an abstract interface for an arbitrary single 
 * component in a beam line.
 * A component is the basic element in the accelerator model, such as a dipole, 
 * quadrupole, etc.
 */-----------------------------------------------------------------------

#include "AbsBeamline/Component.h"
#include "Utilities/LogicalError.h"

const std::vector<double> Component::defaultAperture_m = 
    std::vector<double>({1e6, 1e6, 1.0});

Component::Component() : Component("") {
}

Component::Component(const Component& right)
    : ElementBase(right), 
    exit_face_slope_m(right.exit_face_slope_m), 
    online_m(right.online_m) {
}

Component::Component(const std::string& name)
    : ElementBase(name), 
    exit_face_slope_m(0.0), 
    RefPartBunch_m(nullptr), online_m(false) {
    setAperture(ApertureType::ELLIPTICAL, defaultAperture_m);
}

Component::~Component() {
}

const ElementBase& Component::getDesign() const {
    return *this;
}

void Component::trackBunch(
    PartBunch_t*, 
    const PartData&, 
    bool, 
    bool) const {
    throw LogicalError("Component::trackBunch()", 
        "Called for component \"" + getName() + "\".");
}

void Component::trackMap(
    FVps<double, 6>&, 
    const PartData&, 
    bool, 
    bool) const {
    throw LogicalError("Component::trackMap()", 
        "Called for component \"" + getName() + "\".");
}

void Component::goOnline(const double&) {
    online_m = true;
}

void Component::goOffline() {
    online_m = false;
}

bool Component::Online() {
    return online_m;
}

ElementType Component::getType() const {
    return ElementType::ANY;
}


/* ============================== Apply Functions =========================== */
/**
 * @brief Apply to all particles. Kernel launch moved inside the function. 
 * 
 * @returns true if particle is out-of-bounds (lost), false otherwise
 */
bool Component::apply(){
    return false;
}

/**
 * @brief Apply to particle i
 * 
 * @param i Particle index
 * @param t Time
 * @param E Electric Field
 * @param B Magnetic Field
 * 
 * @returns true if particle is out-of-bounds (lost), false otherwise
 */
bool Component::apply(
    const size_t& i, 
    const double&, 
    Vector_t<double, 3>&, 
    Vector_t<double, 3>&) {
    /* 
    const Vector_t<double, 3>& R = RefPartBunch_m->R[i];
    if (R(2) >= 0.0 && R(2) < getElementLength()) {
        if (!isInsideTransverse(R))
            return true;
    }
    */
    return false;
}

/**
 * @brief Apply to particle with position R and momentum P
 * 
 * @param R Position 
 * @param P Momentum
 * @param t Time
 * @param E Electric Field
 * @param B Magnetic Field
 * 
 * @returns true if particle is out-of-bounds (lost), false otherwise
 */
bool Component::apply(
    const Vector_t<double, 3>& R, 
    const Vector_t<double, 3>& /*P*/, 
    const double& /*t*/,
    Vector_t<double, 3>& /*E*/, 
    Vector_t<double, 3>& /*B*/) {
    if (R(2) >= 0.0 && R(2) < getElementLength()) {
        if (!isInsideTransverse(R))
            return true;
    }
    return false;
}

/**
 * @brief Apply to reference particle with position R and momemtum P
 * 
 * @param R Position 
 * @param P Momentum
 * @param t Time
 * @param E Electric Field
 * @param B Magnetic Field
 * 
 * @returns true if particle is out-of-bounds (lost), false otherwise
 */
bool Component::applyToReferenceParticle(
    const Vector_t<double, 3>& R, 
    const Vector_t<double, 3>& /*P*/, 
    const double& /*t*/,
    Vector_t<double, 3>& /*E*/, 
    Vector_t<double, 3>& /*B*/) {
    if (R(2) >= 0.0 && R(2) < getElementLength()) {
        if (!isInsideTransverse(R))
            return true;
    }
    return false;
}
/* ========================================================================== */