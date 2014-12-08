#ifndef OPAL_SECTION_H
#define OPAL_SECTION_H

#include <vector>
#include <memory>

#include "Solvers/SurfacePhysicsHandler.hh"
#include "AbsBeamline/Component.h"

class WakeFunction;
class SurfacePhysicsHandler;
class BoundaryGeometry;

typedef std::vector<std::shared_ptr<Component> > CompVec;

class OpalSection {
public:
    OpalSection(const CompVec &, const double &, const double &);
    ~OpalSection();

    double getStart(const double &, const double &) const;
    void setStart(const double &);
    double getEnd(const double &, const double &) const;
    double getEnd() const {return end_m;}
    void setEnd(const double &);

    void setOrientation(const Vector_t &);
    double getOrientation(const int &) const;
    const Vector_t &getOrientation() const;

    void setStatus(const bool &);
    const bool &getStatus() const;

    WakeFunction *getWakeFunction();
    std::shared_ptr<const ElementBase> getWakeFunctionOwner();

    SurfacePhysicsHandler *getSurfacePhysicsHandler();
    BoundaryGeometry *getBoundaryGeometry();

    const bool &doesBend() const;
    const bool &hasWake() const;
    const bool &hasBoundaryGeometry() const;
    const bool &hasSurfacePhysics() const;

    void push_back(std::shared_ptr<Component>);
    bool find(const std::shared_ptr<Component>) const;
    CompVec &getElements();

    void print(Inform &) const;

    void previous_is_glued();
    void glue_to(OpalSection *);
    bool is_glued_to(const OpalSection *) const;

    static bool SortAsc(const OpalSection &sle1, const OpalSection &sle2) {
        return (sle1.start_m < sle2.start_m);
    }

    struct OrientationCache {
        double u_factor;
        double v_factor;
    };

    OrientationCache getStartCache() const {
        return StartCache_m;
    }

    OrientationCache getEndCache() const {
        return EndCache_m;
    }

private:
    CompVec elements_m;
    double start_m;
    double end_m;
    bool bends_m;
    bool has_wake_m;
    bool has_boundarygeometry_m;
    bool has_surface_physics_m;
    bool is_live_m;
    WakeFunction *wakefunction_m;
    std::shared_ptr<const ElementBase> wakeFunctionOwner_m;
    SurfacePhysicsHandler *sphys_handler_m;
    BoundaryGeometry *boundarygeometry_m;

    Vector_t orientation_m;
    double exit_face_angle_m;
    bool previous_is_glued_m;
    OpalSection *glued_to_m;

    // Cache computation for calls to getStart and getEnd, because all these
    // trigonometric function calls for every particle take too much time.
    // If certain member variables are changed, the cache needs to be rebuilt.
    // The cache is managed by this class, and its existence should be transparent for callers from
    // outside (except the speed gains, of course :-)). See also functions
    // updateGetStartCache, updateGetEndCache.
    OrientationCache StartCache_m;
    OrientationCache EndCache_m;

    void updateStartCache();

    void updateEndCache();
};

typedef std::vector<OpalSection> SectionList;

inline double OpalSection::getStart(const double &u, const double &v) const {
    return start_m + StartCache_m.u_factor * u + StartCache_m.v_factor * v;
}

inline void OpalSection::setStart(const double &start) {
    start_m = start;
}

inline double OpalSection::getEnd(const double &u, const double &v) const {
    return end_m + EndCache_m.u_factor * u + EndCache_m.v_factor * v;
}

inline void OpalSection::setEnd(const double &end) {
    end_m = end;
}

inline double OpalSection::getOrientation(const int &i) const {
    if(i == 0 || i == 1) {
        return orientation_m(i);
    } else {
        return 0.;
    }
}

inline const Vector_t &OpalSection::getOrientation() const {
    return orientation_m;
}

inline void OpalSection::setStatus(const bool &status) {
    is_live_m = status;
}

inline const bool &OpalSection::getStatus() const {
    return is_live_m;
}

inline const bool &OpalSection::hasWake() const {
    return has_wake_m;
}

inline const bool &OpalSection::hasBoundaryGeometry() const {
    return has_boundarygeometry_m;
}

inline WakeFunction *OpalSection::getWakeFunction() {
    return wakefunction_m;
}

inline std::shared_ptr<const ElementBase> OpalSection::getWakeFunctionOwner() {
    return wakeFunctionOwner_m;
}

inline BoundaryGeometry *OpalSection::getBoundaryGeometry() {
    return boundarygeometry_m;
}

inline const bool &OpalSection::hasSurfacePhysics() const {
    return has_surface_physics_m;
}

inline SurfacePhysicsHandler *OpalSection::getSurfacePhysicsHandler() {
    return sphys_handler_m;
}

inline const bool &OpalSection::doesBend() const {
    return bends_m;
}

inline void OpalSection::push_back(std::shared_ptr<Component> cmp) {
    elements_m.push_back(cmp);
}

inline CompVec &OpalSection::getElements() {
    return elements_m;
}

inline void OpalSection::previous_is_glued() {
    previous_is_glued_m = true;
    updateStartCache();
}

inline void OpalSection::glue_to(OpalSection *sec) {
    glued_to_m = sec;
    updateEndCache();
}

inline bool OpalSection::is_glued_to(const OpalSection *sec) const {
    return sec == glued_to_m;
}

#endif //OPAL_SECTION_H