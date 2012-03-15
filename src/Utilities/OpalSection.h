#ifndef OPAL_SECTION_H
#define OPAL_SECTION_H

#include <vector>

#include "AbsBeamline/Component.h"
#include "Algorithms/PartBunch.h"

class WakeFunction;
class SurfacePhysicsHandler;
class BoundaryGeometry;

typedef vector<Component *> CompVec;

class OpalSection {
public:
    OpalSection(CompVec &, const double &, const double &);
    ~OpalSection();

    double getStart(const double &, const double &) const;
    void setStart(const double &);
    double getEnd(const double &, const double &) const;
    double getEnd() {return end_m;}
    void setEnd(const double &);
    
    void setOrientation(const Vector_t &);
    double getOrientation(const int &) const;
    const Vector_t &getOrientation() const;

    void setStatus(const bool &);
    const bool &getStatus() const;

    WakeFunction *getWakeFunction();
    SurfacePhysicsHandler *getSurfacePhysicsHandler();
    BoundaryGeometry *getBoundaryGeometry();

    const bool &doesBend() const;
    const bool &hasWake() const;
    const bool &hasBoundaryGeometry() const;
    const bool &hasSurfacePhysics() const;

    void push_back(Component *);
    bool find(const Component *) const;
    CompVec &getElements();

    void print(Inform &) const;

    void previous_is_glued();
    void glue_to(OpalSection *);
    bool is_glued_to(const OpalSection *) const;

    static bool SortAsc(const OpalSection &sle1, const OpalSection &sle2) {
        return (sle1.start_m < sle2.start_m);
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
    SurfacePhysicsHandler *sphys_handler_m;
    BoundaryGeometry *boundarygeometry_m;

    Vector_t orientation_m;
    double exit_face_angle_m;
    bool previous_is_glued_m;
    OpalSection *glued_to_m;
};

inline void OpalSection::setStart(const double &start) {
    start_m = start;
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

inline void OpalSection::push_back(Component *cmp) {
    elements_m.push_back(cmp);
}

inline CompVec &OpalSection::getElements() {
    return elements_m;
}

inline void OpalSection::previous_is_glued() {
    previous_is_glued_m = true;
}

inline void OpalSection::glue_to(OpalSection *sec) {
    glued_to_m = sec;
}

inline bool OpalSection::is_glued_to(const OpalSection *sec) const {
    return sec == glued_to_m;
}

#endif //OPAL_SECTION_H
