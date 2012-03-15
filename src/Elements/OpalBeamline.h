#ifndef OPAL_BEAMLINE_H
#define OPAL_BEAMLINE_H

#include <list>
#include <limits>
#include <vector>

#include "Algorithms/Tracker.h"
#include "Algorithms/PartBunch.h"

#include "Beamlines/Beamline.h"
#include "AbsBeamline/AlignWrapper.h"
#include "AbsBeamline/BeamBeam.h"
#include "AbsBeamline/Collimator.h"
#include "AbsBeamline/Corrector.h"
#include "AbsBeamline/Diagnostic.h"
#include "AbsBeamline/Drift.h"
#include "AbsBeamline/ElementBase.h"
#include "AbsBeamline/Lambertson.h"
#include "AbsBeamline/Marker.h"
#include "AbsBeamline/Monitor.h"
#include "AbsBeamline/Multipole.h"
#include "AbsBeamline/RBend.h"
#include "AbsBeamline/RFCavity.h"
#include "AbsBeamline/TravelingWave.h"
#include "AbsBeamline/RFQuadrupole.h"
#include "AbsBeamline/SBend.h"
#include "AbsBeamline/Separator.h"
#include "AbsBeamline/Septum.h"
#include "AbsBeamline/Solenoid.h"
#include "AbsBeamline/ParallelPlate.h"

#include "Solvers/SurfacePhysicsHandler.hh"
#include "Structure/SurfacePhysics.h"
#include "Structure/BoundaryGeometry.h"

#include "Utilities/OpalSection.h"
#include "Utilities/OpalField.h"

class WakeFunction;
class BoundaryGeometry;

#define BEAMLINE_EOL  0x80000000   // end of line
#define BEAMLINE_OOB  0x40000000   // out of bounds
#define BEAMLINE_GEOM 0x30000000   // has geometry
#define BEAMLINE_WAKE 0x20000000   // has wake
#define BEAMLINE_BEND 0x10000000   // bends
#define BEAMLINE_SURFACEPHYSICS 0x08000000 // has surface physics


class OpalBeamline {

public:
    OpalBeamline();
    ~OpalBeamline();

    CompVec &getPredecessors(const Component *);
    CompVec &getSuccessors(const Component *);
    OpalSection &getSectionAt(const Vector_t &, long &);
    OpalSection &getSection(const unsigned int &);
    void getSectionIndexAt(const Vector_t &, long &) const;
    double getSectionStart(const long &) const;
    double getSectionEnd(const unsigned int &) const;
    double getSectionEnd(const Vector_t &, long);

    double getStart(const Vector_t &) const;
    double getEnd(const Vector_t &) const;

    void setOrientation(const Vector_t &, const Vector_t &);
    void setOrientation(const Vector_t &, const unsigned int &);
    void updateOrientation(const Vector_t &, const Vector_t &, const double &, const double &);

    const Vector_t &getOrientation(const Vector_t &) const;
    const Vector_t &getOrientation(const long &) const;

    void resetStatus();
    void setStatus(const unsigned int &, const bool &);
    const bool &getStatus(const unsigned int &) const;

    void switchElements(const double &, const double &, const bool &nomonitors = false);
    void switchAllElements();

    void switchElementsOff(const double &, const string &eltype = "All");
    void switchElementsOff();

    WakeFunction *getWakeFunction(const unsigned int &);

    SurfacePhysicsHandler *getSurfacePhysicsHandler(const unsigned int &);

    BoundaryGeometry *getBoundaryGeometry(const unsigned int &);

    void getKFactors(const unsigned int &index, const Vector_t &pos, const long &sindex, const double &t, Vector_t &KR, Vector_t &KT);
    unsigned long getFieldAt(const unsigned int &, const Vector_t &, const long &, const double &, Vector_t &, Vector_t &);
    unsigned long getFieldAt(const Vector_t &, const Vector_t &, const double &, Vector_t &, Vector_t &);

    template<class T>
    void visit(const T &, Tracker &, PartBunch *);

    void prepareSections();
    void print(Inform &) const;

    bool section_is_glued_to(const long &i, const long &j) const;

    FieldList getElementByType(const std::string type);

    // need this for autophasing in case we have multiple tracks
    double calcBeamlineLenght();
    SectionList sections_m;

    void removeElement(const string &ElName);

private:
    FieldList elements_m;
    //SectionList sections_m;
    bool prepared_m;
    int *online_sections_m;
    bool *online_secs_m;

    static CompVec dummy_list_m;
    static OpalSection dummy_section_m;
};

inline void OpalBeamline::resetStatus() {
    for(unsigned int i = 0; i < sections_m.size(); ++ i) {
        sections_m[i].setStatus(false);
    }
}

inline void OpalBeamline::setStatus(const unsigned int &index, const bool &status) {
    if(index < sections_m.size()) {
        sections_m[index].setStatus(status);
    }
}

inline const bool &OpalBeamline::getStatus(const unsigned int &index) const {
    static const bool default_value = false;
    if(index < sections_m.size()) {
        return sections_m[index].getStatus();
    }
    return default_value;
}

inline WakeFunction *OpalBeamline::getWakeFunction(const unsigned int &index) {
    if(index < sections_m.size()) {
        return sections_m[index].getWakeFunction();
    }
    return 0;
}

inline BoundaryGeometry *OpalBeamline::getBoundaryGeometry(const unsigned int &index) {
    if(index < sections_m.size()) {
        return sections_m[index].getBoundaryGeometry();
    }
    return 0;
}

inline SurfacePhysicsHandler *OpalBeamline::getSurfacePhysicsHandler(const unsigned int &index) {
    if(index < sections_m.size()) {
        return sections_m[index].getSurfacePhysicsHandler();
    }
    return 0;
}

inline double OpalBeamline::getSectionStart(const long &index) const {
    if(index > -1 && (unsigned long) index < sections_m.size()) 
        return sections_m[index].getStart(0., 0.);

    return std::numeric_limits<double>::max();
}

inline double OpalBeamline::getSectionEnd(const unsigned int &index) const {
    if(index < sections_m.size()) 
        return sections_m[index].getEnd(0., 0.);

    return std::numeric_limits<double>::min();
}

inline double OpalBeamline::getSectionEnd(const Vector_t &pos, long sindex) {
    if(sindex > -1 && (unsigned int) sindex < sections_m.size()) {
        OpalSection &section = getSectionAt(pos, sindex);
        if(!(sindex & BEAMLINE_EOL || sindex < 0)) {
            return section.getEnd(pos(0), pos(2));
        }
    }
    return -1.0;
}

inline void OpalBeamline::setOrientation(const Vector_t &angle, const Vector_t &pos) {
    long index = 0;
    getSectionIndexAt(pos, index);
    if(index > -1 && (unsigned long) index < sections_m.size()) {
        sections_m[index].setOrientation(angle);
    }
}

inline void OpalBeamline::setOrientation(const Vector_t &angle, const unsigned int &index) {
    if(index < sections_m.size()) {
        sections_m[index].setOrientation(angle);
    }
}

inline void OpalBeamline::updateOrientation(const Vector_t &angle, const Vector_t &centroid, const double &smin, const double &smax) {

    // Determine which sections are online.
    int j = 0;
    for(int i = 0; i < 5 * Ippl::getNodes(); ++ i) {
        online_sections_m[i] = 0;
    }

    for(unsigned int i = 0; i < sections_m.size(); ++ i) {
        if(sections_m[i].getStatus()) {
            if(i > 0 && sections_m[i - 1].getStatus() && sections_m[i - 1].is_glued_to(&sections_m[i])) {
                continue;
            } else {
                online_sections_m[j + 5 * Ippl::myNode()] = i + 1;
                ++ j;
            }
        }
        if(j == 5) {
            break;
        }
    }

    reduce(online_sections_m, online_sections_m + 5 * Ippl::getNodes(), online_sections_m, OpAddAssign());
    for(int i = 0; i < 5 * Ippl::getNodes(); ++ i) {
        if(online_sections_m[i] > 0)
            online_secs_m[online_sections_m[i] - 1] = true;
    }

    // Update orientation of online sections.
    for(unsigned int i = 0; i < sections_m.size(); ++ i) {
        if(online_secs_m[i]) {
            online_secs_m[i] = false;
            Vector_t oldAngle = sections_m[i].getOrientation();

            // Rotate about z axis before changing section's Euler angles.
            Vector_t newAngle;
            double cosc = cos(oldAngle(2));
            double sinc = -sin(oldAngle(2));

            newAngle(0) = oldAngle(0) + cosc * angle(0) - sinc * angle(1);
            newAngle(1) = oldAngle(1) + sinc * angle(0) + cosc * angle(1);
            newAngle(2) = oldAngle(2);

            sections_m[i].setOrientation(newAngle);

            if(smin < sections_m[i].getStart(0.0, 0.0)) {
                double dist = sections_m[i].getStart(0.0, 0.0) - centroid(2);
                double newdist = dist / (1. - tan(oldAngle(0)) * angle(0));  // * cos(oldAngle(0)) / cos(newAngle(0));
                sections_m[i].setStart(centroid(2) + newdist);
            } else if(smax > sections_m[i].getEnd(0.0, 0.0)) {
                double dist = sections_m[i].getEnd(0.0, 0.0) - centroid(2);
                double newdist = dist / (1. - tan(oldAngle(0)) * angle(0)); // * cos(oldAngle(0)) / cos(newAngle(0));
                sections_m[i].setEnd(centroid(2) + newdist);
                if(i + 1 < sections_m.size() && sections_m[i].is_glued_to(&sections_m[i + 1])) {
                    sections_m[i + 1].setStart(centroid(2) + newdist);
                }
            }
        }
    }
}

inline const Vector_t &OpalBeamline::getOrientation(const Vector_t &pos) const {
    static const Vector_t dummy(0.0);
    long index = 0;
    getSectionIndexAt(pos, index);
    if(index > -1 && (unsigned long) index < sections_m.size()) {
        return sections_m[index].getOrientation();
    } else {
        return dummy;
    }
}

inline const Vector_t &OpalBeamline::getOrientation(const long &i) const {
    static const Vector_t dummy(0.0);
    if(i > -1 && (unsigned long) i < sections_m.size()) {
        return sections_m[i].getOrientation();
    } else {
        return dummy;
    }
}

inline bool OpalBeamline::section_is_glued_to(const long &i, const long &j) const {
    if(i > -1 && i < j && (unsigned long) j < sections_m.size()) {
        for(int k = i; k < j; ++ k) {
            if(!sections_m[k].is_glued_to(&sections_m[k + 1])) {
                return false;
            }
        }
        return true;
    }
    return false;
}

template<class T> inline
void OpalBeamline::visit(const T &element, Tracker &, PartBunch *bunch) {
    Inform msg("OPAL ");
    T *elptr = dynamic_cast<T *>(element.clone()->removeWrappers());
    if(!elptr->hasAttribute("ELEMEDGE")) {
        msg << elptr->getType() << ": no position of the element given!" << endl;
        return;
    }
    //     if (elptr->hasWake()) {
    //         const Wake *wf = elptr->getWake();
    //     }

    double startField = elptr->getAttribute("ELEMEDGE");
    double endField;
    elptr->initialise(bunch, startField, endField, 1.0);
    elements_m.push_back(OpalField(elptr, startField, endField));
}

template<> inline
void OpalBeamline::visit<AlignWrapper>(const AlignWrapper &wrap, Tracker &aTracker, PartBunch *) {
    if(wrap.getType() == "beamline") {
        Beamline *bl = dynamic_cast<Beamline *>(wrap.getElement());
        bl->iterate(aTracker, false);
    } else {
        wrap.getElement()->accept(aTracker);
    }
}

template<> inline
void OpalBeamline::visit<Corrector>(const Corrector &element, Tracker &, PartBunch *) {
    Inform msg("OPAL ");
    msg << element.getType() << " not implemented yet!" << endl;
}

template<> inline
void OpalBeamline::visit<BeamBeam>(const BeamBeam &element, Tracker &, PartBunch *) {
    Inform msg("OPAL ");
    msg << element.getType() << " not implemented yet!" << endl;
}

template<> inline
void OpalBeamline::visit<Diagnostic>(const Diagnostic &element, Tracker &, PartBunch *) {
    Inform msg("OPAL ");
    msg << element.getType() << " not implemented yet!" << endl;
}

template<> inline
void OpalBeamline::visit<Lambertson>(const Lambertson &element, Tracker &, PartBunch *) {
    Inform msg("OPAL ");
    msg << element.getType() << " not implemented yet!" << endl;
}

template<> inline
void OpalBeamline::visit<Marker>(const Marker &element, Tracker &, PartBunch *) {
    WARNMSG(element.getType() << " not implemented yet!" << endl);
}

template<> inline
void OpalBeamline::visit<RFQuadrupole>(const RFQuadrupole &element, Tracker &, PartBunch *) {
    Inform msg("OPAL ");
    msg << element.getType() << " not implemented yet!" << endl;
}

template<> inline
void OpalBeamline::visit<Separator>(const Separator &element, Tracker &, PartBunch *) {
    Inform msg("OPAL ");
    msg << element.getType() << " not implemented yet!" << endl;
}

template<> inline
void OpalBeamline::visit<Septum>(const Septum &element, Tracker &, PartBunch *) {
    Inform msg("OPAL ");
    msg << element.getType() << " not implemented yet!" << endl;
}

inline
void OpalBeamline::removeElement(const string &ElName) {
    for(FieldList::iterator flit = elements_m.begin(); flit != elements_m.end(); ++ flit) {
        if(flit->getElement()->getName() == ElName) {
            flit->setStart(-flit->getEnd());
            flit->setEnd(-flit->getEnd());
        }
    }
}

#endif // OPAL_BEAMLINE_H
