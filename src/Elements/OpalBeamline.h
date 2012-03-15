#ifndef OPAL_BEAMLINE_H
#define OPAL_BEAMLINE_H

#include <list>
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

#include "Solvers/WakeFunction.hh"
#include "Structure/Wake.h"

#include "Utilities/OpalSection.h"
#include "Utilities/OpalField.h"

#define BEAMLINE_EOL  0x80000000   // end of line
#define BEAMLINE_OOB  0x40000000   // out of bounds
#define BEAMLINE_WAKE 0x20000000   // has wake
#define BEAMLINE_BEND 0x10000000   // bends

typedef list<OpalField> FieldList;
typedef vector<OpalSection> SectionList;

class OpalBeamline {

public:
    OpalBeamline();
    ~OpalBeamline();
    
    CompVec& getPredecessors(const Component*);
    CompVec& getSuccessors(const Component*);
    OpalSection& getSectionAt(const Vector_t&, unsigned int&);
    OpalSection& getSection(const unsigned int&);
    void getSectionIndexAt(const Vector_t&, unsigned int&) const;
    double getSectionStart(const unsigned int&) const;
    double getSectionEnd(const unsigned int&) const;

    double getStart(const Vector_t&) const;
    double getEnd(const Vector_t&) const;

    void setOrientation(const Vector_t&, const Vector_t&);
    void setOrientation(const Vector_t&, const unsigned int&);
    void updateOrientation(const Vector_t&);

    const Vector_t& getOrientation(const Vector_t&) const;
    const Vector_t& getOrientation(const unsigned int&) const;

    void resetStatus();
    void setStatus(const unsigned int&, const bool&);
    const bool& getStatus(const unsigned int&) const;

    void switchElements(const double&, const double&);

    void switchElementsOff(const double&, const double&);
    void switchElementsOff();

    WakeFunction* getWake(const unsigned int&);

    const unsigned int& getFieldAt(const unsigned int&, const Vector_t&, const unsigned int&, const double&, Vector_t&, Vector_t&);
    const unsigned int& getFieldAt(const Vector_t&, const double&, Vector_t&, Vector_t&);

    template<class T>
    void visit(const T&, Tracker&, PartBunch*);

    void prepareSections();
    void print(Inform&) const;

private:
    FieldList elements_m;
    SectionList sections_m;
    bool prepared_m;

    static CompVec dummy_list_m;
    static OpalSection dummy_section_m;
};

inline void OpalBeamline::resetStatus()
{
    for (int i = 0; i < sections_m.size(); ++ i) {
        sections_m[i].setStatus(false);
    }
}

inline void OpalBeamline::setStatus(const unsigned int& index, const bool& status)
{
    if (index < sections_m.size()) {
        sections_m[index].setStatus(status);
    }
}

inline const bool& OpalBeamline::getStatus(const unsigned int& index) const
{
    if (index < sections_m.size()) {
        return sections_m[index].getStatus();
    }
}

inline WakeFunction* OpalBeamline::getWake(const unsigned int& index)
{
    if (index < sections_m.size()) {
        return sections_m[index].getWake();
    }
}

inline double OpalBeamline::getSectionStart(const unsigned int& index) const
{
    if (index < sections_m.size()) {
        return sections_m[index].getStart(0.,0.);
    } else {
        return numeric_limits<double>::max();
    }
}

inline double OpalBeamline::getSectionEnd(const unsigned int& index) const
{
    if (index < sections_m.size()) {
        return sections_m[index].getEnd(0.,0.);
    } else {
        return numeric_limits<double>::min();
    }
}

inline void OpalBeamline::setOrientation(const Vector_t& angle, const Vector_t& pos)
{
    unsigned int index = 0;
    getSectionIndexAt(pos,index);
    sections_m[index].setOrientation(angle);
}

inline void OpalBeamline::setOrientation(const Vector_t& angle, const unsigned int& index)
{
    sections_m[index].setOrientation(angle);
}

inline void OpalBeamline::updateOrientation(const Vector_t& angle)
{
    for (SectionList::iterator slit = sections_m.begin(); slit != sections_m.end(); ++ slit) {
        if ((*slit).getStatus()) {
            (*slit).setOrientation((*slit).getOrientation() + angle);
        }
    }
}

inline const Vector_t& OpalBeamline::getOrientation(const Vector_t& pos) const
{
    unsigned int index = 0;
    getSectionIndexAt(pos,index);
    return sections_m[index].getOrientation();
}

inline const Vector_t& OpalBeamline::getOrientation(const unsigned int& i) const
{
    static const Vector_t dummy(0.0);
    if (i < sections_m.size()) {
        return sections_m[i].getOrientation();
    } else {
        return dummy;
    }
}

template<class T> inline
void OpalBeamline::visit(const T& element, Tracker&, PartBunch* bunch)
{
    T *elptr = dynamic_cast<T*>(element.clone()->removeWrappers());
    if (!elptr->hasAttribute("ELEMEDGE")) {
        *gmsg << elptr->getType() << ": no position of the element given!" << endl;
        return;
    }

    double startField = elptr->getAttribute("ELEMEDGE");
    double endField;
    elptr->initialise(bunch, startField, endField, 1.0);
    elements_m.push_back(OpalField(elptr, startField, endField));
}

template<> inline
void OpalBeamline::visit<AlignWrapper>(const AlignWrapper& wrap, Tracker& aTracker, PartBunch*)
{
    if (wrap.getType() == "beamline") {
        Beamline *bl = dynamic_cast<Beamline*>(wrap.getElement());
        bl->iterate(aTracker, false);
    }
    else {
        wrap.getElement()->accept(aTracker);
    }
}

template<> inline
void OpalBeamline::visit<Corrector>(const Corrector& element, Tracker&, PartBunch*)
{
    *gmsg << element.getType() << " not implemented yet!" << endl;
}

template<> inline
void OpalBeamline::visit<BeamBeam>(const BeamBeam& element, Tracker&, PartBunch*)
{
    *gmsg << element.getType() << " not implemented yet!" << endl;
}

template<> inline
void OpalBeamline::visit<Diagnostic>(const Diagnostic& element, Tracker&, PartBunch*)
{
    *gmsg << element.getType() << " not implemented yet!" << endl;
}

template<> inline
void OpalBeamline::visit<Lambertson>(const Lambertson& element, Tracker&, PartBunch*)
{
    *gmsg << element.getType() << " not implemented yet!" << endl;
}

template<> inline
void OpalBeamline::visit<Marker>(const Marker& element, Tracker&, PartBunch*)
{
    *gmsg << element.getType() << " not implemented yet!" << endl;
}

template<> inline
void OpalBeamline::visit<RFQuadrupole>(const RFQuadrupole& element, Tracker&, PartBunch*)
{
    *gmsg << element.getType() << " not implemented yet!" << endl;
}

template<> inline
void OpalBeamline::visit<SBend>(const SBend& element, Tracker&, PartBunch*)
{
    *gmsg << element.getType() << " not implemented yet!" << endl;
}

template<> inline
void OpalBeamline::visit<Separator>(const Separator& element, Tracker&, PartBunch*)
{
    *gmsg << element.getType() << " not implemented yet!" << endl;
}

template<> inline
void OpalBeamline::visit<Septum>(const Septum& element, Tracker&, PartBunch*)
{
    *gmsg << element.getType() << " not implemented yet!" << endl;
}

#endif // OPAL_BEAMLINE_H
