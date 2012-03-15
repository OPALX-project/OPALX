#ifndef OPAL_SECTION_H
#define OPAL_SECTION_H

#include <vector>

#include "AbsBeamline/Component.h"
#include "Solvers/WakeFunction.hh"
#include "Algorithms/PartBunch.h"

typedef vector<Component*> CompVec;

class OpalSection {
public:
    OpalSection(CompVec&, const double&, const double&);
    ~OpalSection();

    double getStart(const double&, const double&) const;
    double getEnd(const double&, const double&) const;

    void setOrientation(const Vector_t&);
    double getOrientation(const int&) const;
    const Vector_t& getOrientation() const;

    void setStatus(const bool&);
    const bool& getStatus() const;

    WakeFunction* getWake();

    const bool& doesBend() const;
    const bool& hasWake() const;

    void push_back(Component*);
    bool find(const Component*) const;
    CompVec& getElements();

    void print(Inform&) const;


    static bool SortAsc(const OpalSection& sle1, const OpalSection& sle2) 
    {
        return (sle1.start_m < sle2.start_m);
    }
        
        
private:
    CompVec elements_m;
    double start_m;
    double end_m;
    bool bends_m;
    bool has_wake_m;
    bool is_live_m;
    WakeFunction *wake_m;
    Vector_t orientation_m;
};

inline double OpalSection::getOrientation(const int& i) const
{
    if (i == 0 || i == 1) {
        return orientation_m(i);
    } else {
        return 0.;
    }
}

inline const Vector_t& OpalSection::getOrientation() const
{
    return orientation_m;
}

inline void OpalSection::setStatus(const bool& status)
{
    is_live_m = status;
}

inline const bool& OpalSection::getStatus() const
{
    return is_live_m;
}

inline const bool& OpalSection::hasWake() const
{
    return has_wake_m;
}

inline WakeFunction* OpalSection::getWake()
{
    return wake_m;
}

inline const bool& OpalSection::doesBend() const
{
    return bends_m;
}

inline void OpalSection::push_back(Component* cmp)
{
    elements_m.push_back(cmp);
}

inline CompVec& OpalSection::getElements()
{
    return elements_m;
}
#endif //OPAL_SECTION_H
