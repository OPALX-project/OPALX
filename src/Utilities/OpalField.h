#ifndef OPAL_FIELD_H
#define OPAL_FIELD_H

#include "AbsBeamline/Component.h"

class OpalField {
public:
    OpalField(Component *, const double &, const double &);
    ~OpalField();
    Component *getElement();
    double getLength() const;
    const double &getStart() const;
    const double &getEnd() const;
    const bool &isOn() const;
    void setOn();
    void setOff();

    static bool SortAsc(const OpalField &fle1, const OpalField &fle2) {
        return (fle1.start_m < fle2.start_m);
    }

    static bool ZeroLength(const OpalField &fle) {
        return (fle.getLength() < 1.e-6);
    }


private:
    Component *element_m;
    double start_m;
    double end_m;
    bool is_on_m;
};


inline Component *OpalField::getElement() {
    return element_m;
}

inline double OpalField::getLength() const {
    return end_m - start_m;
}

inline const double &OpalField::getStart() const {
    return start_m;
}

inline const double &OpalField::getEnd() const {
    return end_m;
}

inline const bool &OpalField::isOn() const {
    return is_on_m;
}

#endif // OPAL_FIELD_H
