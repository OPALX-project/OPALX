#ifndef OPAL_FIELD_H
#define OPAL_FIELD_H

#include <vector>
#include <list>
#include <memory>
#include "AbsBeamline/Component.h"

class OpalField {
public:
    OpalField(std::shared_ptr<Component>, const double &, const double &);
    ~OpalField();
    std::shared_ptr<Component> getElement();
    double getLength() const;
    const double &getStart() const;
    const double &getEnd() const;
    void setStart(const double & z);
    void setEnd(const double & z);
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
    std::shared_ptr<Component> element_m;
    double start_m;
    double end_m;
    bool is_on_m;
};

typedef std::list<OpalField> FieldList;

inline std::shared_ptr<Component> OpalField::getElement() {
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

inline void OpalField::setStart(const double & z) {
    start_m = z;
}

inline void OpalField::setEnd(const double & z) {
    end_m = z;
}
#endif // OPAL_FIELD_H
