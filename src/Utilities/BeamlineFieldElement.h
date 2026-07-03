#ifndef OPALX_BeamlineFieldElement_H
#define OPALX_BeamlineFieldElement_H

#include <list>
#include <memory>
#include "AbsBeamline/Component.h"
#include "Structure/BoundingBox.h"

/**
 * @brief Beamline component together with its active longitudinal field interval.
 *
 * BeamlineFieldElement is the tracking-side record used by OpalBeamline to keep
 * a component pointer, its start/end field positions, and its online/offline
 * state during tracking.
 */
class BeamlineFieldElement {
public:
    BeamlineFieldElement(std::shared_ptr<Component>, const double&, const double&);
    ~BeamlineFieldElement();
    std::shared_ptr<Component> getElement();
    std::shared_ptr<const Component> getElement() const;
    double getLength() const;
    const double& getStart() const;
    const double& getEnd() const;
    void setStart(const double& z);
    void setEnd(const double& z);
    const bool& isOn() const;
    void setOn(const double& kinematicEnergy);
    void setOff();

    static bool SortAsc(const BeamlineFieldElement& fle1, const BeamlineFieldElement& fle2) {
        return (fle1.start_m < fle2.start_m
                || (fle1.start_m == fle2.start_m
                    && fle1.element_m->getName() < fle2.element_m->getName()));
    }

    static bool ZeroLength(const BeamlineFieldElement& fle) { return (fle.getLength() < 1.e-6); }

    BoundingBox getBoundingBoxInLabCoords() const;

    unsigned int order_m;

private:
    std::shared_ptr<Component> element_m;
    double start_m;
    double end_m;
    bool is_on_m;
};

using FieldList = std::list<BeamlineFieldElement>;

inline std::shared_ptr<Component> BeamlineFieldElement::getElement() { return element_m; }

inline std::shared_ptr<const Component> BeamlineFieldElement::getElement() const {
    return element_m;
}

inline double BeamlineFieldElement::getLength() const { return end_m - start_m; }

inline const double& BeamlineFieldElement::getStart() const { return start_m; }

inline const double& BeamlineFieldElement::getEnd() const { return end_m; }

inline const bool& BeamlineFieldElement::isOn() const { return is_on_m; }

inline void BeamlineFieldElement::setStart(const double& z) { start_m = z; }

inline void BeamlineFieldElement::setEnd(const double& z) { end_m = z; }

inline BoundingBox BeamlineFieldElement::getBoundingBoxInLabCoords() const {
    return element_m->getBoundingBoxInLabCoords();
}
#endif  // OPALX_BeamlineFieldElement_H
