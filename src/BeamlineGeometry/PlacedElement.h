#ifndef OPALX_PlacedElement_HH
#define OPALX_PlacedElement_HH

#include "BeamlineGeometry/ElementGeometry.h"
#include "BeamlineGeometry/Misalignment.h"
#include "BeamlineGeometry/PlacementPose.h"
#include "BeamlineGeometry/SupportPlacement.h"

class ElementBase;

/**
 * @class PlacedElement
 * @brief Geometric placement record for an element instance.
 *
 * The placed element combines:
 * - the nominal parent-to-body placement,
 * - the local nominal-to-actual correction,
 * - the canonical body-relative ports,
 * - and an optional support/frame offset.
 *
 * All transforms are represented with CoordinateSystemTrafo in the current
 * implementation. The actual body transform is therefore
 * \f[
 * T_{\mathrm{parent}\rightarrow\mathrm{actual}} =
 * T_{\mathrm{nominal}\rightarrow\mathrm{actual}}
 * \cdot
 * T_{\mathrm{parent}\rightarrow\mathrm{nominal}}.
 * \f]
 *
 * The bridge stage keeps both nominal and actual transforms available. This is
 * important because current OPALX runtime code still queries nominal placement
 * and applies misalignment separately in some paths.
 */
class PlacedElement {
public:
    explicit PlacedElement(
            const ElementBase* element, const PlacementPose& nominal = PlacementPose(),
            const Misalignment& correction  = Misalignment(),
            const ElementGeometry& geometry = ElementGeometry(),
            const SupportPlacement& support = SupportPlacement())
        : element_m(element),
          nominal_m(nominal),
          correction_m(correction),
          geometry_m(geometry),
          support_m(support) {}

    const ElementBase* getElement() const { return element_m; }
    const PlacementPose& getPlacementPose() const { return nominal_m; }
    const Misalignment& getMisalignment() const { return correction_m; }
    const ElementGeometry& getGeometry() const { return geometry_m; }
    const SupportPlacement& getSupportPlacement() const { return support_m; }

    CoordinateSystemTrafo getNominalBodyTransform() const { return nominal_m.getParentToNominal(); }

    CoordinateSystemTrafo getNominalEntryTransform() const {
        return geometry_m.getEntry().getBodyToPort() * getNominalBodyTransform();
    }

    CoordinateSystemTrafo getNominalBodyPortTransform() const {
        return geometry_m.getBody().getBodyToPort() * getNominalBodyTransform();
    }

    CoordinateSystemTrafo getNominalExitTransform() const {
        return geometry_m.getExit().getBodyToPort() * getNominalBodyTransform();
    }

    CoordinateSystemTrafo getNominalSupportTransform() const {
        return support_m.getBodyToSupport() * getNominalBodyTransform();
    }

    CoordinateSystemTrafo getActualBodyTransform() const {
        return correction_m.getNominalToActual() * nominal_m.getParentToNominal();
    }

    CoordinateSystemTrafo getEntryTransform() const {
        return geometry_m.getEntry().getBodyToPort() * getActualBodyTransform();
    }

    CoordinateSystemTrafo getBodyTransform() const {
        return geometry_m.getBody().getBodyToPort() * getActualBodyTransform();
    }

    CoordinateSystemTrafo getExitTransform() const {
        return geometry_m.getExit().getBodyToPort() * getActualBodyTransform();
    }

    CoordinateSystemTrafo getSupportTransform() const {
        return support_m.getBodyToSupport() * getActualBodyTransform();
    }

private:
    const ElementBase* element_m;
    PlacementPose nominal_m;
    Misalignment correction_m;
    ElementGeometry geometry_m;
    SupportPlacement support_m;
};

#endif
