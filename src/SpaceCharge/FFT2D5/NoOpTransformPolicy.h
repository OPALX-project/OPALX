/** @file NoOpTransformPolicy.h @brief Conservative FFT2D5 callback-frame filtering. */
#ifndef OPALX_FFT2D5_NO_OP_TRANSFORM_POLICY_H
#define OPALX_FFT2D5_NO_OP_TRANSFORM_POLICY_H

#include <typeinfo>

#include "Algorithms/CoordinateSystemTrafo.h"
#include "BeamlineCore/DriftRep.h"
#include "BeamlineCore/MarkerRep.h"

namespace opalx::spacecharge::fft2d5_detail {
    enum class FrameOperation { All, Fields, Monitors, Apertures };

    /** The caller supplies the existing MPI-agreed coasting eligibility. Unknown concrete
     * field types remain active: a subclass may override apply() even if its type is DRIFT.
     */
    inline bool callbackNeedsFrame(
            const ElementBase& element, FrameOperation operation, bool spatialRing, bool enabled) {
        if (!enabled) return true;
        switch (operation) {
            case FrameOperation::Monitors:
                return element.getType() == ElementType::MONITOR;
            case FrameOperation::Apertures:
                // This flag's false value is the virtual aperture callback's no-op contract.
                return element.getFlagDeleteOnTransverseExit();
            case FrameOperation::Fields:
                if (spatialRing && element.getType() == ElementType::MONITOR) return false;
                return typeid(element) != typeid(Drift) && typeid(element) != typeid(DriftRep)
                       && typeid(element) != typeid(Marker) && typeid(element) != typeid(MarkerRep);
            case FrameOperation::All:
                return true;
        }
        return true;
    }

    /** Run an unchanged callback in its element frame, unless it is a proven no-op.
     * Selection, bounds queries, occurrence order and frame composition stay with the tracker.
     */
    template <class Container, class Callback>
    bool inElementFrame(
            const ElementBase& element, Container& particles, const CoordinateSystemTrafo& toLocal,
            FrameOperation operation, bool spatialRing, bool enabled, Callback&& callback) {
        if (!callbackNeedsFrame(element, operation, spatialRing, enabled)) return false;
        const auto fromLocal = toLocal.inverted();
        particles.transformBunch(toLocal);
        callback();
        particles.transformBunch(fromLocal);
        return true;
    }
}  // namespace opalx::spacecharge::fft2d5_detail

#endif
