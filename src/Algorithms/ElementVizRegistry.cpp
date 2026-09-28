#include "Algorithms/ElementVizRegistry.h"

#include "AbsBeamline/ElementBase.h"
#include "Algorithms/Matrix.h"
#include "Elements/OpalBeamline.h"
#include "PartBunch/PartBunch.h"
#include "Utilities/OpalException.h"

#include <Kokkos_Core.hpp>

namespace opalx {

    void ElementVizRegistry::createFields(OpalBeamline& beamline, PartBunch_t& bunch) {
        auto& domain     = bunch.cartesianDomain();
        const int nghost = 1;

        const auto elements = beamline.getElements();
        for (const auto& element : elements) {
            if (!element)
                continue;

            const auto type = element->getType();
            // Skip elements that carry no field.
            if (type == ElementType::DRIFT || type == ElementType::MARKER
                || type == ElementType::MONITOR || type == ElementType::COLLIMATOR
                || type == ElementType::LASER || type == ElementType::PROBE
                || type == ElementType::VACUUM || type == ElementType::SOURCE
                || type == ElementType::ANY) {
                continue;
            }

            ElementField ef;
            ef.name    = "E_ext_" + element->getName();
            ef.element = element;
            ef.field.initialize(domain.mesh(), domain.layout(), nghost);
            ef.field = 0.0;

            elementFields_m.push_back(std::move(ef));
        }
    }

    void ElementVizRegistry::sampleDynamic(
            OpalBeamline& beamline, PartBunch_t& bunch,
            const Vector_t<double, 3>& bunchBBoxMin,
            const Vector_t<double, 3>& bunchBBoxMax) {

        if (elementFields_m.empty())
            return;

        const auto& domain = bunch.cartesianDomain();
        const auto& mesh   = domain.mesh();
        const auto& layout = domain.layout();
        const auto origin  = mesh.getOrigin();
        const auto h       = mesh.getMeshSpacing();

        const auto localND = layout.getLocalNDIndex();
        const int nghost   = 1;

        auto pc = bunch.getParticleContainer(0);
        if (!pc || !pc->getReference())
            return;

        const auto& toLab = pc->getToLabTrafo();
        const double t    = bunch.getT();

        for (auto& ef : elementFields_m) {
            // Zero the field before accumulating.
            ef.field = 0.0;

            const auto& element = ef.element;
            if (!element)
                continue;

            // Skip elements whose longitudinal bounds don't overlap the bunch.
            // Element bounds are in element-local coordinates; we compare
            // against the bunch bounding box in the reference frame.
            // This is a coarse cull — the per-point apply() does exact bounds.
            double zBegin = 0.0, zEnd = 0.0;
            element->getFieldExtent(zBegin, zEnd);
            const double elemZMin = zBegin;
            const double elemZMax = zEnd;
            if (elemZMax < bunchBBoxMin[2] || elemZMin > bunchBBoxMax[2])
                continue;

            // Transform from reference frame to element-local CS.
            const CoordinateSystemTrafo refToLocalCSTrafo =
                    (beamline.getMisalignment(element)
                     * (beamline.getCSTrafoLab2Local(element) * toLab));
            const CoordinateSystemTrafo localToRefCSTrafo = refToLocalCSTrafo.inverted();

            auto view = ef.field.getView();
            auto host = Kokkos::create_mirror_view(view);

            for (int i = localND[0].first() + nghost; i <= localND[0].last() + nghost; ++i) {
                const double x = (i - nghost) * h[0] + origin[0];
                for (int j = localND[1].first() + nghost; j <= localND[1].last() + nghost; ++j) {
                    const double y = (j - nghost) * h[1] + origin[1];
                    for (int k = localND[2].first() + nghost; k <= localND[2].last() + nghost;
                         ++k) {
                        const double z = (k - nghost) * h[2] + origin[2];

                        const Vector_t<double, 3> refPoint(x, y, z);
                        const Vector_t<double, 3> localR = refToLocalCSTrafo.transformTo(refPoint);

                        Vector_t<double, 3> localE(0.0);
                        Vector_t<double, 3> localB(0.0);
                        element->apply(localR, Vector_t<double, 3>(0.0), t, localE, localB);

                        const Vector_t<double, 3> refE = localToRefCSTrafo.rotateTo(localE);
                        host(i, j, k) = refE;
                    }
                }
            }

            Kokkos::deep_copy(view, host);
        }
    }

}  // namespace opalx
