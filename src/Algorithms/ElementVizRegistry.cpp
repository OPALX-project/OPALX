#include "Algorithms/ElementVizRegistry.h"

#include "AbsBeamline/ElementBase.h"
#include "AbsBeamline/ElementTypeTable.h"
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
            if (!viz::carriesField(type))
                continue;

            ElementField ef;
            ef.name    = element->getName();
            ef.element = element;
            ef.E.initialize(domain.mesh(), domain.layout(), nghost);
            ef.E = 0.0;
            ef.B.initialize(domain.mesh(), domain.layout(), nghost);
            ef.B = 0.0;

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
            // Zero both fields before accumulating.
            ef.E = 0.0;
            ef.B = 0.0;

            const auto& element = ef.element;
            if (!element)
                continue;

            // Transform from reference frame to element-local CS.
            const CoordinateSystemTrafo refToLocalCSTrafo =
                    (beamline.getMisalignment(element)
                     * (beamline.getCSTrafoLab2Local(element) * toLab));
            const CoordinateSystemTrafo localToRefCSTrafo = refToLocalCSTrafo.inverted();

            auto Eview = ef.E.getView();
            auto Bview = ef.B.getView();
            auto Ehost = Kokkos::create_mirror_view(Eview);
            auto Bhost = Kokkos::create_mirror_view(Bview);

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
                        const Vector_t<double, 3> refB = localToRefCSTrafo.rotateTo(localB);
                        Ehost(i, j, k) = refE;
                        Bhost(i, j, k) = refB;
                    }
                }
            }

            Kokkos::deep_copy(Eview, Ehost);
            Kokkos::deep_copy(Bview, Bhost);
        }
    }

}  // namespace opalx
