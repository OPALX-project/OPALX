#ifndef OPALX_CATALYST_ELEMENT_VIZ_REGISTRY_H
#define OPALX_CATALYST_ELEMENT_VIZ_REGISTRY_H

#include "OPALTypes.h"
#include "PartBunch/PartBunchFwd.h"
#include "PartBunch/ParticleContainerTypes.h"

#include <memory>
#include <string>
#include <vector>

class OpalBeamline;
class ElementBase;

namespace ippl { class CatalystAdaptor; }

namespace opalx {

    /**
     * @brief Owns external field arrays on the bunch mesh and samples
     *        beamline element fields onto them each step.
     *
     * When ExternalFieldsOnBunchMesh is active, one ippl::Field<Vector_t, 3>
     * is created per element with a field.  Each field is registered as
     * an array on the "mesh" block of the particle_bunch_N multimesh.
     * Every step, sampleDynamic() iterates the active elements (those near
     * the bunch) and calls element->apply(R, P, t, E, B) at each local
     * grid point, zeroing contributions outside the element bounds.
     */
    class ElementVizRegistry {
    public:
        using Mesh_t       = ippl::UniformCartesian<double, 3>;
        using Layout_t     = ippl::FieldLayout<3>;
        using Centering_t  = Mesh_t::DefaultCentering;
        using VectorField  = ippl::Field<ippl::Vector<double, 3>, 3, Mesh_t, Centering_t>;

        /**
         * @brief One-time creation of external field arrays on the bunch mesh.
         *
         * Iterates the beamline, creates a VectorField on the bunch mesh
         * for each element with a field.  The caller is responsible for
         * registering each field as a Catalyst array via addMeshBlock().
         *
         * @param beamline The OPAL beamline (provides elements + transforms).
         * @param bunch The PartBunch (provides mesh + layout).
         */
        void createFields(OpalBeamline& beamline, PartBunch_t& bunch);

        /**
         * @brief Iterate over all registered element fields.
         *
         * @param fn Callback receiving (name, field) for each element.
         */
        template <typename Fn>
        void forEachField(Fn&& fn) {
            for (auto& ef : elementFields_m) {
                fn("E_ext_" + ef.name, ef.E);
                fn("B_ext_" + ef.name, ef.B);
            }
        }

        /**
         * @brief Per-step sampling of external fields onto the bunch mesh.
         *
         * Iterates active elements (those whose bounds intersect the bunch
         * bounding box) and calls element->apply(R, P, t, E, B) at each
         * local grid point.  Fields are zeroed before accumulation, so
         * inactive elements contribute zero.
         *
         * @param beamline The OPAL beamline (provides element transforms).
         * @param bunch The PartBunch (provides mesh, time, reference).
         * @param bunchBBoxMin Bounding box minimum in reference frame.
         * @param bunchBBoxMax Bounding box maximum in reference frame.
         */
        void sampleDynamic(
                OpalBeamline& beamline, PartBunch_t& bunch,
                const Vector_t<double, 3>& bunchBBoxMin,
                const Vector_t<double, 3>& bunchBBoxMax);

        /// @return Number of registered element field arrays.
        size_t size() const { return elementFields_m.size(); }

    private:
        /// Owned field on the bunch mesh, one per element with a field.
        struct ElementField {
            std::string name;
            std::shared_ptr<ElementBase> element;
            VectorField E;  ///< Electric field sampled on bunch mesh
            VectorField B;  ///< Magnetic field sampled on bunch mesh
        };

        std::vector<ElementField> elementFields_m;
    };

}  // namespace opalx

#endif  // OPALX_CATALYST_ELEMENT_VIZ_REGISTRY_H
