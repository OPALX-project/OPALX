/**
 * @file PicWorkspace.tpp
 * @brief Implements persistent Cartesian PIC workspace storage.
 */

#ifndef OPALX_SPACE_CHARGE_PIC_WORKSPACE_TPP
#define OPALX_SPACE_CHARGE_PIC_WORKSPACE_TPP

#include "Utilities/OpalException.h"

#include <limits>

namespace opalx::spacecharge {

    template <typename T, unsigned Dim>
    PicWorkspace<T, Dim>::PicWorkspace(
            Vector spacing, Vector lower, Vector upper, std::array<bool, Dim> decomposition,
            ippl::NDIndex<Dim> domain, Vector origin, bool allPeriodic)
        : spacing_m(spacing),
          lower_m(lower),
          upper_m(upper),
          decomposition_m(decomposition),
          mesh_m(domain, spacing, origin),
          layout_m(MPI_COMM_WORLD, domain, decomposition, allPeriodic),
          allPeriodic_m(allPeriodic),
          accumulatedElectricField_m(std::make_shared<VectorField>()),
          accumulatedMagneticField_m(std::make_shared<VectorField>()),
          flippedZSlabField_m(std::make_shared<VectorField>()) {}

    template <typename T, unsigned Dim>
    typename PicWorkspace<T, Dim>::Extents PicWorkspace<T, Dim>::layoutExtents() const {
        Extents extents{};
        const auto& domain = layout_m.getDomain();
        for (unsigned d = 0; d < Dim; ++d) {
            extents[d] = static_cast<std::size_t>(domain[d].length());
        }
        return extents;
    }

    template <typename T, unsigned Dim>
    bool PicWorkspace<T, Dim>::rebuildGlobalLayoutInPlace(
            const Extents& extents, const std::array<bool, Dim>& decomposition) {
        // Keep the resize policy even when the current extents already match. A later correction
        // transition may rebuild the layout with this stored decomposition.
        decomposition_m = decomposition;
        if (layoutExtents() == extents) {
            return false;
        }

        ippl::NDIndex<Dim> domain;
        for (unsigned d = 0; d < Dim; ++d) {
            if (extents[d] == 0
                || extents[d] > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
                throw OpalException(
                        "PicWorkspace::rebuildGlobalLayoutInPlace",
                        "Every PIC layout extent must fit in a positive IPPL Index.");
            }
            domain[d] = ippl::Index(static_cast<int>(extents[d]));
        }

        layout_m.initialize(domain, decomposition_m, allPeriodic_m);
        updateFieldLayoutsAfterLayoutChange();
        return true;
    }

    template <typename T, unsigned Dim>
    void PicWorkspace<T, Dim>::setGeometry(
            Vector lower, Vector upper, Vector spacing, Vector origin) {
        lower_m   = lower;
        upper_m   = upper;
        spacing_m = spacing;
        mesh_m.setMeshSpacing(spacing_m);
        mesh_m.setOrigin(origin);
    }

    template <typename T, unsigned Dim>
    void PicWorkspace<T, Dim>::initializeFields(std::string_view solverType) {
        Inform m("PicWorkspace::initializeFields");
        if (fieldsInitialized_m) {
            throw OpalException(
                    "PicWorkspace::initializeFields",
                    "The workspace fields are already initialized.");
        }

        m << level3 << "Mesh spacing = " << mesh_m.getMeshSpacing() << endl;
        m << level3 << "Origin       = " << mesh_m.getOrigin() << endl;
        m << level3 << "FL           = " << layout_m << endl;

        electricField_m.initialize(mesh_m, layout_m);
        chargeDensity_m.initialize(mesh_m, layout_m);
        potentialInitialized_m = solverType == "CG";
        if (potentialInitialized_m) {
            potential_m.initialize(mesh_m, layout_m);
        }

        accumulatedElectricField_m->initialize(mesh_m, layout_m);
        accumulatedMagneticField_m->initialize(mesh_m, layout_m);
        flippedZSlabField_m->initialize(mesh_m, layout_m, electricField_m.getNghost());
        fieldsInitialized_m = true;
        m << level3 << "Persistent PIC fields and scratch initialized." << endl;
    }

    template <typename T, unsigned Dim>
    void PicWorkspace<T, Dim>::updateFieldLayoutsAfterLayoutChange() {
        if (!fieldsInitialized_m) {
            throw OpalException(
                    "PicWorkspace::updateFieldLayoutsAfterLayoutChange",
                    "The workspace fields must be initialized before a layout refresh.");
        }

        electricField_m.updateLayout(layout_m);
        chargeDensity_m.updateLayout(layout_m);
        if (potentialInitialized_m) {
            potential_m.updateLayout(layout_m);
            potential_m = 0.0;
            potential_m.setFieldBC(potential_m.getFieldBC());
        }
        accumulatedElectricField_m->updateLayout(layout_m);
        accumulatedMagneticField_m->updateLayout(layout_m);
        flippedZSlabField_m->updateLayout(layout_m);
    }

    template <typename T, unsigned Dim>
    typename PicWorkspace<T, Dim>::VectorField& PicWorkspace<T, Dim>::mirrorScratchFor(
            const VectorField& source) {
        const bool incompatible = &flippedZSlabField_m->getLayout() != &source.getLayout()
                                  || flippedZSlabField_m->getNghost() != source.getNghost();
        if (incompatible) {
            throw OpalException(
                    "PicWorkspace::mirrorScratchFor",
                    "Persistent mirror scratch does not match the source field layout.");
        }
        return *flippedZSlabField_m;
    }

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_PIC_WORKSPACE_TPP
