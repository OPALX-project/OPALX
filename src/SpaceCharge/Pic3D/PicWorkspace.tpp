/**
 * @file PicWorkspace.tpp
 * @brief Implements persistent Cartesian PIC workspace storage.
 */

#ifndef OPALX_SPACE_CHARGE_PIC3D_WORKSPACE_TPP
#define OPALX_SPACE_CHARGE_PIC3D_WORKSPACE_TPP

#include "Utilities/OpalException.h"

namespace opalx::spacecharge {

    template <typename T, unsigned Dim>
    PicWorkspace<T, Dim>::PicWorkspace(Domain& domain) : domain_m(domain) {}

    template <typename T, unsigned Dim>
    void PicWorkspace<T, Dim>::initializeFields(std::string_view solverType) {
        Inform m("PicWorkspace::initializeFields");
        if (fieldsInitialized_m) {
            throw OpalException(
                    "PicWorkspace::initializeFields",
                    "The workspace fields are already initialized.");
        }

        m << level3 << "Mesh spacing = " << mesh().getMeshSpacing() << endl;
        m << level3 << "Origin       = " << mesh().getOrigin() << endl;
        m << level3 << "FL           = " << layout() << endl;

        electricField_m.initialize(mesh(), layout());
        chargeDensity_m.initialize(mesh(), layout());
        potentialInitialized_m = solverType == "CG";
        if (potentialInitialized_m) {
            potential_m.initialize(mesh(), layout());
        }

        accumulatedElectricField_m.initialize(mesh(), layout());
        accumulatedMagneticField_m.initialize(mesh(), layout());
        flippedZSlabField_m.initialize(mesh(), layout(), electricField_m.getNghost());
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

        electricField_m.updateLayout(layout());
        chargeDensity_m.updateLayout(layout());
        if (potentialInitialized_m) {
            potential_m.updateLayout(layout());
            potential_m = 0.0;
            potential_m.setFieldBC(potential_m.getFieldBC());
        }
        accumulatedElectricField_m.updateLayout(layout());
        accumulatedMagneticField_m.updateLayout(layout());
        flippedZSlabField_m.updateLayout(layout());
    }

    template <typename T, unsigned Dim>
    typename PicWorkspace<T, Dim>::VectorField& PicWorkspace<T, Dim>::mirrorScratchFor(
            const VectorField& source) {
        const bool incompatible = &flippedZSlabField_m.getLayout() != &source.getLayout()
                                  || flippedZSlabField_m.getNghost() != source.getNghost();
        if (incompatible) {
            throw OpalException(
                    "PicWorkspace::mirrorScratchFor",
                    "Persistent mirror scratch does not match the source field layout.");
        }
        return flippedZSlabField_m;
    }

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_PIC3D_WORKSPACE_TPP
