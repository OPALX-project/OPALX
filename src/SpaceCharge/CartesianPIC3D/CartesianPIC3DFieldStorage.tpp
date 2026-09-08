/**
 * @file CartesianPIC3DFieldStorage.tpp
 * @brief Implements persistent CartesianPIC3D field storage.
 */

#ifndef OPALX_SPACE_CHARGE_CARTESIAN_PIC3D_FIELD_STORAGE_TPP
#define OPALX_SPACE_CHARGE_CARTESIAN_PIC3D_FIELD_STORAGE_TPP

#include "Utilities/OpalException.h"

#include <Kokkos_Core.hpp>

namespace opalx::spacecharge {

    template <typename T, unsigned Dim>
    CartesianPIC3DFieldStorage<T, Dim>::CartesianPIC3DFieldStorage(Domain& domain)
        : domain_m(domain) {}

    template <typename T, unsigned Dim>
    void CartesianPIC3DFieldStorage<T, Dim>::initializeFields(PoissonSolverType solverType) {
        Inform m("CartesianPIC3DFieldStorage::initializeFields");
        if (fieldsInitialized_m) {
            throw OpalException(
                    "CartesianPIC3DFieldStorage::initializeFields",
                    "The CartesianPIC3D fields are already initialized.");
        }

        m << level3 << "Mesh spacing = " << mesh().getMeshSpacing() << endl;
        m << level3 << "Origin       = " << mesh().getOrigin() << endl;
        m << level3 << "FL           = " << layout() << endl;

        electricField_m.initialize(mesh(), layout());
        chargeDensity_m.initialize(mesh(), layout());
        potentialInitialized_m = solverType == PoissonSolverType::ConjugateGradient;
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
    void CartesianPIC3DFieldStorage<T, Dim>::updateFieldLayoutsAfterLayoutChange() {
        if (!fieldsInitialized_m) {
            throw OpalException(
                    "CartesianPIC3DFieldStorage::updateFieldLayoutsAfterLayoutChange",
                    "The CartesianPIC3D fields must be initialized before a layout refresh.");
        }

        // updateLayout() reallocates Kokkos views. Complete work using the previous field storage
        // before any of those device allocations can be released.
        Kokkos::fence();
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
    typename CartesianPIC3DFieldStorage<T, Dim>::VectorField&
    CartesianPIC3DFieldStorage<T, Dim>::mirrorScratchFor(const VectorField& source) {
        const bool incompatible = &flippedZSlabField_m.getLayout() != &source.getLayout()
                                  || flippedZSlabField_m.getNghost() != source.getNghost();
        if (incompatible) {
            throw OpalException(
                    "CartesianPIC3DFieldStorage::mirrorScratchFor",
                    "Persistent mirror scratch does not match the source field layout.");
        }
        return flippedZSlabField_m;
    }

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_CARTESIAN_PIC3D_FIELD_STORAGE_TPP
