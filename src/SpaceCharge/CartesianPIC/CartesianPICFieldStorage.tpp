/**
 * @file CartesianPICFieldStorage.tpp
 * @brief Implements persistent Cartesian PIC field storage.
 */

#ifndef OPALX_SPACE_CHARGE_CARTESIAN_PIC_FIELD_STORAGE_TPP
#define OPALX_SPACE_CHARGE_CARTESIAN_PIC_FIELD_STORAGE_TPP

#include "Utilities/OpalException.h"

#include <Kokkos_Core.hpp>

namespace opalx::spacecharge {

    template <typename T, unsigned Dim>
    CartesianPICFieldStorage<T, Dim>::CartesianPICFieldStorage(Domain& domain) : domain_m(domain) {}

    template <typename T, unsigned Dim>
    void CartesianPICFieldStorage<T, Dim>::initializeFields(std::string_view solverType) {
        Inform m("CartesianPICFieldStorage::initializeFields");
        if (fieldsInitialized_m) {
            throw OpalException(
                    "CartesianPICFieldStorage::initializeFields",
                    "The Cartesian PIC fields are already initialized.");
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
    void CartesianPICFieldStorage<T, Dim>::updateFieldLayoutsAfterLayoutChange() {
        if (!fieldsInitialized_m) {
            throw OpalException(
                    "CartesianPICFieldStorage::updateFieldLayoutsAfterLayoutChange",
                    "The Cartesian PIC fields must be initialized before a layout refresh.");
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
    typename CartesianPICFieldStorage<T, Dim>::VectorField&
    CartesianPICFieldStorage<T, Dim>::mirrorScratchFor(const VectorField& source) {
        const bool incompatible = &flippedZSlabField_m.getLayout() != &source.getLayout()
                                  || flippedZSlabField_m.getNghost() != source.getNghost();
        if (incompatible) {
            throw OpalException(
                    "CartesianPICFieldStorage::mirrorScratchFor",
                    "Persistent mirror scratch does not match the source field layout.");
        }
        return flippedZSlabField_m;
    }

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_CARTESIAN_PIC_FIELD_STORAGE_TPP
