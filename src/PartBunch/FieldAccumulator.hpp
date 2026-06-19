#ifndef OPALX_FIELD_ACCUMULATOR_HPP
#define OPALX_FIELD_ACCUMULATOR_HPP

#include <memory>

#include "Ippl.h"

#include "PartBunch/FieldContainer.hpp"
#include "PartBunch/FieldMirror.hpp"
#include "PartBunch/FieldOps.hpp"
#include "Physics/Physics.h"
#include "Utilities/OpalException.h"

template <typename T, unsigned Dim>
class FieldAccumulator {
public:
    void bind(
            std::shared_ptr<VField_t<T, Dim>> electricField,
            std::shared_ptr<VField_t<T, Dim>> magneticField) {
        if (!electricField) {
            throw OpalException(
                    "FieldAccumulator::bind",
                    "Temporary E field (Etmp) is not initialized.");
        }
        if (!magneticField) {
            throw OpalException(
                    "FieldAccumulator::bind",
                    "Temporary B field (Btmp) is not initialized.");
        }

        electricField_m = electricField;
        magneticField_m = magneticField;
    }

    void clear() {
        requireBound("FieldAccumulator::clear");
        opalx::fieldops::setVectorField(*electricField_m, Vector_t<T, Dim>(0.0));
        opalx::fieldops::setVectorField(*magneticField_m, Vector_t<T, Dim>(0.0));
    }

    void accumulate(
            const VField_t<T, Dim>& ePrime, double gammaBin,
            const Vector_t<double, Dim>& pmean, double bFieldSign = 1.0, int flipAxis = -1) {
        requireBound("FieldAccumulator::accumulate");

        const double invGamma         = (gammaBin > 0.0) ? (1.0 / gammaBin) : 0.0;
        const Vector_t<double, Dim> v = Physics::c * pmean * invGamma;
        const double vNorm            = Kokkos::sqrt(v.dot(v));
        const Vector_t<double, Dim> w = (vNorm > 0.0) ? (v / vNorm) : Vector_t<double, Dim>(0.0);
        const double gammaMinusOne    = gammaBin - 1.0;
        const double gammaOverCSq     = gammaBin / (Physics::c * Physics::c);

        auto ePrimeView = ePrime.getView();
        auto eTmpView   = electricField_m->getView();
        auto bTmpView   = magneticField_m->getView();

        if (flipAxis < 0) {
            ippl::parallel_for(
                    "FieldAccumulator::accumulate", ePrime.getFieldRangePolicy(),
                    KOKKOS_LAMBDA(const ippl::RangePolicy<Dim>::index_array_type& idx) {
                        Vector_t<T, Dim> eFrame = apply(ePrimeView, idx);
                        const T eDotW           = eFrame.dot(w);
                        Vector_t<T, Dim> eLab =
                                gammaBin * eFrame - gammaMinusOne * eDotW * w;
                        Vector_t<T, Dim> bLab =
                                bFieldSign * gammaOverCSq * cross(v, eFrame);
                        Vector_t<T, Dim> eTotal = apply(eTmpView, idx);
                        Vector_t<T, Dim> bTotal = apply(bTmpView, idx);
                        eTotal += eLab;
                        bTotal += bLab;
                        apply(eTmpView, idx) = eTotal;
                        apply(bTmpView, idx) = bTotal;
                    });
            return;
        }

        const int flipAxisCap = flipAxis;
        if (flipAxisCap != static_cast<int>(Dim) - 1) {
            throw OpalException(
                    "FieldAccumulator::accumulate",
                    "flipAxis != Dim-1 not supported (only z-axis flip implemented).");
        }

        buildFlippedField(ePrime);
        auto flippedView = flippedField_m->getView();

        ippl::parallel_for(
                "FieldAccumulator::accumulate[flipped]", ePrime.getFieldRangePolicy(),
                KOKKOS_LAMBDA(const ippl::RangePolicy<Dim>::index_array_type& idx) {
                    Vector_t<T, Dim> eFrame = flippedView(idx[0], idx[1], idx[2]);

                    for (unsigned d = 0; d < Dim; ++d) {
                        if (static_cast<int>(d) != flipAxisCap) {
                            eFrame[d] = -eFrame[d];
                        }
                    }

                    const T eDotW = eFrame.dot(w);
                    Vector_t<T, Dim> eLab =
                            gammaBin * eFrame - gammaMinusOne * eDotW * w;
                    Vector_t<T, Dim> bLab = bFieldSign * gammaOverCSq * cross(v, eFrame);
                    Vector_t<T, Dim> eTotal = apply(eTmpView, idx);
                    Vector_t<T, Dim> bTotal = apply(bTmpView, idx);
                    eTotal += eLab;
                    bTotal += bLab;
                    apply(eTmpView, idx) = eTotal;
                    apply(bTmpView, idx) = bTotal;
                });
    }

    template <typename ParticleContainer>
    void gatherToParticles(ParticleContainer& pc) {
        requireBound("FieldAccumulator::gatherToParticles");
        gather(pc.E, *electricField_m, pc.R);
        gather(pc.B, *magneticField_m, pc.R);
    }

private:
    void requireBound(const char* where) const {
        if (!electricField_m || !magneticField_m) {
            throw OpalException(where, "Temporary field accumulator is not bound.");
        }
    }

    void buildFlippedField(const VField_t<T, Dim>& src) {
        auto& layout         = src.getLayout();
        auto& mesh           = src.get_mesh();
        const int srcNghost  = src.getNghost();
        const bool needsInit = !flippedField_m || &flippedField_m->getLayout() != &layout
                               || flippedField_m->getNghost() != srcNghost;

        if (!flippedField_m) {
            flippedField_m = std::make_shared<VField_t<T, Dim>>();
        }
        if (needsInit) {
            flippedField_m->initialize(mesh, layout, srcNghost);
        }

        opalx::detail::mirrorField(src, *flippedField_m, Dim - 1);
    }

    std::shared_ptr<VField_t<T, Dim>> electricField_m;
    std::shared_ptr<VField_t<T, Dim>> magneticField_m;
    std::shared_ptr<VField_t<T, Dim>> flippedField_m;
};

#endif  // OPALX_FIELD_ACCUMULATOR_HPP
