#ifndef OPALX_FIELD_OPS_HPP
#define OPALX_FIELD_OPS_HPP

#include "PartBunch/FieldContainer.hpp"

namespace opalx::fieldops {

    template <unsigned Dim>
    void setScalarField(Field_t<Dim>& field, double value) {
        auto view = field.getView();
        Kokkos::deep_copy(view, value);
    }

    template <unsigned Dim>
    void scaleAndShiftScalarField(Field_t<Dim>& field, double scale, double shift) {
        auto view = field.getView();

        ippl::parallel_for(
                "opalx::fieldops::scaleAndShiftScalarField", field.getFieldRangePolicy(),
                KOKKOS_LAMBDA(const typename ippl::RangePolicy<Dim>::index_array_type& idx) {
                    apply(view, idx) = apply(view, idx) * scale + shift;
                });
    }

    template <typename T, unsigned Dim>
    void setVectorField(VField_t<T, Dim>& field, const Vector_t<T, Dim>& value) {
        auto view = field.getView();
        Kokkos::deep_copy(view, value);
    }

}  // namespace opalx::fieldops

#endif  // OPALX_FIELD_OPS_HPP
