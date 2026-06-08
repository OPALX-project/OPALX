/*
 *  Copyright (c) 2025, Jon Thompson
 *  All rights reserved.
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *  1. Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer.
 *  2. Redistributions in binary form must reproduce the above copyright notice,
 *     this list of conditions and the following disclaimer in the documentation
 *     and/or other materials provided with the distribution.
 *  3. Neither the name of STFC nor the names of its contributors may be used to
 *     endorse or promote products derived from this software without specific
 *     prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef OPALX_SUBFIELD_H
#define OPALX_SUBFIELD_H

#include "Utility/TypeUtils.h"

#include "Field/BConds.h"
#include "Field/BareField.h"

#include "Meshes/UniformCartesian.h"

template <typename T, unsigned Dim, class Mesh, class Centering, class... ViewArgs>
class SubField : public ippl::BareField<T, Dim, ViewArgs...> {
    template <typename... Props>
    using base_type = SubField<T, Dim, Mesh, Centering, Props...>;

public:
    using Mesh_t      = Mesh;
    using Centering_t = Cell;
    using Layout_t    = ippl::FieldLayout<Dim>;
    using BareField_t = ippl::BareField<T, Dim, ViewArgs...>;
    using view_type   = typename BareField_t::view_type;
    using BConds_t    = ippl::BConds<SubField<T, Dim, Mesh, Centering, ViewArgs...>, Dim>;

    using uniform_type = typename ippl::detail::CreateUniformType<
            base_type, typename view_type::uniform_type>::type;

    // A default constructor, which should be used only if the user calls the
    // 'initialize' function before doing anything else.  There are no special
    // checks in the rest of the SubField methods to check that the SubField has
    // been properly initialized.
    SubField() : BareField_t(), mesh_m(nullptr), bc_m() {}

    SubField(const SubField&) = default;

    /*!
     * Creates a new SubField with the same properties and contents
     * @return A deep copy of the field
     */
    SubField deepCopy() const {
        SubField<T, Dim, Mesh, Centering, ViewArgs...> copy(
                *mesh_m, this->getLayout(), this->getNghost());
        Kokkos::deep_copy(copy.getView(), this->getView());
        return copy;
    }

    virtual ~SubField() = default;

    // Constructors including a Mesh object as argument:
    SubField(Mesh_t& m, Layout_t& l, int nghost = 1) { initialize(m, l, nghost); }

    // Initialize the SubField, also specifying a mesh
    void initialize(Mesh_t& m, Layout_t& l, int nghost = 1) {
        BareField_t::initialize(l, nghost);
        mesh_m = &m;
        for (unsigned int face = 0; face < 2 * Dim; ++face) {
            bc_m[face] = std::make_shared<ippl::NoBcFace<SubField<T, Dim, Mesh, Centering, ViewArgs...>>>(
                    face);
        }
    }

    // ML
    void updateLayout(Layout_t& l, int nghost = 1) { BareField_t::updateLayout(l, nghost); }

    void setFieldBC(BConds_t& bc) {
        bc_m = bc;
        bc_m.findBCNeighbors(*this);
    }

    // Access to the mesh
    KOKKOS_INLINE_FUNCTION Mesh_t& get_mesh() const { return *mesh_m; }

    /*!
     * Use the midpoint rule to calculate the field's volume integral
     * @return Integral of the field over its domain
     */
    T getVolumeIntegral() const {
        typename Mesh::value_type dV = mesh_m->getCellVolume();
        return this->sum() * dV;
    }

    /*!
     * Use the midpoint rule to calculate the field's volume average
     * @return Integral of the field divided by the mesh volume
     */
    T getVolumeAverage() const { return getVolumeIntegral() / mesh_m->getMeshVolume(); }

    BConds_t& getFieldBC() { return bc_m; }
    // Assignment from constants and other arrays.
    using ippl::BareField<T, Dim, ViewArgs...>::operator=;

private:
    // The Mesh object, and a flag indicating if we constructed it
    Mesh_t* mesh_m;

    // The boundary conditions.
    BConds_t bc_m;
};

template <typename T, unsigned Dim, class Mesh, class Centering, class... ViewArgs>
struct ippl::detail::isExpression<SubField<T, Dim, Mesh, Centering, ViewArgs...>> : std::true_type {};


#endif  // OPALX_SUBFIELD_H
