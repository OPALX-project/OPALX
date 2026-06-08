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
#ifndef OPALX_SOLVE2D5_HPP_
#define OPALX_SOLVE2D5_HPP_
#include "AbstractObjects/OpalData.h"
#include "Fields/Interpolation/MMatrix.h"
#include "Utilities/Util.h"

template <typename T>
Solve2d5<T>::Solve2d5(
        std::string solver, Field_t<3>* rho, VField_t<T, 3>* E, Field_t<3>* phi,
        std::shared_ptr<BCHandler_t> bcHandler)
    : Base(solver, rho, E, phi, bcHandler, 0, true), phi_m(phi), rho_m(rho), E_m(E) {}

template <typename T>
void Solve2d5<T>::initSolver() {
    // The grid dimensions
    ippl::NDIndex<2> ndIndex2d({rho_m->getDomain()[0], rho_m->getDomain()[1]});
    auto numSlices = rho_m->getDomain()[2].length();
    // The slice mesh
    const auto* mesh3d = &rho_m->get_mesh();
    Vector2D_t spacing({mesh3d->getMeshSpacing()[0], mesh3d->getMeshSpacing()[1]});
    Vector2D_t origin({mesh3d->getOrigin()[0], mesh3d->getOrigin()[1]});
    sliceMesh_m = std::make_shared<Mesh2D_t>(ndIndex2d, spacing, origin);
    // The slice layout
    auto* layout3d = &rho_m->getLayout();
    std::array isParallel{layout3d->isParallel()[0], layout3d->isParallel()[1]};
    sliceLayout_m = std::make_shared<Layout2D_t>(layout3d->comm, ndIndex2d, isParallel);
    // Solver parameters
    ippl::ParameterList params;
    params.add("use_heffte_defaults", false);
    params.add("use_pencils", true);
    params.add("use_gpu_aware", true);
    params.add("comm", ippl::a2av);
    params.add("r2c_direction", 0);
    params.add("algorithm", OpenSolver2D_t::HOCKNEY);
    params.add("output_type", OpenSolver2D_t::SOL_AND_GRAD);
    // Create the slices and their solvers
    // Note that here we create new Kokkos arrays for the slices so we are
    // going to have to copy data into and out of them during the solve calls.
    // If instead the Field type would accept a Kokkos subview rather than always
    // creating its own Kokkos array, we could avoid the copy operations.
    twoDSolvers_m.resize(numSlices);
    for (size_t z = 0; z < numSlices; ++z) {
        twoDSolvers_m[z].E_m      = std::make_shared<VField_t<T, 2>>(*sliceMesh_m, *sliceLayout_m);
        twoDSolvers_m[z].rho_m    = std::make_shared<Field_t<2>>(*sliceMesh_m, *sliceLayout_m);
        twoDSolvers_m[z].solver_m = std::make_shared<OpenSolver2D_t>(
                *twoDSolvers_m[z].E_m, *twoDSolvers_m[z].rho_m, params);
    }
}

template <typename T>
void Solve2d5<T>::runSolver() {
    // std::get<OpenSolver_t<double, 3>>(this->getSolver()).solve();  // For now.

    // This is the copying required
#if 0
        Kokkos::deep_copy(
                twoDSolvers_m[z].rho_m,
                Kokkos::subview(phi_m->getView(), Kokkos::ALL(), Kokkos::ALL(), z));
        Kokkos::deep_copy(
                Kokkos::subview(rho_m->getView(), Kokkos::ALL(), Kokkos::ALL(), z),
                twoDSolvers_m[z].rho_m);
        Kokkos::deep_copy(
                Kokkos::subview(E_m->getView(), Kokkos::ALL(), Kokkos::ALL(), z),
                twoDSolvers_m[z].E_m);
#endif
}

template <typename T>
void Solve2d5<T>::loadReferencePath() {
    // The path name of the file created by the OrbitThreader
    auto opal            = OpalData::getInstance();
    std::string fileName = Util::combineFilePath(
            {opal->getAuxiliaryOutputDirectory(), opal->getInputBasename() + "_DesignPath.dat"});
    // Open the file
    if (!std::filesystem::exists(fileName)) {
        throw OpalException("Solve2d5::loadReferencePath", "File does not exist: " + fileName);
    }
    std::ifstream file(fileName);
    std::string line;
    std::getline(file, line);  // Skip the header line
    // Now read the remaining lines
    std::vector<Vector3D_t> referencePath;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        double s, rx, ry, rz, px, py, pz, ex, ey, ez, bx, by, bz, ke, t;
        iss >> s >> rx >> ry >> rz >> px >> py >> pz >> ex >> ey >> ez >> bx >> by >> bz >> ke >> t;
        if (iss.fail()) {
            throw OpalException("Solve2d5::loadReferencePath", "Failed to parse line: " + line);
        }
        referencePath.emplace_back();
        referencePath.back().data_m[0] = rx;  // Ugh 'cos IPPL fails to implement
        referencePath.back().data_m[1] = ry;  // the constructor that takes the
        referencePath.back().data_m[2] = rz;  // std::array<>
    }
    // Transfer into a device Kokkos view
    referencePath_m = ReferenceView_t("ref", referencePath.size());
    auto hostView   = Kokkos::create_mirror_view(referencePath_m);
    for (size_t i = 0; i < referencePath.size(); ++i) {
        hostView(i) = referencePath[i];
    }
    Kokkos::deep_copy(referencePath_m, hostView);
}

template <typename T>
template <typename DiagnosticPolicy>
void Solve2d5<T>::scatterToGrid(const PartBunch_t& bunch, DiagnosticPolicy diagnostic) {
    if (referencePath_m.extent(0) > 1) {
        for (auto& pcs = bunch.getParticleContainers(); const auto& pc : pcs) {
            pc->updateMoments();
            pc->scaleDtByCharge();  // Work out a way to include this in the parallel for
            const auto& r       = pc->R.getView();
            const auto& p       = pc->P.getView();
            const auto& ref     = referencePath_m;
            const auto meanPs   = pc->getMeanP().data_m[2];
            const auto& dt      = pc->dt.getView();
            const auto& invalid = pc->InvalidMask.getView();
            const auto& mesh    = rho_m->get_mesh();
            const auto& dr      = mesh.getMeshSpacing();
            const auto invDr    = 1.0 / dr;
            const int nghost    = rho_m->getNghost();
            const auto& layout  = rho_m->getLayout();
            const auto& lDom    = layout.getLocalNDIndex();
            const auto rhoView  = rho_m->getView();
            const auto& origin = mesh.getOrigin();
            Kokkos::deep_copy(rhoView, 0.0);
            Kokkos::parallel_for(
                    "Solve2d5::scatterToGrid()", pc->getLocalNum(), KOKKOS_LAMBDA(const size_t n) {
                        doScatterToGrid(
                                n, r, p, ref, meanPs, dt, invalid, invDr, nghost, lDom, rhoView,
                                origin, diagnostic);
                    });
            Kokkos::fence();
            diagnostic.scatterCharge(rhoView);
            pc->unscaleDtByCharge();  // Work out a way to include this in the parallel for
        }
    }
}

template <typename T>
template <typename DiagnosticPolicy>
KOKKOS_FUNCTION void Solve2d5<T>::doScatterToGrid(
        const size_t n, const VectorView_t& r, const VectorView_t& p, const ReferenceView_t& ref,
        const T meanPs, const ScalarView_t& dt, const BooleanView_t& invalid, Vector3D_t invDr,
        const int nghost, const ippl::NDIndex<3> lDom, ScalarGridView3D_t rho, Vector3D_t origin,
        DiagnosticPolicy diagnostic) {
    if (!invalid(n)) {
        // Into Frenet-Serret coordinates
        Vector3D_t fsR, fsP, bUnit, nUnit, tUnit;
        convertToFrenetSerret(n, r, p, ref, fsR, fsP, bUnit, nUnit, tUnit);
        diagnostic.frenetSerret(n, fsR, fsP, invalid(n));
        // Boost into the beam frame
        boostToBeamFrame(meanPs, fsP);
        diagnostic.boostToBeam(n, fsR, fsP, invalid(n));
        // CiC scatter the charge to the 3D rho grid
        scatterToRho(n, fsR, dt, invDr, nghost, lDom, rho, origin);
    }
}

template <typename T>
KOKKOS_FUNCTION void Solve2d5<T>::convertToFrenetSerret(
        const size_t n, const VectorView_t& r, const VectorView_t& p, const ReferenceView_t& ref,
        Vector3D_t& fsR, Vector3D_t& fsP, Vector3D_t& bUnit, Vector3D_t& nUnit, Vector3D_t& tUnit) {
    // Find the segment with the shortest normal to the point
    T bestDist2 = std::numeric_limits<T>::max();
    size_t bestI{};
    T bestU{};
    Vector3D_t bestDi{};
    Vector3D_t bestRc{};
    for (size_t i = 0; i < ref.size() - 1; ++i) {
        Vector3D_t di = ref(i + 1) - ref(i);
        const T di2   = di.dot(di);
        if (di2 > 0.0) {
            const T u       = Kokkos::clamp(di.dot(r(n) - ref(i)) / di2, 0.0, 1.0);
            Vector3D_t rc   = ref(i) + u * di;
            Vector3D_t diff = r(n) - rc;
            const T dist2   = diff.dot(diff);
            if (dist2 < bestDist2) {
                bestI     = i;
                bestU     = u;
                bestDi    = di;
                bestDist2 = dist2;
                bestRc    = rc;
            }
        }
    }
    // The basis unit vectors
    tUnit = bestDi / euclidean_norm(bestDi);
    nUnit = {0, 1, 0};
    bUnit = ippl::cross(nUnit, tUnit);
    // Calculate the S coordinate
    fsR.data_m[2] = bestU * euclidean_norm(bestDi);
    for (size_t i = 1; i <= bestI; ++i) {
        Vector3D_t diff = ref(i) - ref(i - 1);
        fsR.data_m[2] += euclidean_norm(diff);
    }
    // Remaining position coordinates
    const Vector3D_t rDiff = r(n) - bestRc;
    fsR.data_m[0]          = rDiff.dot(bUnit);
    fsR.data_m[1]          = rDiff.dot(nUnit);
    // The momentum coordinates
    fsP.data_m[0] = p(n).dot(bUnit);
    fsP.data_m[1] = p(n).dot(nUnit);
    fsP.data_m[2] = p(n).dot(tUnit);
}

template <typename T>
KOKKOS_FUNCTION void Solve2d5<T>::boostToBeamFrame(const T meanPs, Vector3D_t& fsP) {
    // Transform the longitudinal momentum coordinate into the beam reference frame
    auto gammaB   = Kokkos::sqrt(1.0 + meanPs * meanPs);
    auto gamma    = Kokkos::sqrt(1.0 + fsP.data_m[2] * fsP.data_m[2]);
    fsP.data_m[2] = gammaB * fsP.data_m[2] - meanPs * gamma;
}

template <typename T>
KOKKOS_FUNCTION void Solve2d5<T>::scatterToRho(
        const size_t n, Vector3D_t fsR, const ScalarView_t& dt, Vector3D_t invDr, const int nghost,
        const ippl::NDIndex<3>& lDom, ScalarGridView3D_t rho, Vector3D_t origin) {
    // CiC scatter the charge to the 3D rho grid
    const auto l                 = (fsR - origin) * invDr + 0.5;
    ippl::Vector<int, Dim> index = l;
    ippl::Vector<T, Dim> whi     = l - index;
    ippl::Vector<T, Dim> wlo     = 1.0 - whi;
    ippl::Vector<int, Dim> args  = index - lDom.first() + nghost;
    bool inBounds                = true;
    for (unsigned d = 0; d < Dim; ++d) {
        // CIC touches args[d] and args[d] - 1, so valid args are
        // [1, extent - 1]. Anything outside would underflow or
        // overrun the field view on the device.
        inBounds = inBounds && args[d] > 0 && args[d] < static_cast<int>(rho.extent(d));
    }
    if (inBounds) {
        ippl::Vector<size_t, Dim> viewArgs = args;
        ippl::detail::scatterToField(
                std::make_index_sequence<1 << Dim>{}, rho, wlo, whi, viewArgs, dt(n));
    }
}

template <typename T>
template <typename DiagnosticPolicy>
void Solve2d5<T>::solvePoissons(DiagnosticPolicy /*diagnostic*/) {
    // Copy the 3D charge density grid into the array of 2D grids, solve the 2D poisson on each,
    // then copy the E field results into the 3D E field grid.
    for (size_t z = 0; z < twoDSolvers_m.size(); ++z) {
        auto& s = twoDSolvers_m[z];
        // Do the 2D solve to get phi, Ex and Ey
        Kokkos::deep_copy(
            s.rho_m->getView(), Kokkos::subview(rho_m->getView(), Kokkos::ALL(), Kokkos::ALL(), z));
        s.solver_m->solve();   // 2D E is now filled in, the 2D rho array now contains 2D phi
        Kokkos::fence();
        // Copy into the 3D E field grid
        using Policy = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;
        auto e3d = E_m->getView();
        auto e2d = s.E_m->getView();
        Kokkos::parallel_for(
            "Zero Z",
            Policy({0,0}, {e2d.extent(0), e2d.extent(1)}),
            KOKKOS_LAMBDA(const int i, const int j) {
                e3d(i,j,z)[0] = e2d(i,j)[0];
                e3d(i,j,z)[1] = e2d(i,j)[1];
                e3d(i,j,z)[2] = 0.0;
            });
        Kokkos::fence();
    }
}


#endif
