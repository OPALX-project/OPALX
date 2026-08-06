//
// Copyright (c) 2008 - 2026, Paul Scherrer Institut, Villigen PSI, Switzerland
//
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
#ifndef OPALX_SOLVE2D5_HPP_
#define OPALX_SOLVE2D5_HPP_
#include <Kokkos_NumericTraits.hpp>
#include <filesystem>
#include "AbstractObjects/OpalData.h"
#include "Fields/Interpolation/MMatrix.h"
#include "Utilities/Util.h"

template <typename T>
Solve2d5<T>::Solve2d5(
        PartBunch_t* partBunch, std::string solver, Field_t<3U>* rho, VField_t<T, 3U>* E,
        Field_t<3U>* phi, std::shared_ptr<BCHandler_t> bcHandler, const Vector<int, 3U>& nR,
        const LongitudinalFieldMode longitudinalFieldMode, const T pipeSizeX, const T pipeSizeY,
        const T beamRadius, const bool closedRing, bool scatterChargeLongitudinally,
        const std::string& refPathFileName)
    : Base(solver, rho, E, phi, bcHandler, 0, true),
      partBunch_m(partBunch),
      beamRadius_m(beamRadius),
      longitudinalFieldMode_m(longitudinalFieldMode),
      closedRing_m(closedRing),
      scatterChargeLongitudinally_m(scatterChargeLongitudinally),
      nR_m(nR),
      referencePathFileName_m(refPathFileName),
      solver_m(solver),
      pipeSizeX_m(pipeSizeX),
      pipeSizeY_m(pipeSizeY) {}

template <typename T>
void Solve2d5<T>::orbitThreadersReady() {
    if (ippl::Comm->size() != 1) {
        throw OpalException(
                "Solve2d5::orbitThreadersReady",
                "FFT2D5 currently supports only one MPI rank. Distributed fields and ORB load "
                "balancing are not implemented for this solver.");
    }

    // Load the reference path to determine the Frenet-Serret domain dimensions
    auto pathLength = loadReferencePath();
    sizer_m         = {pipeSizeX_m, pipeSizeY_m, pathLength};
    originr_m       = {-pipeSizeX_m / 2, -pipeSizeY_m / 2, 0};
    hr_m            = sizer_m / nR_m;
    for (unsigned i = 0; i < Dim; i++) {
        domain_m[i] = ippl::Index(nR_m[i]);
    }
    // Create the 3D fields and field container
    // Why are the first three parameters to the FieldContainer_t constructor references?
    Vector3D_t rmax = originr_m + sizer_m;
    partBunch_m->setFieldContainer(
            std::make_shared<FieldContainer_t>(
                    hr_m, originr_m, rmax, std::array{false, false, false}, domain_m, originr_m,
                    true));
    partBunch_m->getFieldContainer()->initializeFields(solver_m);
    rho_m = &partBunch_m->getFieldContainer()->getRho();
    E_m   = &partBunch_m->getFieldContainer()->getE();
    // Set the field solver base class members to reflect changes to the field container
    FieldSolver<T, 3U>::setRho(rho_m);
    FieldSolver<T, 3U>::setE(E_m);
    FieldSolver<T, 3U>::setPhi(&partBunch_m->getFieldContainer()->getPhi());
    // The grid dimensions
    ippl::NDIndex ndIndex2d(Vector<unsigned, 2U>{nR_m.data_m[0], nR_m.data_m[1]});
    auto numSlices = nR_m.data_m[2];
    // The slice mesh
    Vector2D_t spacing(hr_m.data_m[0], hr_m.data_m[1]);
    Vector2D_t origin(originr_m.data_m[0], originr_m.data_m[1]);
    sliceMesh_m = std::make_shared<Mesh2D_t>(ndIndex2d, spacing, origin);
    // The slice layout
    auto* layout3d = &rho_m->getLayout();
    std::array isParallel{layout3d->isParallel()[0], layout3d->isParallel()[1]};
    sliceLayout_m = std::make_shared<Layout2D_t>(layout3d->comm, ndIndex2d, isParallel);
    // Solver parameters
    solverParams_m = std::make_unique<ippl::ParameterList>();
    solverParams_m->add("use_heffte_defaults", false);
    solverParams_m->add("use_pencils", true);
    solverParams_m->add("use_gpu_aware", true);
    solverParams_m->add("comm", ippl::a2av);
    solverParams_m->add("r2c_direction", 0);
    solverParams_m->add("algorithm", OpenSolver2D_t::HOCKNEY);
    solverParams_m->add("output_type", OpenSolver2D_t::SOL_AND_GRAD);
    // Create the slices and their solvers
    // Note that here we create new Kokkos arrays for the slices so we are
    // going to have to copy data into and out of them during the solve calls.
    // If instead the Field type would accept a Kokkos subview rather than always
    // creating its own Kokkos array, we could avoid the copy operations.
    twoDSolvers_m.resize(numSlices);
    for (size_t z = 0; z < numSlices; ++z) {
        twoDSolvers_m[z].E_m      = std::make_shared<VField_t<T, 2U>>(*sliceMesh_m, *sliceLayout_m);
        twoDSolvers_m[z].rho_m    = std::make_shared<Field_t<2U>>(*sliceMesh_m, *sliceLayout_m);
        twoDSolvers_m[z].solver_m = std::make_shared<OpenSolver2D_t>(
                *twoDSolvers_m[z].E_m, *twoDSolvers_m[z].rho_m, *solverParams_m);
    }
    lineDensity_m         = LineDensityView_t("lineDensity", numSlices + LineDensityGhostCells);
    lineDensityGradient_m = LineDensityView_t("lineDensityGradient", numSlices);
}

template <typename T>
template <typename DiagnosticPolicy>
void Solve2d5<T>::doRunSolver(DiagnosticPolicy diagnostic) {
    scatterToGrid<DiagnosticPolicy>(*partBunch_m, diagnostic);
    solvePoissons<DiagnosticPolicy>(diagnostic);
    calculateLineDensity<DiagnosticPolicy>(diagnostic);
    gatherFromGrid<DiagnosticPolicy>(*partBunch_m, diagnostic);
}

template <typename T>
template <typename DiagnosticPolicy>
std::unique_ptr<DiagnosticPolicy> Solve2d5<T>::createDiagnostic(
        typename NullDiagnostic::Kind kind) {
    auto result = std::make_unique<DiagnosticPolicy>(kind);
    result->initialise(*partBunch_m, *rho_m, lineDensity_m, lineDensityGradient_m, *E_m);
    return result;
}

template <typename T>
T Solve2d5<T>::loadReferencePath() {
    // Return the total length of the path
    T length{};
    // The path name of the file created by the OrbitThreader
    auto opal            = OpalData::getInstance();
    std::string fileName = referencePathFileName_m;
    if (fileName.empty()) {
        fileName = Util::combineFilePath(
                {opal->getAuxiliaryOutputDirectory(),
                 opal->getInputBasename() + "_DesignPath.dat"});
    }
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
        referencePath.back().data_m[0] = rx;  // Because IPPL fails to implement
        referencePath.back().data_m[1] = ry;  // an appropriate constructor
        referencePath.back().data_m[2] = rz;  //
        // Accumulate the segment lengths
        if (referencePath.size() > 1) {
            auto index       = referencePath.size() - 1;
            Vector3D_t delta = referencePath[index] - referencePath[index - 1];
            length += delta.Pnorm();
        }
    }
    // Transfer into a device Kokkos view
    referencePath_m = ReferenceView_t("ref", referencePath.size());
    auto hostView   = Kokkos::create_mirror_view(referencePath_m);
    for (size_t i = 0; i < referencePath.size(); ++i) {
        hostView(i) = referencePath[i];
    }
    Kokkos::deep_copy(referencePath_m, hostView);
    // Return the total length of the path
    return length;
}

// Place the charge from all the particles into the charge density grid.
// Multiple particle containers in the bunch are supported.
// The particle container array dt is temporarily used to contain the charge to deposit.
// The real work is handled by the doScatterToGrid Kokkos kernel.
// The periodic boundary conditions are handled when the ring is marked as closed.
template <typename T>
template <typename DiagnosticPolicy>
void Solve2d5<T>::scatterToGrid(const PartBunch_t& bunch, DiagnosticPolicy diagnostic) {
    if (referencePath_m.extent(0) > 1) {
        const auto& ref    = referencePath_m;
        const auto invDr   = 1.0 / hr_m;
        const int nghost   = rho_m->getNghost();
        const auto& layout = rho_m->getLayout();
        const auto& lDom   = layout.getLocalNDIndex();
        const auto rhoView = rho_m->getView();
        const auto& origin = originr_m;
        // Do the scattering for all the particle containers
        Kokkos::deep_copy(rhoView, 0.0);
        for (auto& pcs = bunch.getParticleContainers(); const auto& pc : pcs) {
            pc->updateMoments();
            pc->scaleDtByCharge();
            const auto& r       = pc->R.getView();
            const auto& p       = pc->P.getView();
            const auto meanPs   = pc->getMeanP().data_m[2];
            const auto& dt      = pc->dt.getView();
            const auto& invalid = pc->InvalidMask.getView();
            if (scatterChargeLongitudinally_m) {
                Kokkos::parallel_for(
                        "Solve2d5::scatterToGrid::scatter", pc->getLocalNum(),
                        KOKKOS_LAMBDA(const size_t n) {
                            doScatterToGrid<true>(
                                    n, r, p, ref, meanPs, dt, invalid, invDr, nghost, lDom, rhoView,
                                    origin, diagnostic);
                        });
            } else {
                Kokkos::parallel_for(
                        "Solve2d5::scatterToGrid::scatter", pc->getLocalNum(),
                        KOKKOS_LAMBDA(const size_t n) {
                            doScatterToGrid<false>(
                                    n, r, p, ref, meanPs, dt, invalid, invDr, nghost, lDom, rhoView,
                                    origin, diagnostic);
                        });
            }
            Kokkos::fence();
            diagnostic.scatterCharge(rhoView);
            pc->unscaleDtByCharge();
        }
        // Handle the closed ring periodic boundary condition
        if (closedRing_m) {
            using Policy2D_t          = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;
            constexpr auto firstRealZ = LineDensityFirstRealCell;
            const auto lastRealZ      = nR_m.data_m[2] - 1 + LineDensityFirstRealCell;
            Kokkos::parallel_for(
                    "Solve2d5::scatterToGrid::boundaries",
                    Policy2D_t({0, 0}, {nR_m.data_m[0], nR_m.data_m[1]}),
                    KOKKOS_LAMBDA(const size_t i, const size_t j) {
                        rhoView(i, j, firstRealZ) += rhoView(i, j, lastRealZ + 1);
                        rhoView(i, j, lastRealZ) += rhoView(i, j, firstRealZ - 1);
                        rhoView(i, j, lastRealZ + 1)  = 0;
                        rhoView(i, j, firstRealZ - 1) = 0;
                    });
            Kokkos::fence();
        }
        // Now scale by volume and time step to get the proper charge density.
        // ippl::apply is a function that accesses a view using indices in an array like structure
        // and is from IpplOperations.h
        const auto cellVolume =
                std::reduce(hr_m.begin(), hr_m.end(), 1.0, std::multiplies<double>());
        const auto scale = bunch.getdT() * cellVolume;
        ippl::parallel_for(
                "Solve2d5::scatterToGrid::scale", rho_m->getFieldRangePolicy(),
                KOKKOS_LAMBDA(const ippl::RangePolicy<Dim>::index_array_type& idx) {
                    ippl::apply(rhoView, idx) = ippl::apply(rhoView, idx) / scale;
                });
        Kokkos::fence();
        diagnostic.scatterChargeDensity(rhoView);
    }
}

// This is a Kokkos kernel that deposits the charge for a single particle into the
// charge density grid.  The particle's coordinates are first translated into
// the Frenet-Serret coordinate system using the previously loaded reference path.
// The particle is then boosted into the beam's reference frame before the charge
// is scattered using one of two CiC algorithms depending on the ScatterLongitudinally
// template parameter.   When true, trilinear CiC is used.  When false, the particle
// is assigned to a single z slice and bilinear CiC is used in that slice.
template <typename T>
template <bool ScatterLongitudinally, typename DiagnosticPolicy>
KOKKOS_FUNCTION void Solve2d5<T>::doScatterToGrid(
        const size_t n, const VectorView_t& r, const VectorView_t& p, const ReferenceView_t& ref,
        const T meanPs, const ScalarView_t& dt, const BooleanView_t& invalid, Vector3D_t invDr,
        const int nghost, const ippl::NDIndex<3U> lDom, ScalarGridView3D_t rho, Vector3D_t origin,
        DiagnosticPolicy diagnostic) {
    if (!invalid(n)) {
        // Into Frenet-Serret coordinates
        Vector3D_t fsR, fsP, bUnit, nUnit, tUnit;
        convertToFrenetSerret(n, r, p, ref, fsR, fsP, bUnit, nUnit, tUnit);
        diagnostic.frenetSerretScatter(n, fsR, fsP, invalid(n));
        // Boost into the beam frame
        boostToBeamFrame(meanPs, fsP);
        diagnostic.boostToBeam(n, fsR, fsP, invalid(n));
        // CiC scatter the charge to the 3D rho grid
        scatterToRho<ScatterLongitudinally>(n, fsR, dt, invDr, nghost, lDom, rho, origin);
    }
}

template <typename T>
KOKKOS_FUNCTION void Solve2d5<T>::convertToFrenetSerret(
        const size_t n, const VectorView_t& r, const VectorView_t& p, const ReferenceView_t& ref,
        Vector3D_t& fsR, Vector3D_t& fsP, Vector3D_t& bUnit, Vector3D_t& nUnit, Vector3D_t& tUnit) {
    // Find the segment with the shortest normal to the point
    T bestDist2 = Kokkos::Experimental::finite_max_v<T>;
    size_t bestI{};
    T bestU{};
    Vector3D_t bestDi{};
    Vector3D_t bestRc{};
    auto rn       = r(n);
    auto segments = ref.extent(0);
    for (size_t i = 0; i < segments - 1; ++i) {
        auto refi     = ref(i);
        auto refip1   = ref(i + 1);
        Vector3D_t di = refip1 - refi;
        const T di2   = di.dot(di);
        if (di2 > 0.0) {
            const T u       = Kokkos::clamp(di.dot(rn - refi) / di2, 0.0, 1.0);
            Vector3D_t rc   = refi + u * di;
            Vector3D_t diff = rn - rc;
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
template <typename ViewType>
KOKKOS_FUNCTION bool Solve2d5<T>::makeWeights(
        Vector3D_t fsR, Vector3D_t origin, Vector3D_t invDr, int nghost,
        const ippl::NDIndex<3U>& lDom, const ViewType& view, ippl::Vector<T, 3U>& whi,
        ippl::Vector<T, 3U>& wlo, ippl::Vector<int, 3U>& args) {
    const auto l                = (fsR - origin) * invDr + 0.5;
    ippl::Vector<int, 3U> index = l;
    whi                         = l - index;
    wlo                         = 1.0 - whi;
    args                        = index - lDom.first() + nghost;
    // CIC touches args[d] and args[d] - 1, so valid args are
    // [1, extent - 1]. Anything outside would underflow or
    // overrun the field view on the device.
    bool inBounds = args[0] > 0 && args[0] < static_cast<int>(view.extent(0));
    inBounds      = inBounds && args[1] > 0 && args[1] < static_cast<int>(view.extent(1));
    inBounds      = inBounds && args[2] > 0 && args[2] < static_cast<int>(view.extent(2));
    return inBounds;
}

template <typename T>
template <bool ScatterLongitudinally>
KOKKOS_FUNCTION void Solve2d5<T>::scatterToRho(
        const size_t n, Vector3D_t fsR, const ScalarView_t& dt, Vector3D_t invDr, const int nghost,
        const ippl::NDIndex<3U>& lDom, ScalarGridView3D_t rho, Vector3D_t origin) {
    // CiC scatter the charge to the 3D rho grid
    ippl::Vector<T, Dim> whi, wlo;
    ippl::Vector<int, Dim> args;
    if (makeWeights(fsR, origin, invDr, nghost, lDom, rho, whi, wlo, args)) {
        if constexpr (ScatterLongitudinally) {
            scatter3D(rho, wlo, whi, args[0], args[1], args[2], dt(n));
        } else {
            scatter2D(rho, wlo, whi, args[0], args[1], args[2], dt(n));
        }
    }
}

template <typename T>
template <typename DiagnosticPolicy>
void Solve2d5<T>::calculateLineDensity(DiagnosticPolicy diagnostic) {
    if (longitudinalFieldMode_m != LongitudinalFieldMode::None) {
        using Policy2D_t       = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;
        const auto rho3d       = rho_m->getView();
        auto deviceLineDensity = lineDensity_m;
        auto hostLineDensity   = Kokkos::create_mirror_view(lineDensity_m);
        const auto numSlices   = nR_m.data_m[2];
        // Calculate the total charge density for each z slice
        for (size_t k = 0; k < numSlices + LineDensityGhostCells; ++k) {
            T sum{};
            Kokkos::parallel_reduce(
                    Policy2D_t({0, 0}, {rho3d.extent(0), rho3d.extent(1)}),
                    KOKKOS_LAMBDA(const size_t i, const size_t j, T& localSum) {
                        localSum += rho3d(i, j, k);
                    },
                    sum);
            hostLineDensity(k) = sum;
        }
        // Set the ghost cells to the boundary conditions
        if (closedRing_m) {
            hostLineDensity(0) = hostLineDensity(numSlices + LineDensityGhostCells - 2);
            hostLineDensity(numSlices + LineDensityGhostCells - 1) =
                    hostLineDensity(LineDensityFirstRealCell);
        } else {
            hostLineDensity(0)                                     = 0.0;
            hostLineDensity(numSlices + LineDensityGhostCells - 1) = 0.0;
        }
        Kokkos::deep_copy(deviceLineDensity, hostLineDensity);
        diagnostic.totalDensity(lineDensity_m);
        // Convert this to line density
        auto dx = hr_m[0];
        auto dy = hr_m[1];
        Kokkos::parallel_for(
                "Solve2d5::calculateLineDensity::convert", numSlices + LineDensityGhostCells,
                KOKKOS_LAMBDA(const size_t k) { deviceLineDensity(k) *= dx * dy; });
        Kokkos::fence();
        diagnostic.lineDensity(lineDensity_m);
        // Find the gradient
        auto lineDensityGradient = lineDensityGradient_m;
        const auto dz            = rho_m->get_mesh().getMeshSpacing().data_m[2];
        Kokkos::parallel_for(
                "Solve2d5::calculateLineDensity::gradient", numSlices,
                KOKKOS_LAMBDA(const size_t k) {
                    lineDensityGradient(k) = (deviceLineDensity(k + LineDensityFirstRealCell + 1)
                                              - deviceLineDensity(k + LineDensityFirstRealCell - 1))
                                             / (2.0 * dz);
                });
        Kokkos::fence();
        diagnostic.lineDensityGradient(lineDensityGradient_m);
    } else {
        Kokkos::deep_copy(lineDensity_m, 0);
        diagnostic.totalDensity(lineDensity_m);
        diagnostic.lineDensity(lineDensity_m);
        Kokkos::deep_copy(lineDensityGradient_m, 0);
        diagnostic.lineDensityGradient(lineDensityGradient_m);
    }
}

template <typename T>
template <typename DiagnosticPolicy>
void Solve2d5<T>::solvePoissons(DiagnosticPolicy diagnostic) {
    using Policy = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;
    auto e3d     = E_m->getView();
    auto nghost  = E_m->getNghost();
    // Copy the 3D charge density grid into the array of 2D grids, solve the 2D poisson on each,
    // then copy the E field results into the 3D E field grid.
    for (size_t z = 0; z < twoDSolvers_m.size(); ++z) {
        auto& s = twoDSolvers_m[z];
        // Do the 2D solve to get Ex and Ey
        auto rho2d = s.rho_m->getView();
        Kokkos::deep_copy(
                rho2d, Kokkos::subview(rho_m->getView(), Kokkos::ALL(), Kokkos::ALL(), z + nghost));
        // Scale by the coupling constant then solve
        Kokkos::parallel_for(
                "Solve2d5::solvePoissons::coupling",
                Policy({0, 0}, {rho2d.extent(0), rho2d.extent(1)}),
                KOKKOS_LAMBDA(const size_t i, const size_t j) {
                    rho2d(i, j) *= 1 / Physics::epsilon_0;
                });
        s.solver_m->solve();
        Kokkos::fence();
        auto e2d = s.E_m->getView();
        // Copy the 2D E field into the 3D E field grid
        Kokkos::parallel_for(
                "Solve2d5::solvePoissons::copy", Policy({0, 0}, {e2d.extent(0), e2d.extent(1)}),
                KOKKOS_LAMBDA(const size_t i, const size_t j) {
                    e3d(i, j, z + nghost)[0] = e2d(i, j)[0];
                    e3d(i, j, z + nghost)[1] = e2d(i, j)[1];
                    e3d(i, j, z + nghost)[2] = 0.0;
                });
        Kokkos::fence();
    }
    // Set the ghost slices
    constexpr auto leftGhost = 0;
    auto rightGhost          = e3d.extent(2) - 1;
    if (closedRing_m) {
        Kokkos::parallel_for(
                "Solve2d5::solvePoissons::ghostClosed",
                Policy({0, 0}, {e3d.extent(0), e3d.extent(1)}),
                KOKKOS_LAMBDA(const size_t i, const size_t j) {
                    e3d(i, j, leftGhost)[0]  = e3d(i, j, rightGhost - 1)[0];
                    e3d(i, j, leftGhost)[1]  = e3d(i, j, rightGhost - 1)[1];
                    e3d(i, j, leftGhost)[2]  = e3d(i, j, rightGhost - 1)[2];
                    e3d(i, j, rightGhost)[0] = e3d(i, j, leftGhost + 1)[0];
                    e3d(i, j, rightGhost)[1] = e3d(i, j, leftGhost + 1)[1];
                    e3d(i, j, rightGhost)[2] = e3d(i, j, leftGhost + 1)[2];
                });
        Kokkos::fence();
    } else {
        Kokkos::parallel_for(
                "Solve2d5::solvePoissons::ghostOpen",
                Policy({0, 0}, {e3d.extent(0), e3d.extent(1)}),
                KOKKOS_LAMBDA(const size_t i, const size_t j) {
                    e3d(i, j, leftGhost)[0]  = 0;
                    e3d(i, j, leftGhost)[1]  = 0;
                    e3d(i, j, leftGhost)[2]  = 0;
                    e3d(i, j, rightGhost)[0] = 0;
                    e3d(i, j, rightGhost)[1] = 0;
                    e3d(i, j, rightGhost)[2] = 0;
                });
        Kokkos::fence();
    }
    diagnostic.eField(e3d);
}

template <typename T>
template <typename DiagnosticPolicy>
void Solve2d5<T>::gatherFromGrid(const PartBunch_t& bunch, DiagnosticPolicy diagnostic) {
    if (referencePath_m.extent(0) > 1) {
        for (auto& pcs = bunch.getParticleContainers(); const auto& pc : pcs) {
            const auto& r                  = pc->R.getView();
            const auto& p                  = pc->P.getView();
            const auto& ref                = referencePath_m;
            const auto meanPs              = pc->getMeanP().data_m[2];
            const auto& e                  = pc->E.getView();
            const auto& b                  = pc->B.getView();
            const auto& invalid            = pc->InvalidMask.getView();
            const auto& dr                 = hr_m;
            const auto invDr               = 1.0 / dr;
            const int nghost               = rho_m->getNghost();
            const auto& layout             = rho_m->getLayout();
            const auto& lDom               = layout.getLocalNDIndex();
            const auto& origin             = originr_m;
            const auto gammaB              = Kokkos::sqrt(1.0 + meanPs * meanPs);
            const auto betaB               = meanPs / gammaB;
            const auto lineDensityGradient = lineDensityGradient_m;
            const auto eField              = E_m->getView();
            const auto pipeRadius          = std::min(sizer_m.data_m[0], sizer_m.data_m[1]);
            T gBy4PiEpsilon0;
            if (longitudinalFieldMode_m == LongitudinalFieldMode::Cylindrical) {
                gBy4PiEpsilon0 = CircularPipeG0 + 2 * Kokkos::log(pipeRadius / beamRadius_m);
            } else if (longitudinalFieldMode_m == LongitudinalFieldMode::Plates) {
                gBy4PiEpsilon0 = ParallelPlatesG0
                                 + 2 * Kokkos::log(4 * pipeRadius / Physics::pi / beamRadius_m);
            } else {
                gBy4PiEpsilon0 = OpenG0;
            }
            gBy4PiEpsilon0 /= 4 * Physics::pi * Physics::epsilon_0;
            Kokkos::parallel_for(
                    "Solve2d5::gatherFromGrid", pc->getLocalNum(), KOKKOS_LAMBDA(const size_t n) {
                        doGatherFromGrid(
                                n, r, p, ref, gammaB, betaB, e, b, invalid, invDr, nghost, lDom,
                                eField, origin, gBy4PiEpsilon0, lineDensityGradient, diagnostic);
                    });
            Kokkos::fence();
        }
    }
}

template <typename T>
template <typename DiagnosticPolicy>
KOKKOS_FUNCTION void Solve2d5<T>::doGatherFromGrid(
        const size_t n, const VectorView_t& r, const VectorView_t& p, const ReferenceView_t& ref,
        const T beamGamma, const T beamBeta, const VectorView_t& e, const VectorView_t& b,
        const BooleanView_t& invalid, Vector3D_t invDr, const int nghost,
        const ippl::NDIndex<3U> lDom, VectorGridView3D_t eField, Vector3D_t origin,
        T gBy4PiEpsilon0, LineDensityView_t lineDensityGradient, DiagnosticPolicy diagnostic) {
    if (!invalid(n)) {
        // Into Frenet-Serret coordinates
        Vector3D_t fsR, fsP, bUnit, nUnit, tUnit;
        convertToFrenetSerret(n, r, p, ref, fsR, fsP, bUnit, nUnit, tUnit);
        diagnostic.frenetSerretGather(n, fsR, fsP, invalid(n));
        // CiC Gather the boosted E field
        gatherFromEField(n, fsR, e, invDr, nghost, lDom, eField, origin);
        diagnostic.gatherEField(n, e(n), b(n), invalid(n));
        // Unboost from the beam frame
        unboostFromBeamFrame(n, beamGamma, beamBeta, e, b);
        diagnostic.deboostFromBeam(n, e(n), b(n), invalid(n));
        // Calculate the longitudinal E field in Frenet-Serret coordinates
        const int index = (fsR.data_m[2] - origin.data_m[2]) * invDr.data_m[2] + 0.5;
        const auto ldg  = lineDensityGradient(index);
        e(n).data_m[2]  = gBy4PiEpsilon0 * ldg / beamGamma / beamGamma;
        diagnostic.longitudinalField(n, e(n), b(n), invalid(n));
        // And finally back into lab coordinates
        convertFromFrenetSerret(n, bUnit, nUnit, tUnit, e, b);
        diagnostic.labFrameFields(n, e(n), b(n), invalid(n));
    }
}

template <typename T>
KOKKOS_FUNCTION void Solve2d5<T>::gatherFromEField(
        const size_t n, Vector3D_t fsR, const VectorView_t& e, Vector3D_t invDr, const int nghost,
        const ippl::NDIndex<3U>& lDom, VectorGridView3D_t eField, Vector3D_t origin) {
    // CiC gather the boosted E field to the 3D E field grid
    ippl::Vector<T, Dim> whi, wlo;
    ippl::Vector<int, Dim> args;
    if (makeWeights(fsR, origin, invDr, nghost, lDom, eField, whi, wlo, args)) {
        e(n) = gather2D(eField, wlo, whi, args[0], args[1], args[2] - 1)
               + gather2D(eField, wlo, whi, args[0], args[1], args[2]);
    }
}

template <typename T>
KOKKOS_FUNCTION Solve2d5<T>::Vector3D_t Solve2d5<T>::gather2D(
        VectorGridView3D_t eField, const ippl::Vector<T, 3U>& wlo, const ippl::Vector<T, 3U>& whi,
        int x, int y, int z) {
    Vector3D_t result;
    result = wlo[0] * wlo[1] * eField(x - 1, y - 1, z);
    result += whi[0] * wlo[1] * eField(x, y - 1, z);
    result += wlo[0] * whi[1] * eField(x - 1, y, z);
    result += whi[0] * whi[1] * eField(x, y, z);
    return result;
}

template <typename T>
KOKKOS_FUNCTION void Solve2d5<T>::scatter3D(
        ScalarGridView3D_t rho, const ippl::Vector<T, 3U>& wlo, const ippl::Vector<T, 3U>& whi,
        int x, int y, int z, T charge) {
    rho(x - 1, y - 1, z - 1) = wlo[0] * wlo[1] * wlo[2] * charge;
    rho(x, y - 1, z - 1)     = whi[0] * wlo[1] * wlo[2] * charge;
    rho(x - 1, y, z - 1)     = wlo[0] * whi[1] * wlo[2] * charge;
    rho(x, y, z - 1)         = whi[0] * whi[1] * wlo[2] * charge;
    rho(x - 1, y - 1, z)     = wlo[0] * wlo[1] * whi[2] * charge;
    rho(x, y - 1, z)         = whi[0] * wlo[1] * whi[2] * charge;
    rho(x - 1, y, z)         = wlo[0] * whi[1] * whi[2] * charge;
    rho(x, y, z)             = whi[0] * whi[1] * whi[2] * charge;
}

template <typename T>
KOKKOS_FUNCTION void Solve2d5<T>::scatter2D(
        ScalarGridView3D_t rho, const ippl::Vector<T, 3U>& wlo, const ippl::Vector<T, 3U>& whi,
        int x, int y, int z, T charge) {
    rho(x - 1, y - 1, z) = wlo[0] * wlo[1] * charge;
    rho(x, y - 1, z)     = whi[0] * wlo[1] * charge;
    rho(x - 1, y, z)     = wlo[0] * whi[1] * charge;
    rho(x, y, z)         = whi[0] * whi[1] * charge;
}

template <typename T>
KOKKOS_FUNCTION void Solve2d5<T>::unboostFromBeamFrame(
        size_t n, T beamGamma, T beamBeta, const VectorView_t& e, const VectorView_t& b) {
    // Transform the E field from the boosted beam frame into E and B fields
    e(n).data_m[0] *= beamGamma;
    e(n).data_m[1] *= beamGamma;
    b(n).data_m[0] = beamBeta * e(n).data_m[1] / Physics::c;
    b(n).data_m[1] = -beamBeta * e(n).data_m[0] / Physics::c;
    b(n).data_m[2] = 0.0;
}

template <typename T>
KOKKOS_FUNCTION void Solve2d5<T>::convertFromFrenetSerret(
        size_t n, const Vector3D_t& bUnit, const Vector3D_t& nUnit, const Vector3D_t& tUnit,
        const VectorView_t& e, const VectorView_t& b) {
    e(n) = e(n).data_m[0] * bUnit + e(n).data_m[1] * nUnit + e(n).data_m[2] * tUnit;
    b(n) = b(n).data_m[0] * bUnit + b(n).data_m[1] * nUnit + b(n).data_m[2] * tUnit;
}

#endif
