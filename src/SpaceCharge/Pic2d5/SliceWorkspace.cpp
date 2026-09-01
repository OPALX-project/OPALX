/**
 * @file SliceWorkspace.cpp
 * @brief Implements persistent FFT2D5 field and slice-solver construction.
 */

#include "SpaceCharge/Pic2d5/SliceWorkspace.h"

#include "Utilities/OpalException.h"

#include <array>

namespace opalx::spacecharge {
    namespace {

        ippl::NDIndex<3> makeDomain3(const std::array<std::size_t, 3>& meshSize) {
            ippl::NDIndex<3> domain;
            for (unsigned dimension = 0; dimension < 3; ++dimension) {
                domain[dimension] = ippl::Index(static_cast<int>(meshSize[dimension]));
            }
            return domain;
        }

        ippl::NDIndex<2> makeDomain2(const std::array<std::size_t, 3>& meshSize) {
            ippl::NDIndex<2> domain;
            domain[0] = ippl::Index(static_cast<int>(meshSize[0]));
            domain[1] = ippl::Index(static_cast<int>(meshSize[1]));
            return domain;
        }

        SliceWorkspace::Vector3 makeSize(const Pic2d5Config& config, double pathLength) {
            SliceWorkspace::Vector3 size;
            size[0] = config.pipeSizeX();
            size[1] = config.pipeSizeY();
            size[2] = pathLength;
            return size;
        }

        SliceWorkspace::Vector3 makeOrigin(const Pic2d5Config& config) {
            SliceWorkspace::Vector3 origin;
            origin[0] = -0.5 * config.pipeSizeX();
            origin[1] = -0.5 * config.pipeSizeY();
            origin[2] = 0.0;
            return origin;
        }

        SliceWorkspace::Vector3 makeSpacing(const Pic2d5Config& config, double pathLength) {
            SliceWorkspace::Vector3 spacing;
            const auto& meshSize = config.meshSize();
            spacing[0]           = config.pipeSizeX() / static_cast<double>(meshSize[0]);
            spacing[1]           = config.pipeSizeY() / static_cast<double>(meshSize[1]);
            spacing[2]           = pathLength / static_cast<double>(meshSize[2]);
            return spacing;
        }

    }  // namespace

    SliceWorkspace::SliceWorkspace(const Pic2d5Config& config, double pathLength)
        : meshSize_m(config.meshSize()),
          spacing_m(makeSpacing(config, pathLength)),
          origin_m(makeOrigin(config)),
          size_m(makeSize(config, pathLength)),
          domain_m(makeDomain3(meshSize_m)),
          mesh_m(domain_m, spacing_m, origin_m),
          layout_m(MPI_COMM_WORLD, domain_m, std::array<bool, 3>{false, false, false}, true),
          sliceDomain_m(makeDomain2(meshSize_m)),
          sliceMesh_m(
                  sliceDomain_m, Vector2(spacing_m[0], spacing_m[1]),
                  Vector2(origin_m[0], origin_m[1])),
          sliceLayout_m(MPI_COMM_WORLD, sliceDomain_m, std::array<bool, 2>{false, false}) {
        if (!(pathLength > 0.0)) {
            throw OpalException(
                    "SliceWorkspace::SliceWorkspace",
                    "The FFT2D5 reference-path length must be positive.");
        }

        chargeDensity_m.initialize(mesh_m, layout_m);
        electricField_m.initialize(mesh_m, layout_m);

        solverParameters_m.add("use_heffte_defaults", false);
        solverParameters_m.add("use_pencils", true);
        solverParameters_m.add("use_gpu_aware", true);
        solverParameters_m.add("comm", ippl::a2av);
        solverParameters_m.add("r2c_direction", 0);
        solverParameters_m.add("algorithm", OpenSolver2::HOCKNEY);
        solverParameters_m.add("output_type", OpenSolver2::SOL_AND_GRAD);

        slices_m.resize(meshSize_m[2]);
        for (Slice& slice : slices_m) {
            slice.electricField = std::make_shared<VectorField2>(sliceMesh_m, sliceLayout_m);
            slice.chargeDensity = std::make_shared<ScalarField2>(sliceMesh_m, sliceLayout_m);
            slice.solver        = std::make_shared<OpenSolver2>(
                    *slice.electricField, *slice.chargeDensity, solverParameters_m);
        }
    }

}  // namespace opalx::spacecharge
