/**
 * @file FieldContainer.hpp
 * @brief Temporary compatibility facade over the solver-owned PIC workspace.
 */

#ifndef OPAL_FIELD_CONTAINER_H
#define OPAL_FIELD_CONTAINER_H

#include <array>
#include <memory>
#include <string>

#include "SpaceCharge/Pic/PicWorkspace.h"
#include "Utilities/OpalException.h"

template <unsigned Dim>
using Mesh_t = ippl::UniformCartesian<double, Dim>;

template <typename T, unsigned Dim>
using PLayout_t = typename ippl::ParticleSpatialLayout<T, Dim, Mesh_t<Dim>>;

template <unsigned Dim>
using Centering_t = typename Mesh_t<Dim>::DefaultCentering;

template <unsigned Dim>
using FieldLayout_t = ippl::FieldLayout<Dim>;

template <typename T, unsigned Dim, class... ViewArgs>
using Field = ippl::Field<T, Dim, Mesh_t<Dim>, Centering_t<Dim>, ViewArgs...>;

template <unsigned Dim, class... ViewArgs>
using Field_t = Field<double, Dim, ViewArgs...>;

template <typename T, unsigned Dim>
using Vector_t = ippl::Vector<T, Dim>;

template <typename T, unsigned Dim, class... ViewArgs>
using VField_t = Field<Vector_t<T, Dim>, Dim, ViewArgs...>;

/**
 * @brief Forwards legacy callers to one stable @c PicWorkspace object.
 *
 * The facade and transitional solver bridge share ownership so existing samplers, particle
 * layouts, and the temporary IPPL PicManager surface continue to see the same mesh, layout, and
 * field addresses. The compatibility share is removed with the legacy manager path.
 */
template <typename T, unsigned Dim = 3>
class FieldContainer final {
public:
    using Workspace = opalx::spacecharge::PicWorkspace<T, Dim>;

    FieldContainer(
            Vector_t<T, Dim>& spacing, Vector_t<T, Dim>& lower, Vector_t<T, Dim>& upper,
            std::array<bool, Dim> decomposition, ippl::NDIndex<Dim> domain, Vector_t<T, Dim> origin,
            bool allPeriodic)
        : workspace_m(
                  std::make_shared<Workspace>(
                          spacing, lower, upper, decomposition, domain, origin, allPeriodic)) {}

    FieldContainer(const FieldContainer&)            = delete;
    FieldContainer& operator=(const FieldContainer&) = delete;

    [[nodiscard]] Workspace& workspace() { return requireWorkspace(); }
    [[nodiscard]] const Workspace& workspace() const { return requireWorkspace(); }

    [[nodiscard]] std::shared_ptr<Workspace> sharedWorkspace() {
        static_cast<void>(requireWorkspace());
        return workspace_m;
    }

    [[nodiscard]] VField_t<T, Dim>& getE() { return requireWorkspace().electricField(); }
    [[nodiscard]] Field_t<Dim>& getRho() { return requireWorkspace().chargeDensity(); }
    [[nodiscard]] Field<T, Dim>& getPhi() { return requireWorkspace().potential(); }

    [[nodiscard]] Mesh_t<Dim>& getMesh() { return requireWorkspace().mesh(); }
    [[nodiscard]] FieldLayout_t<Dim>& getFL() { return requireWorkspace().layout(); }

    [[nodiscard]] std::shared_ptr<VField_t<T, Dim>> getTempEField() {
        return {workspace_m, &requireWorkspace().accumulatedElectricField()};
    }
    [[nodiscard]] std::shared_ptr<VField_t<T, Dim>> getTempBField() {
        return {workspace_m, &requireWorkspace().accumulatedMagneticField()};
    }
    [[nodiscard]] std::shared_ptr<VField_t<T, Dim>> getFlippedZSlabField() {
        return {workspace_m, &requireWorkspace().flippedZSlabField()};
    }
    [[nodiscard]] std::shared_ptr<VField_t<T, Dim>> getOrCreateFlippedZSlabField(
            const VField_t<T, Dim>& source) {
        return {workspace_m, &requireWorkspace().mirrorScratchFor(source)};
    }

    void initializeFields(const std::string& solverType = "") {
        requireWorkspace().initializeFields(solverType);
    }

private:
    [[nodiscard]] Workspace& requireWorkspace() {
        if (!workspace_m) {
            throw OpalException("FieldContainer", "PIC workspace is not available.");
        }
        return *workspace_m;
    }

    [[nodiscard]] const Workspace& requireWorkspace() const {
        if (!workspace_m) {
            throw OpalException("FieldContainer", "PIC workspace is not available.");
        }
        return *workspace_m;
    }

    std::shared_ptr<Workspace> workspace_m;
};

#endif  // OPAL_FIELD_CONTAINER_H
