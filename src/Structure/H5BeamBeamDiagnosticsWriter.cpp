#include "Structure/H5BeamBeamDiagnosticsWriter.h"

#include "Utilities/OpalException.h"
#include "Utility/PAssert.h"

#include <filesystem>

namespace {
    constexpr char kWhere[] = "H5BeamBeamDiagnosticsWriter";
}

H5BeamBeamDiagnosticsWriter::H5BeamBeamDiagnosticsWriter(const std::string& fileName)
    : file_m(0), fileName_m(fileName), stepCounter_m(0), headerWritten_m(false) {
    open();
}

H5BeamBeamDiagnosticsWriter::~H5BeamBeamDiagnosticsWriter() { close(); }

void H5BeamBeamDiagnosticsWriter::close() {
    preparedStep_m.reset();
    if (file_m) {
        ippl::Comm->barrier();
        reportOnError(H5CloseFile(file_m), "close");
        file_m = 0;
    }
}

void H5BeamBeamDiagnosticsWriter::beginStep(
        const StepMetadata& meta, const Field_t<3>& rhoDensity) {
    PreparedStep step;
    step.meta       = meta;
    step.localIndex = rhoDensity.getLayout().getLocalNDIndex();
    step.rho        = flattenScalarField(rhoDensity, step.localIndex);
    preparedStep_m  = std::move(step);
}

void H5BeamBeamDiagnosticsWriter::endStep(
        const Field_t<3>& phiField, const VField_t<double, 3>& efield) {
    if (!preparedStep_m.has_value()) {
        throw OpalException(kWhere, "endStep() called without a matching beginStep().");
    }

    if (!headerWritten_m) {
        writeHeader();
    }

    const PreparedStep& step            = *preparedStep_m;
    const auto& meta                    = step.meta;
    const std::vector<h5_float64_t> phi = flattenScalarField(phiField, step.localIndex);

    reportOnError(H5SetStep(file_m, stepCounter_m), "H5SetStep");
    writeStepMetadata(meta);

    writeScalarField("rho", step.rho, step.localIndex, meta.origin, meta.spacing);
    writeScalarField("phi", phi, step.localIndex, meta.origin, meta.spacing);

    std::vector<h5_float64_t> ex;
    std::vector<h5_float64_t> ey;
    std::vector<h5_float64_t> ez;
    flattenVectorField(efield, step.localIndex, ex, ey, ez);

    writeScalarField("Ex", ex, step.localIndex, meta.origin, meta.spacing);
    writeScalarField("Ey", ey, step.localIndex, meta.origin, meta.spacing);
    writeScalarField("Ez", ez, step.localIndex, meta.origin, meta.spacing);

    ++stepCounter_m;
    preparedStep_m.reset();
}

void H5BeamBeamDiagnosticsWriter::reportOnError(h5_int64_t rc, const char* where) {
    if (rc == H5_SUCCESS) {
        return;
    }

    throw OpalException(kWhere, std::string(where) + " failed with H5 rc=" + std::to_string(rc));
}

std::vector<h5_float64_t> H5BeamBeamDiagnosticsWriter::flattenScalarField(
        const Field_t<3>& field, const ippl::NDIndex<3>& localIndex) {
    const int nghost = field.getNghost();
    auto fieldView   = field.getView();
    auto fieldHost   = field.getHostMirror();
    Kokkos::deep_copy(fieldHost, fieldView);
    ippl::Comm->barrier();
    const std::size_t nx = localIndex[0].length();
    const std::size_t ny = localIndex[1].length();
    const std::size_t nz = localIndex[2].length();

    std::vector<h5_float64_t> data;
    data.reserve(nx * ny * nz);

    for (int iz = localIndex[2].first(); iz <= localIndex[2].last(); ++iz) {
        for (int iy = localIndex[1].first(); iy <= localIndex[1].last(); ++iy) {
            for (int ix = localIndex[0].first(); ix <= localIndex[0].last(); ++ix) {
                data.push_back(fieldHost(ix + nghost, iy + nghost, iz + nghost));
            }
        }
    }

    return data;
}

void H5BeamBeamDiagnosticsWriter::flattenVectorField(
        const VField_t<double, 3>& field, const ippl::NDIndex<3>& localIndex,
        std::vector<h5_float64_t>& x, std::vector<h5_float64_t>& y, std::vector<h5_float64_t>& z) {
    const int nghost = field.getNghost();
    auto fieldView   = field.getView();
    auto fieldHost   = field.getHostMirror();
    Kokkos::deep_copy(fieldHost, fieldView);

    const std::size_t nx = localIndex[0].length();
    const std::size_t ny = localIndex[1].length();
    const std::size_t nz = localIndex[2].length();
    const std::size_t n  = nx * ny * nz;

    x.clear();
    y.clear();
    z.clear();
    x.reserve(n);
    y.reserve(n);
    z.reserve(n);

    for (int iz = localIndex[2].first(); iz <= localIndex[2].last(); ++iz) {
        for (int iy = localIndex[1].first(); iy <= localIndex[1].last(); ++iy) {
            for (int ix = localIndex[0].first(); ix <= localIndex[0].last(); ++ix) {
                const auto value = fieldHost(ix + nghost, iy + nghost, iz + nghost);
                x.push_back(value[0]);
                y.push_back(value[1]);
                z.push_back(value[2]);
            }
        }
    }
}

void H5BeamBeamDiagnosticsWriter::open() {
    if (file_m) {
        return;
    }

    if (ippl::Comm->rank() == 0) {
        std::filesystem::path outputPath(fileName_m);
        const std::filesystem::path parent = outputPath.parent_path();
        if (!parent.empty()) {
            std::filesystem::create_directories(parent);
        }
        if (std::filesystem::exists(outputPath)) {
            std::filesystem::remove(outputPath);
        }
    }
    ippl::Comm->barrier();

    h5_prop_t props = H5CreateFileProp();
    MPI_Comm comm   = ippl::Comm->getCommunicator();
    h5_err_t h5err  = H5SetPropFileMPIOCollective(props, &comm);
#if defined(NDEBUG)
    (void)h5err;
#endif
    PAssert(h5err != H5_ERR);

    file_m = H5OpenFile(fileName_m.c_str(), H5_O_WRONLY, props);
    PAssert(file_m != (h5_file_t)H5_ERR);
    H5CloseProp(props);
}

void H5BeamBeamDiagnosticsWriter::writeHeader() {
    reportOnError(H5WriteFileAttribString(file_m, "OPAL_flavour", "OPALX"), "OPAL_flavour");
    reportOnError(
            H5WriteFileAttribString(file_m, "diagnostics_kind", "beambeam_field_diagnostics"),
            "diagnostics_kind");
    headerWritten_m = true;
}

void H5BeamBeamDiagnosticsWriter::writeStepMetadata(const StepMetadata& meta) {
    const h5_int64_t globalStep = static_cast<h5_int64_t>(meta.globalStep);
    const h5_int64_t active     = meta.interactionWindowActive ? 1 : 0;

    reportOnError(H5WriteStepAttribInt64(file_m, "global_step", &globalStep, 1), "global_step");
    reportOnError(H5WriteStepAttribFloat64(file_m, "time", &meta.time, 1), "time");
    reportOnError(
            H5WriteStepAttribFloat64(file_m, "path_length_s", &meta.pathLengthS, 1),
            "path_length_s");
    reportOnError(
            H5WriteStepAttribInt64(file_m, "shape", meta.shape.data(), meta.shape.size()), "shape");
    reportOnError(
            H5WriteStepAttribFloat64(file_m, "mesh_origin", meta.origin.data(), meta.origin.size()),
            "mesh_origin");
    reportOnError(
            H5WriteStepAttribFloat64(
                    file_m, "mesh_spacing", meta.spacing.data(), meta.spacing.size()),
            "mesh_spacing");
    reportOnError(
            H5WriteStepAttribString(file_m, "coordinate_frame", meta.coordinateFrame.c_str()),
            "coordinate_frame");
    reportOnError(
            H5WriteStepAttribString(file_m, "snapshot_kind", meta.snapshotKind.c_str()),
            "snapshot_kind");
    reportOnError(
            H5WriteStepAttribInt64(file_m, "interaction_window_active", &active, 1),
            "interaction_window_active");
    reportOnError(
            H5WriteStepAttribFloat64(file_m, "interaction_point_s", &meta.interactionPointS, 1),
            "interaction_point_s");
    reportOnError(
            H5WriteStepAttribFloat64(
                    file_m, "ip_element_s_range", meta.beamBeamSRange.data(),
                    meta.beamBeamSRange.size()),
            "ip_element_s_range");
    reportOnError(
            H5WriteStepAttribFloat64(file_m, "particle_total_charge", &meta.particleTotalCharge, 1),
            "particle_total_charge");
    reportOnError(
            H5WriteStepAttribFloat64(
                    file_m, "particle_mean_r", meta.particleMeanR.data(),
                    meta.particleMeanR.size()),
            "particle_mean_r");
    reportOnError(
            H5WriteStepAttribFloat64(file_m, "particle_mean_s", &meta.particleMeanS, 1),
            "particle_mean_s");
}

void H5BeamBeamDiagnosticsWriter::writeScalarField(
        const std::string& name, const std::vector<h5_float64_t>& data,
        const ippl::NDIndex<3>& localIndex, const std::array<double, 3>& origin,
        const std::array<double, 3>& spacing) {
    reportOnError(
            H5Block3dSetView(
                    file_m, localIndex[0].min(), localIndex[0].max(), localIndex[1].min(),
                    localIndex[1].max(), localIndex[2].min(), localIndex[2].max()),
            "H5Block3dSetView");

    reportOnError(
            H5Block3dWriteScalarFieldFloat64(file_m, name.c_str(), data.data()), name.c_str());
    reportOnError(
            H5Block3dSetFieldOrigin(file_m, name.c_str(), origin[0], origin[1], origin[2]),
            "H5Block3dSetFieldOrigin");
    reportOnError(
            H5Block3dSetFieldSpacing(file_m, name.c_str(), spacing[0], spacing[1], spacing[2]),
            "H5Block3dSetFieldSpacing");
}
