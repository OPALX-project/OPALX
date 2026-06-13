#include "Structure/H5BeamBeamDiagnosticsWriter.h"

#include "Utilities/OpalException.h"
#include "Utility/PAssert.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <limits>
#include <set>

namespace {
    constexpr char kWhere[] = "H5BeamBeamDiagnosticsWriter";

    std::set<std::string>& initializedQuantityFiles() {
        static std::set<std::string> files;
        return files;
    }

    h5_file_t openQuantityFile(const std::string& fileName, bool& isFirstWrite) {
        auto& initialized = initializedQuantityFiles();
        isFirstWrite      = initialized.find(fileName) == initialized.end();

        if (ippl::Comm->rank() == 0) {
            std::filesystem::path outputPath(fileName);
            const std::filesystem::path parent = outputPath.parent_path();
            if (!parent.empty()) {
                std::filesystem::create_directories(parent);
            }
            if (isFirstWrite && std::filesystem::exists(outputPath)) {
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

        const h5_int32_t flags = isFirstWrite ? H5_O_WRONLY : H5_O_APPENDONLY;
        h5_file_t file         = H5OpenFile(fileName.c_str(), flags, props);
        PAssert(file != (h5_file_t)H5_ERR);
        H5CloseProp(props);

        initialized.insert(fileName);
        return file;
    }
}  // namespace

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

double H5BeamBeamDiagnosticsWriter::normalizeParticleMeanS(const StepMetadata& meta) {
    const double rzPhysicalS = meta.particleMeanR[2];
    const double legacyValue = meta.pathLengthS + rzPhysicalS;

    if (!std::isfinite(meta.particleMeanS) || !std::isfinite(meta.pathLengthS)
        || !std::isfinite(rzPhysicalS)) {
        return meta.particleMeanS;
    }

    const double scale = std::max(
            {1.0, std::abs(meta.particleMeanS), std::abs(legacyValue), std::abs(rzPhysicalS)});
    const double tolerance = 64.0 * std::numeric_limits<double>::epsilon() * scale;

    if (std::abs(meta.particleMeanS - legacyValue) <= tolerance) {
        return rzPhysicalS;
    }

    return meta.particleMeanS;
}

std::string H5BeamBeamDiagnosticsWriter::makeOutputFileName(
        const std::string& basename, const std::string& what, const std::string& type,
        const std::string& tag) {
    std::filesystem::path file("data");
    std::string name = basename + "-" + what + "_" + type;
    if (!tag.empty()) {
        name += "-" + tag;
    }
    name += ".h5";
    file /= name;
    return file.string();
}

void H5BeamBeamDiagnosticsWriter::writeScalarQuantity(
        const std::string& fileName, const std::string& fieldName, const StepMetadata& meta,
        const Field_t<3>& field) {
    bool isFirstWrite = false;
    h5_file_t file    = openQuantityFile(fileName, isFirstWrite);

    if (isFirstWrite) {
        writeHeader(file);
    }

    const h5_ssize_t step = isFirstWrite ? 0 : H5GetNumSteps(file);
    reportOnError(H5SetStep(file, step), "H5SetStep");
    writeStepMetadata(file, meta);

    const ippl::NDIndex<3> localIndex = field.getLayout().getLocalNDIndex();
    const auto data                   = flattenScalarField(field, localIndex);
    writeScalarField(file, fieldName, data, localIndex, meta.origin, meta.spacing);

    ippl::Comm->barrier();
    reportOnError(H5CloseFile(file), "close scalar quantity");
}

void H5BeamBeamDiagnosticsWriter::writeVectorQuantity(
        const std::string& fileName, const std::string& fieldName, const StepMetadata& meta,
        const VField_t<double, 3>& field) {
    bool isFirstWrite = false;
    h5_file_t file    = openQuantityFile(fileName, isFirstWrite);

    if (isFirstWrite) {
        writeHeader(file);
    }

    const h5_ssize_t step = isFirstWrite ? 0 : H5GetNumSteps(file);
    reportOnError(H5SetStep(file, step), "H5SetStep");
    writeStepMetadata(file, meta);

    const ippl::NDIndex<3> localIndex = field.getLayout().getLocalNDIndex();
    std::vector<h5_float64_t> x;
    std::vector<h5_float64_t> y;
    std::vector<h5_float64_t> z;
    flattenVectorField(field, localIndex, x, y, z);
    writeVectorField(file, fieldName, x, y, z, localIndex, meta.origin, meta.spacing);

    ippl::Comm->barrier();
    reportOnError(H5CloseFile(file), "close vector quantity");
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
    const int nghost     = field.getNghost();
    auto fieldView       = field.getView();
    const std::size_t nx = localIndex[0].length();
    const std::size_t ny = localIndex[1].length();
    const std::size_t nz = localIndex[2].length();
    const std::size_t n  = nx * ny * nz;

    Kokkos::View<h5_float64_t*> dataDevice("H5BeamBeamDiagnosticsWriter::dataDevice", n);
    Kokkos::parallel_for(
            "H5BeamBeamDiagnosticsWriter::flattenScalarField",
            Kokkos::MDRangePolicy<Kokkos::Rank<3>>(
                    {0, 0, 0}, {static_cast<int>(nz), static_cast<int>(ny), static_cast<int>(nx)}),
            KOKKOS_LAMBDA(const int izLocal, const int iyLocal, const int ixLocal) {
                const std::size_t flatIndex =
                        (static_cast<std::size_t>(izLocal) * ny + static_cast<std::size_t>(iyLocal))
                                * nx
                        + static_cast<std::size_t>(ixLocal);
                // localIndex describes the global H5 block; the Kokkos field view is local and
                // therefore must be accessed in local coordinates plus ghost offset only.
                dataDevice(flatIndex) =
                        fieldView(ixLocal + nghost, iyLocal + nghost, izLocal + nghost);
            });

    auto dataHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), dataDevice);

    std::vector<h5_float64_t> data;
    data.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        data.push_back(dataHost(i));
    }

    return data;
}

void H5BeamBeamDiagnosticsWriter::flattenVectorField(
        const VField_t<double, 3>& field, const ippl::NDIndex<3>& localIndex,
        std::vector<h5_float64_t>& x, std::vector<h5_float64_t>& y, std::vector<h5_float64_t>& z) {
    const int nghost = field.getNghost();
    auto fieldView   = field.getView();

    const std::size_t nx = localIndex[0].length();
    const std::size_t ny = localIndex[1].length();
    const std::size_t nz = localIndex[2].length();
    const std::size_t n  = nx * ny * nz;

    Kokkos::View<h5_float64_t*> xDevice("H5BeamBeamDiagnosticsWriter::xDevice", n);
    Kokkos::View<h5_float64_t*> yDevice("H5BeamBeamDiagnosticsWriter::yDevice", n);
    Kokkos::View<h5_float64_t*> zDevice("H5BeamBeamDiagnosticsWriter::zDevice", n);

    Kokkos::parallel_for(
            "H5BeamBeamDiagnosticsWriter::flattenVectorField",
            Kokkos::MDRangePolicy<Kokkos::Rank<3>>(
                    {0, 0, 0}, {static_cast<int>(nz), static_cast<int>(ny), static_cast<int>(nx)}),
            KOKKOS_LAMBDA(const int izLocal, const int iyLocal, const int ixLocal) {
                const std::size_t flatIndex =
                        (static_cast<std::size_t>(izLocal) * ny + static_cast<std::size_t>(iyLocal))
                                * nx
                        + static_cast<std::size_t>(ixLocal);
                const auto value = fieldView(ixLocal + nghost, iyLocal + nghost, izLocal + nghost);
                xDevice(flatIndex) = value[0];
                yDevice(flatIndex) = value[1];
                zDevice(flatIndex) = value[2];
            });

    auto xHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), xDevice);
    auto yHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), yDevice);
    auto zHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), zDevice);

    x.clear();
    y.clear();
    z.clear();
    x.reserve(n);
    y.reserve(n);
    z.reserve(n);

    for (std::size_t i = 0; i < n; ++i) {
        x.push_back(xHost(i));
        y.push_back(yHost(i));
        z.push_back(zHost(i));
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
    writeHeader(file_m);
    headerWritten_m = true;
}

void H5BeamBeamDiagnosticsWriter::writeHeader(h5_file_t file) {
    reportOnError(H5WriteFileAttribString(file, "OPAL_flavour", "OPALX"), "OPAL_flavour");
    reportOnError(
            H5WriteFileAttribString(file, "diagnostics_kind", "beambeam_field_diagnostics"),
            "diagnostics_kind");
}

void H5BeamBeamDiagnosticsWriter::writeStepMetadata(const StepMetadata& meta) {
    writeStepMetadata(file_m, meta);
}

void H5BeamBeamDiagnosticsWriter::writeStepMetadata(h5_file_t file, const StepMetadata& meta) {
    const h5_int64_t globalStep       = static_cast<h5_int64_t>(meta.globalStep);
    const h5_int64_t active           = meta.interactionWindowActive ? 1 : 0;
    const h5_int64_t particleTotalNum = static_cast<h5_int64_t>(meta.particleTotalNum);
    const double particleMeanS        = normalizeParticleMeanS(meta);

    reportOnError(H5WriteStepAttribInt64(file, "global_step", &globalStep, 1), "global_step");
    reportOnError(H5WriteStepAttribFloat64(file, "time", &meta.time, 1), "time");
    reportOnError(
            H5WriteStepAttribFloat64(file, "path_length_s", &meta.pathLengthS, 1), "path_length_s");
    reportOnError(
            H5WriteStepAttribInt64(file, "shape", meta.shape.data(), meta.shape.size()), "shape");
    reportOnError(
            H5WriteStepAttribFloat64(file, "mesh_origin", meta.origin.data(), meta.origin.size()),
            "mesh_origin");
    reportOnError(
            H5WriteStepAttribFloat64(
                    file, "mesh_spacing", meta.spacing.data(), meta.spacing.size()),
            "mesh_spacing");
    reportOnError(
            H5WriteStepAttribString(file, "coordinate_frame", meta.coordinateFrame.c_str()),
            "coordinate_frame");
    reportOnError(
            H5WriteStepAttribString(file, "snapshot_kind", meta.snapshotKind.c_str()),
            "snapshot_kind");
    reportOnError(
            H5WriteStepAttribInt64(file, "interaction_window_active", &active, 1),
            "interaction_window_active");
    reportOnError(
            H5WriteStepAttribFloat64(file, "interaction_point_s", &meta.interactionPointS, 1),
            "interaction_point_s");
    reportOnError(
            H5WriteStepAttribFloat64(
                    file, "interaction_point_local_z", &meta.interactionPointLocalZ, 1),
            "interaction_point_local_z");
    reportOnError(
            H5WriteStepAttribFloat64(
                    file, "ip_element_s_range", meta.beamBeamSRange.data(),
                    meta.beamBeamSRange.size()),
            "ip_element_s_range");
    reportOnError(
            H5WriteStepAttribInt64(file, "particle_total_num", &particleTotalNum, 1),
            "particle_total_num");
    reportOnError(
            H5WriteStepAttribFloat64(file, "particle_total_charge", &meta.particleTotalCharge, 1),
            "particle_total_charge");
    reportOnError(
            H5WriteStepAttribFloat64(
                    file, "particle_mean_r", meta.particleMeanR.data(), meta.particleMeanR.size()),
            "particle_mean_r");
    reportOnError(
            H5WriteStepAttribFloat64(file, "particle_mean_s", &particleMeanS, 1),
            "particle_mean_s");
    if (!meta.fieldStage.empty()) {
        reportOnError(
                H5WriteStepAttribString(file, "field_stage", meta.fieldStage.c_str()),
                "field_stage");
    }
}

void H5BeamBeamDiagnosticsWriter::writeScalarField(
        const std::string& name, const std::vector<h5_float64_t>& data,
        const ippl::NDIndex<3>& localIndex, const std::array<double, 3>& origin,
        const std::array<double, 3>& spacing) {
    writeScalarField(file_m, name, data, localIndex, origin, spacing);
}

void H5BeamBeamDiagnosticsWriter::writeScalarField(
        h5_file_t file, const std::string& name, const std::vector<h5_float64_t>& data,
        const ippl::NDIndex<3>& localIndex, const std::array<double, 3>& origin,
        const std::array<double, 3>& spacing) {
    reportOnError(
            H5Block3dSetView(
                    file, localIndex[0].min(), localIndex[0].max(), localIndex[1].min(),
                    localIndex[1].max(), localIndex[2].min(), localIndex[2].max()),
            "H5Block3dSetView");

    reportOnError(H5Block3dWriteScalarFieldFloat64(file, name.c_str(), data.data()), name.c_str());
    reportOnError(
            H5Block3dSetFieldOrigin(file, name.c_str(), origin[0], origin[1], origin[2]),
            "H5Block3dSetFieldOrigin");
    reportOnError(
            H5Block3dSetFieldSpacing(file, name.c_str(), spacing[0], spacing[1], spacing[2]),
            "H5Block3dSetFieldSpacing");
}

void H5BeamBeamDiagnosticsWriter::writeVectorField(
        h5_file_t file, const std::string& name, const std::vector<h5_float64_t>& x,
        const std::vector<h5_float64_t>& y, const std::vector<h5_float64_t>& z,
        const ippl::NDIndex<3>& localIndex, const std::array<double, 3>& origin,
        const std::array<double, 3>& spacing) {
    reportOnError(
            H5Block3dSetView(
                    file, localIndex[0].min(), localIndex[0].max(), localIndex[1].min(),
                    localIndex[1].max(), localIndex[2].min(), localIndex[2].max()),
            "H5Block3dSetView");

    reportOnError(
            H5Block3dWriteVector3dFieldFloat64(file, name.c_str(), x.data(), y.data(), z.data()),
            name.c_str());
    reportOnError(
            H5Block3dSetFieldOrigin(file, name.c_str(), origin[0], origin[1], origin[2]),
            "H5Block3dSetFieldOrigin");
    reportOnError(
            H5Block3dSetFieldSpacing(file, name.c_str(), spacing[0], spacing[1], spacing[2]),
            "H5Block3dSetFieldSpacing");
}
