#include "Structure/H5BeamBeamDiagnosticsWriter.h"

#include "Utilities/OpalException.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <limits>
#include <set>
#include <utility>

namespace {
    constexpr char kWhere[] = "H5BeamBeamDiagnosticsWriter";

    // HDF5 operations are collective for this writer. Make reported failures
    // uniform before unwinding handles, so peers do not enter the next write.
    void requireCollectively(bool success, const std::string& operation) {
        int failed    = success ? 0 : 1;
        int anyFailed = 0;
        MPI_Allreduce(&failed, &anyFailed, 1, MPI_INT, MPI_MAX, ippl::Comm->getCommunicator());
        if (anyFailed) {
            throw OpalException(kWhere, operation + " failed on at least one MPI rank.");
        }
    }

    template <typename Handle, auto Close>
    class ScopedH5Handle {
    public:
        explicit ScopedH5Handle(Handle value) : value_m(value) {}
        ScopedH5Handle(const ScopedH5Handle&)            = delete;
        ScopedH5Handle& operator=(const ScopedH5Handle&) = delete;
        ~ScopedH5Handle() noexcept {
            if (close() != H5_SUCCESS) {
                std::fputs("H5BeamBeamDiagnosticsWriter: handle cleanup failed.\n", stderr);
            }
        }
        Handle get() const { return value_m; }
        Handle release() { return std::exchange(value_m, 0); }
        h5_err_t close() {
            const auto value = release();
            return value ? Close(value) : H5_SUCCESS;
        }

    private:
        Handle value_m;
    };

    using ScopedFile       = ScopedH5Handle<h5_file_t, H5CloseFile>;
    using ScopedProperties = ScopedH5Handle<h5_prop_t, H5CloseProp>;

    h5_file_t openDiagnosticFile(const std::string& fileName, bool truncate) {
        bool prepared = true;
        if (ippl::Comm->rank() == 0) {
            const std::filesystem::path outputPath(fileName);
            const auto parent = outputPath.parent_path();
            std::error_code error;
            if (!parent.empty()) {
                std::filesystem::create_directories(parent, error);
            }
            if (!error && truncate) {
                std::filesystem::remove(outputPath, error);
            }
            prepared = !error;
        }
        requireCollectively(prepared, "Preparing diagnostic output " + fileName);

        const h5_prop_t props = H5CreateFileProp();
        const bool validProps = static_cast<std::intptr_t>(props) > 0;
        ScopedProperties properties(validProps ? props : 0);
        requireCollectively(validProps, "H5CreateFileProp");
        MPI_Comm comm = ippl::Comm->getCommunicator();
        requireCollectively(
                H5SetPropFileMPIOCollective(props, &comm) == H5_SUCCESS,
                "H5SetPropFileMPIOCollective");

        const h5_file_t file =
                H5OpenFile(fileName.c_str(), truncate ? H5_O_WRONLY : H5_O_APPENDONLY, props);
        const bool validFile = static_cast<std::intptr_t>(file) > 0;
        ScopedFile opened(validFile ? file : 0);
        requireCollectively(validFile, "H5OpenFile " + fileName);
        requireCollectively(properties.close() == H5_SUCCESS, "H5CloseProp");
        return opened.release();
    }

    std::set<std::string>& initializedQuantityFiles() {
        static std::set<std::string> files;
        return files;
    }

    h5_file_t openQuantityFile(const std::string& fileName, bool& isFirstWrite) {
        auto& initialized = initializedQuantityFiles();
        isFirstWrite      = initialized.find(fileName) == initialized.end();

        return openDiagnosticFile(fileName, isFirstWrite);
    }

    // Keep the Kokkos flatten kernels outside H5BeamBeamDiagnosticsWriter's
    // private static member functions.  With CUDA, KOKKOS_LAMBDA uses extended
    // __host__ __device__ lambdas, and NVCC rejects those lambdas when the
    // enclosing member function has private or protected access.
    Kokkos::View<h5_float64_t*> flattenScalarFieldDevice(
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
                        {0, 0, 0},
                        {static_cast<int>(nz), static_cast<int>(ny), static_cast<int>(nx)}),
                KOKKOS_LAMBDA(const int izLocal, const int iyLocal, const int ixLocal) {
                    const std::size_t flatIndex = (static_cast<std::size_t>(izLocal) * ny
                                                   + static_cast<std::size_t>(iyLocal))
                                                          * nx
                                                  + static_cast<std::size_t>(ixLocal);
                    // localIndex describes the global H5 block; the Kokkos field view is local
                    // and therefore must be accessed in local coordinates plus ghost offset only.
                    dataDevice(flatIndex) =
                            fieldView(ixLocal + nghost, iyLocal + nghost, izLocal + nghost);
                });

        return dataDevice;
    }

    void flattenVectorFieldDevice(
            const VField_t<double, 3>& field, const ippl::NDIndex<3>& localIndex,
            Kokkos::View<h5_float64_t*>& xDevice, Kokkos::View<h5_float64_t*>& yDevice,
            Kokkos::View<h5_float64_t*>& zDevice) {
        const int nghost = field.getNghost();
        auto fieldView   = field.getView();

        const std::size_t nx = localIndex[0].length();
        const std::size_t ny = localIndex[1].length();
        const std::size_t nz = localIndex[2].length();

        Kokkos::parallel_for(
                "H5BeamBeamDiagnosticsWriter::flattenVectorField",
                Kokkos::MDRangePolicy<Kokkos::Rank<3>>(
                        {0, 0, 0},
                        {static_cast<int>(nz), static_cast<int>(ny), static_cast<int>(nx)}),
                KOKKOS_LAMBDA(const int izLocal, const int iyLocal, const int ixLocal) {
                    const std::size_t flatIndex = (static_cast<std::size_t>(izLocal) * ny
                                                   + static_cast<std::size_t>(iyLocal))
                                                          * nx
                                                  + static_cast<std::size_t>(ixLocal);
                    const auto value =
                            fieldView(ixLocal + nghost, iyLocal + nghost, izLocal + nghost);
                    xDevice(flatIndex) = value[0];
                    yDevice(flatIndex) = value[1];
                    zDevice(flatIndex) = value[2];
                });
    }
}  // namespace

H5BeamBeamDiagnosticsWriter::H5BeamBeamDiagnosticsWriter(const std::string& fileName)
    : file_m(0), fileName_m(fileName), stepCounter_m(0), headerWritten_m(false) {
    open();
}

H5BeamBeamDiagnosticsWriter::~H5BeamBeamDiagnosticsWriter() {
    // All ranks must release the writer before MPI finalization. Explicit close()
    // reports errors; destruction only performs best-effort nonthrowing cleanup.
    ScopedFile cleanup(std::exchange(file_m, 0));
}

void H5BeamBeamDiagnosticsWriter::close() {
    preparedStep_m.reset();
    if (file_m) {
        ippl::Comm->barrier();
        const auto file = std::exchange(file_m, 0);
        reportOnError(H5CloseFile(file), "close");
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
    ScopedFile handle(openQuantityFile(fileName, isFirstWrite));
    const h5_file_t file = handle.get();

    if (isFirstWrite) {
        writeHeader(file);
    }

    const h5_ssize_t step = isFirstWrite ? 0 : H5GetNumSteps(file);
    requireCollectively(step >= 0, "H5GetNumSteps");
    reportOnError(H5SetStep(file, step), "H5SetStep");
    writeStepMetadata(file, meta);

    const ippl::NDIndex<3> localIndex = field.getLayout().getLocalNDIndex();
    const auto data                   = flattenScalarField(field, localIndex);
    writeScalarField(file, fieldName, data, localIndex, meta.origin, meta.spacing);

    ippl::Comm->barrier();
    reportOnError(handle.close(), "close scalar quantity");
    initializedQuantityFiles().insert(fileName);
}

void H5BeamBeamDiagnosticsWriter::writeVectorQuantity(
        const std::string& fileName, const std::string& fieldName, const StepMetadata& meta,
        const VField_t<double, 3>& field) {
    bool isFirstWrite = false;
    ScopedFile handle(openQuantityFile(fileName, isFirstWrite));
    const h5_file_t file = handle.get();

    if (isFirstWrite) {
        writeHeader(file);
    }

    const h5_ssize_t step = isFirstWrite ? 0 : H5GetNumSteps(file);
    requireCollectively(step >= 0, "H5GetNumSteps");
    reportOnError(H5SetStep(file, step), "H5SetStep");
    writeStepMetadata(file, meta);

    const ippl::NDIndex<3> localIndex = field.getLayout().getLocalNDIndex();
    std::vector<h5_float64_t> x;
    std::vector<h5_float64_t> y;
    std::vector<h5_float64_t> z;
    flattenVectorField(field, localIndex, x, y, z);
    writeVectorField(file, fieldName, x, y, z, localIndex, meta.origin, meta.spacing);

    ippl::Comm->barrier();
    reportOnError(handle.close(), "close vector quantity");
    initializedQuantityFiles().insert(fileName);
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
    requireCollectively(rc == H5_SUCCESS, where);
}

std::vector<h5_float64_t> H5BeamBeamDiagnosticsWriter::flattenScalarField(
        const Field_t<3>& field, const ippl::NDIndex<3>& localIndex) {
    const std::size_t nx = localIndex[0].length();
    const std::size_t ny = localIndex[1].length();
    const std::size_t nz = localIndex[2].length();
    const std::size_t n  = nx * ny * nz;

    Kokkos::View<h5_float64_t*> dataDevice = flattenScalarFieldDevice(field, localIndex);
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
    const std::size_t nx = localIndex[0].length();
    const std::size_t ny = localIndex[1].length();
    const std::size_t nz = localIndex[2].length();
    const std::size_t n  = nx * ny * nz;

    Kokkos::View<h5_float64_t*> xDevice("H5BeamBeamDiagnosticsWriter::xDevice", n);
    Kokkos::View<h5_float64_t*> yDevice("H5BeamBeamDiagnosticsWriter::yDevice", n);
    Kokkos::View<h5_float64_t*> zDevice("H5BeamBeamDiagnosticsWriter::zDevice", n);

    flattenVectorFieldDevice(field, localIndex, xDevice, yDevice, zDevice);

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

    file_m = openDiagnosticFile(fileName_m, true);
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
