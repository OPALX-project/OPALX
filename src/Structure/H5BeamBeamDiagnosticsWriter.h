#ifndef OPAL_H5_BEAMBEAM_DIAGNOSTICS_WRITER_H
#define OPAL_H5_BEAMBEAM_DIAGNOSTICS_WRITER_H

#include "PartBunch/FieldContainer.hpp"

#include "H5hut.h"

#include <array>
#include <optional>
#include <string>
#include <vector>

/**
 * @brief Write BeamBeam field diagnostics to a single HDF5 file.
 *
 * This writer is intentionally narrow in scope:
 * - one file per run
 * - one HDF5 step per BeamBeam-triggered diagnostic event
 * - scalar charge density and electric-field components on the field mesh
 *
 * The charge density is captured before the solver coupling constant is applied,
 * while the electric field is written after the solve. To preserve that timing
 * without duplicating higher-level file logic in PartBunch, the writer exposes
 * a small two-stage API: beginStep() caches the physical charge density, and
 * endStep() writes the full step once the electric field is available.
 */
class H5BeamBeamDiagnosticsWriter {
public:
    struct StepMetadata {
        long long globalStep = 0;
        double time          = 0.0;
        double pathLengthS   = 0.0;

        std::array<h5_int64_t, 3> shape   = {0, 0, 0};
        std::array<double, 3> origin      = {0.0, 0.0, 0.0};
        std::array<double, 3> spacing     = {0.0, 0.0, 0.0};
        std::string coordinateFrame       = "beam";
        std::string snapshotKind          = "active_beambeam";

        bool interactionWindowActive      = false;
        double interactionPointS          = 0.0;
        std::array<double, 2> beamBeamSRange = {0.0, 0.0};

        double particleTotalCharge        = 0.0;
        std::array<double, 3> particleMeanR = {0.0, 0.0, 0.0};
        double particleMeanS              = 0.0;
        double bunchSRef                  = 0.0;
        double bunchTailS                 = 0.0;
        double bunchHeadS                 = 0.0;
        double beamBeamWindowEndS         = 0.0;
        h5_int64_t leavingBeamBeamWindow  = 0;
    };

    explicit H5BeamBeamDiagnosticsWriter(const std::string& fileName);
    ~H5BeamBeamDiagnosticsWriter();

    void close();

    void beginStep(const StepMetadata& meta, const Field_t<3>& rhoDensity);
    void endStep(const VField_t<double, 3>& efield);

private:
    struct PreparedStep {
        StepMetadata meta;
        ippl::NDIndex<3> localIndex;
        std::vector<h5_float64_t> rho;
    };

    static void reportOnError(h5_int64_t rc, const char* where);
    static std::vector<h5_float64_t> flattenScalarField(
        const Field_t<3>& field,
        const ippl::NDIndex<3>& localIndex);
    static void flattenVectorField(
        const VField_t<double, 3>& field,
        const ippl::NDIndex<3>& localIndex,
        std::vector<h5_float64_t>& x,
        std::vector<h5_float64_t>& y,
        std::vector<h5_float64_t>& z);

    void open();
    void writeHeader();
    void writeStepMetadata(const StepMetadata& meta);
    void writeScalarField(
        const std::string& name,
        const std::vector<h5_float64_t>& data,
        const ippl::NDIndex<3>& localIndex,
        const std::array<double, 3>& origin,
        const std::array<double, 3>& spacing);

    h5_file_t file_m;
    std::string fileName_m;
    h5_ssize_t stepCounter_m;
    bool headerWritten_m;
    std::optional<PreparedStep> preparedStep_m;
};

#endif
