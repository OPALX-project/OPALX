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
 * - scalar charge density, scalar potential, and electric-field components
 *   on the field mesh
 *
 * The charge density is captured before the solver coupling constant is applied,
 * while the scalar potential and electric field are written after the solve.
 * To preserve that timing without duplicating higher-level file logic in
 * PartBunch, the writer exposes a small two-stage API: beginStep() caches the
 * physical charge density, and endStep() writes the full step once the solve
 * outputs are available.
 */
class H5BeamBeamDiagnosticsWriter {
public:
    struct StepMetadata {
        long long globalStep = 0;
        double time          = 0.0;
        double pathLengthS   = 0.0;

        std::array<h5_int64_t, 3> shape = {0, 0, 0};
        std::array<double, 3> origin    = {0.0, 0.0, 0.0};
        std::array<double, 3> spacing   = {0.0, 0.0, 0.0};
        std::string coordinateFrame     = "beam";
        std::string snapshotKind        = "active_beambeam";

        bool interactionWindowActive         = false;
        double interactionPointS             = 0.0;
        double interactionPointLocalZ        = 0.0;
        std::array<double, 2> beamBeamSRange = {0.0, 0.0};

        long long particleTotalNum = 0;
        double particleTotalCharge = 0.0;
        /**
         * @brief Bunch centroid in the diagnostic coordinate frame [m].
         *
         * For the current BeamBeam diagnostics this frame stores the longitudinal
         * component as the physical laboratory/reference-coordinate @f$s@f$ of the
         * bunch center. It must therefore not be shifted again by pathLengthS when
         * producing particleMeanS.
         */
        std::array<double, 3> particleMeanR = {0.0, 0.0, 0.0};
        /**
         * @brief Physical longitudinal bunch centroid @f$\langle s\rangle@f$ [m].
         *
         * This is a position, not an offset from the reference particle. Legacy
         * diagnostics accidentally wrote @f$s_\mathrm{ref} + \langle r_z\rangle@f$
         * even though @f$\langle r_z\rangle@f$ was already in the physical
         * longitudinal coordinate. normalizeParticleMeanS() keeps writer output
         * consistent with this definition.
         */
        double particleMeanS = 0.0;

        std::string fieldStage = "";
    };

    explicit H5BeamBeamDiagnosticsWriter(const std::string& fileName);
    ~H5BeamBeamDiagnosticsWriter();

    void close();

    /**
     * @brief Return the physical @f$\langle s\rangle@f$ that should be written.
     *
     * BeamBeam diagnostics historically had one producer that filled
     * StepMetadata::particleMeanS as @f$s_\mathrm{ref} + \langle r_z\rangle@f$,
     * while StepMetadata::particleMeanR[2] was already the physical longitudinal
     * centroid. That makes the table appear to advance at roughly
     * @f$2\beta c@f$. The writer fixes only that recognizable legacy value and
     * otherwise preserves the caller-provided physical centroid.
     */
    static double normalizeParticleMeanS(const StepMetadata& meta);

    static std::string makeOutputFileName(
            const std::string& basename, const std::string& what, const std::string& type,
            const std::string& tag);

    static void writeScalarQuantity(
            const std::string& fileName, const std::string& fieldName, const StepMetadata& meta,
            const Field_t<3>& field);

    static void writeVectorQuantity(
            const std::string& fileName, const std::string& fieldName, const StepMetadata& meta,
            const VField_t<double, 3>& field);

    void beginStep(const StepMetadata& meta, const Field_t<3>& rhoDensity);
    void endStep(const Field_t<3>& phiField, const VField_t<double, 3>& efield);

private:
    struct PreparedStep {
        StepMetadata meta;
        ippl::NDIndex<3> localIndex;
        std::vector<h5_float64_t> rho;
    };

    static void reportOnError(h5_int64_t rc, const char* where);
    static std::vector<h5_float64_t> flattenScalarField(
            const Field_t<3>& field, const ippl::NDIndex<3>& localIndex);
    static void flattenVectorField(
            const VField_t<double, 3>& field, const ippl::NDIndex<3>& localIndex,
            std::vector<h5_float64_t>& x, std::vector<h5_float64_t>& y,
            std::vector<h5_float64_t>& z);

    static void writeHeader(h5_file_t file);
    static void writeStepMetadata(h5_file_t file, const StepMetadata& meta);
    static void writeScalarField(
            h5_file_t file, const std::string& name, const std::vector<h5_float64_t>& data,
            const ippl::NDIndex<3>& localIndex, const std::array<double, 3>& origin,
            const std::array<double, 3>& spacing);
    static void writeVectorField(
            h5_file_t file, const std::string& name, const std::vector<h5_float64_t>& x,
            const std::vector<h5_float64_t>& y, const std::vector<h5_float64_t>& z,
            const ippl::NDIndex<3>& localIndex, const std::array<double, 3>& origin,
            const std::array<double, 3>& spacing);

    void open();
    void writeHeader();
    void writeStepMetadata(const StepMetadata& meta);
    void writeScalarField(
            const std::string& name, const std::vector<h5_float64_t>& data,
            const ippl::NDIndex<3>& localIndex, const std::array<double, 3>& origin,
            const std::array<double, 3>& spacing);

    h5_file_t file_m;
    std::string fileName_m;
    h5_ssize_t stepCounter_m;
    bool headerWritten_m;
    std::optional<PreparedStep> preparedStep_m;
};

#endif
