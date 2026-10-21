#ifndef OPAL_H5_BEAMBEAM_DIAGNOSTICS_WRITER_H
#define OPAL_H5_BEAMBEAM_DIAGNOSTICS_WRITER_H

#include "Manager/datatypes.h"

#include "H5hut.h"

#include <array>
#include <optional>
#include <string>
#include <vector>

/**
 * @brief Collective HDF5 output of caller-supplied BeamBeam mesh diagnostics.
 *
 * The production BeamBeam path uses writeScalarQuantity(): one file per quantity
 * and snapshot kind, with a new HDF5 step per call. The first successful write
 * to a filename in this process replaces existing output; later writes append.
 * Use a distinct filename for a distinct run within the same process.
 * Existing append targets are checked read-only on rank zero before collective
 * H5hut open, avoiding an upstream failed-open reporting bug for corrupt files.
 * Rejected files are not replaced or truncated. This check is not a guarantee
 * against concurrent external changes to the file between validation and open.
 *
 * Values are written verbatim: this writer applies no field scaling or coordinate
 * transformation. Current BeamBeam rho is the last bin's normalized primary
 * deposition before the Poisson solve, not physical C/m^3 and not the virtual
 * copy. ParticleMeshFieldTransfer applies couplingConstant/(dt*gamma), with
 * additional cell-volume normalization where required by the backend. The
 * producer records this convention in StepMetadata::fieldStage.
 *
 * The separate instance API owns one file: beginStep() caches the supplied rho,
 * then endStep() writes rho, phi and Ex/Ey/Ez. Its caller is responsible for
 * supplying a consistent frame and documenting the units/scaling of each field.
 *
 * File open/write/close operations are collective on ippl::Comm; every rank must
 * call in the same order with matching filenames/metadata. Reported I/O failures
 * are synchronized before throwing OpalException. Destroy/close all writers on
 * all ranks before MPI finalization. Explicit close() reports failures; destructor
 * cleanup is nonthrowing and reports failures to stderr. Not thread-safe.
 */
class H5BeamBeamDiagnosticsWriter {
public:
    /// Metadata copied without frame conversion; coordinates below are in metres.
    struct StepMetadata {
        /// Global tracking-step index (not the sequential HDF5 step index).
        long long globalStep = 0;
        /// Simulation time [s].
        double time = 0.0;
        /// Reference-particle path coordinate s [m].
        double pathLengthS = 0.0;

        /// Global owned-grid dimensions, excluding halos.
        std::array<h5_int64_t, 3> shape = {0, 0, 0};
        /// Mesh origin in coordinateFrame [m].
        std::array<double, 3> origin = {0.0, 0.0, 0.0};
        /// Mesh spacing in coordinateFrame [m].
        std::array<double, 3> spacing = {0.0, 0.0, 0.0};
        /// Producer-defined axes label; the writer does not transform values.
        std::string coordinateFrame = "beam";
        /// Producer-defined event label, e.g. before/after window enlargement.
        std::string snapshotKind = "active_beambeam";

        /// Whether the collective interaction window is active at this event.
        bool interactionWindowActive = false;
        /// Interaction point in the reference path coordinate s [m].
        double interactionPointS = 0.0;
        /// Interaction point in the diagnostic mesh's local z coordinate [m].
        double interactionPointLocalZ = 0.0;
        /// Geometric element start/end in reference s [m], not the fixed solve window.
        std::array<double, 2> beamBeamSRange = {0.0, 0.0};

        /// Global number of physical primary macroparticles (excluding virtual copy/witnesses).
        long long particleTotalNum = 0;
        /// Signed physical primary charge [C], excluding virtual copy/witnesses.
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

        /// Producer-defined normalization/solve-stage label for interpreting field values.
        std::string fieldStage = "";
    };

    /// Collectively create/truncate fileName, creating missing parent directories.
    /// @param fileName Output path shared by every participating rank.
    explicit H5BeamBeamDiagnosticsWriter(const std::string& fileName);
    /// Collectively release an open file without throwing; prefer explicit close().
    ~H5BeamBeamDiagnosticsWriter();

    /// Collectively close the file and discard pending data; idempotent after closing.
    /// @throws OpalException On a reported HDF5 close error.
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
     * @param meta Producer metadata, including the reference path and bunch centroid.
     * @return Corrected physical s [m], or the original value if no legacy match exists.
     */
    static double normalizeParticleMeanS(const StepMetadata& meta);

    /**
     * @brief Form data/basename-what_type[-tag].h5 without opening it.
     * @param basename Input/run basename.
     * @param what Quantity label, e.g. RHO.
     * @param type Scalar/vector label.
     * @param tag Optional event suffix; empty omits it.
     * @return Relative output path.
     */
    static std::string makeOutputFileName(
            const std::string& basename, const std::string& what, const std::string& type,
            const std::string& tag);

    /**
     * @brief Collectively append one scalar grid and metadata, then close the file.
     * @param fileName Output path, consistently used on every rank.
     * @param fieldName HDF5 field name.
     * @param meta Frame, mesh and interpretation metadata.
     * @param field Distributed scalar grid; only owned cells are written, without scaling.
     * @throws OpalException On reported file preparation or HDF5 errors.
     */
    static void writeScalarQuantity(
            const std::string& fileName, const std::string& fieldName, const StepMetadata& meta,
            const Field_t<3>& field);

    /** @copydoc writeScalarQuantity
     * @details Writes the three components of a distributed vector grid.
     */
    static void writeVectorQuantity(
            const std::string& fileName, const std::string& fieldName, const StepMetadata& meta,
            const VField_t<double, 3>& field);

    /**
     * @brief Cache owned rho cells and metadata for the next endStep(), replacing any pending step.
     * @param meta Frame/mesh metadata shared by rho, phi and E.
     * @param rhoDensity Caller-normalized rho; no unit conversion is performed.
     */
    void beginStep(const StepMetadata& meta, const Field_t<3>& rhoDensity);
    /**
     * @brief Collectively write the cached rho and current solve outputs to the instance file.
     * @param phiField Potential on the same mesh/layout as the cached rho.
     * @param efield Electric field on that mesh/layout, in the caller's convention.
     * @throws OpalException If beginStep() was not called or HDF5 reports an error.
     */
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
