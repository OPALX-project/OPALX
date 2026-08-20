/**
 * @file EmittedFromFile.h
 * @brief Defines emitted file distributions with legacy and explicit birth-time semantics.
 */

#ifndef OPALX_EMITTED_FROM_FILE_H
#define OPALX_EMITTED_FROM_FILE_H

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include "Distribution.h"
#include "Ippl.h"
#include "OPALTypes.h"
#include "SamplingBase.hpp"

/**
 * @class EmittedFromFile
 * @brief Reads emitted distributions and emits particles at per-record birth times.
 *
 * Legacy old-OPAL rows are interpreted positionally as
 * `x px y py t pz [optional-bin]`; their time is negated and recentered using
 * the old-OPAL emission-window convention.
 *
 * Explicit spacetime rows require a named header containing
 * `x y z px py pz birth_time` (`t` is accepted as an alias for `birth_time`
 * when `z` is present). Columns are mapped by name. Positions are offsets from
 * the emission-source `R0`, momenta are normalized @f$p/(mc)@f$ offsets from
 * `P0`, and `birth_time` is a time offset in seconds from the source `T0`.
 */
class EmittedFromFile : public SamplingBase {
public:
    using size_type = ippl::detail::size_type;

    /**
     * @brief Constructor for EmittedFromFile.
     *
     * @param pc Shared pointer to ParticleContainer.
     * @param fc Shared pointer to FieldContainer.
     * @param opalDist Borrowed Distribution.
     */
    EmittedFromFile(
            std::shared_ptr<ParticleContainer_t> pc, std::shared_ptr<FieldContainer_t> fc,
            Distribution_t* opalDist);

    /**
     * @brief Convenience constructor that takes the filename directly.
     *
     * This is primarily intended for tests, where constructing a full OPALX
     * Distribution object is unnecessary.
     *
     * @param pc Shared pointer to ParticleContainer.
     * @param fc Shared pointer to FieldContainer.
     * @param filename Path to the emitted particle file.
     */
    EmittedFromFile(
            std::shared_ptr<ParticleContainer_t> pc, std::shared_ptr<FieldContainer_t> fc,
            const std::string& filename);

    /**
     * @brief Reads the selected file records and builds the birth-time inventory.
     *
     * @param numberOfParticles Number of particles requested by the caller.
     * @param nr Number of grid points per direction (not used here).
     */
    void generateParticles(size_t& numberOfParticles, Vector_t<double, 3> nr) override;

    /**
     * @brief Emits all file records whose birth times fall into the current step.
     *
     * @param t Start time of the current tracker step.
     * @param dt Time step size.
     */
    void emitParticles(double t, double dt) override;

    /**
     * @brief Returns whether all file records in the inventory have been emitted.
     *
     * @return True if the inventory exists and the next global index is past its end.
     */
    bool isEmissionDone(double /*t*/) const override {
        return inventoryBuilt_m && nextGlobalIndex_m >= records_m.size();
    }

    /// @copydoc SamplingBase::getEmittedFraction
    double getEmittedFraction() const override {
        if (!inventoryBuilt_m || records_m.empty()) {
            return 0.0;
        }
        return std::clamp(
                static_cast<double>(nextGlobalIndex_m) / static_cast<double>(records_m.size()), 0.0,
                1.0);
    }

    /**
     * @brief Reports whether an initial reference momentum is available.
     *
     * @return True after the file inventory has been built.
     */
    bool hasInitialReferenceMomentum() const override { return inventoryBuilt_m; }

    /**
     * @brief Returns the average initial reference momentum from selected records.
     */
    Vector_t<double, 3> getInitialReferenceMomentum() const override { return initialRefP_m; }

    /**
     * @brief Reports that this sampler provides an initial reference position.
     *
     * @return True after the file inventory has been built.
     */
    bool hasInitialReferencePosition() const override { return inventoryBuilt_m; }

    /**
     * @brief Returns the mean initial position used by the tracker.
     *
     * @return The mean selected file position plus the emission-source R0.
     */
    Vector_t<double, 3> getInitialReferencePosition() const override { return initialRefR_m; }

    /**
     * @brief Returns the global time shift needed to reach the earliest birth.
     *
     * Legacy files retain their centered old-OPAL shift. Explicit files shift
     * only when `T0 + min(birth_time)` is negative.
     *
     * @return Non-negative amount subtracted from the initial tracker time.
     */
    double getGlobalTimeShift() const override { return globalTimeShift_m; }

    /**
     * @brief Returns the preferred emission time step.
     *
     * @return emissionTime divided by the configured number of emission steps.
     */
    double getEmissionTimeStep() const override {
        return emissionSteps_m > 0 && emissionTime_m > 0.0
                       ? emissionTime_m / static_cast<double>(emissionSteps_m)
                       : 0.0;
    }

    /**
     * @brief Returns the total emission time spanned by selected file records.
     */
    double getEmissionTime() const { return emissionTime_m; }

    /**
     * @brief Generates the local particles assigned to this rank for one emission step.
     *
     * Particle coordinates and momentum offsets come from the file. Each particle
     * receives a per-particle fractional step based on its birth time and is
     * drifted for half of that effective step.
     *
     * @param nlocalBefore Local particle count before allocation.
     * @param globalBegin First global record index assigned to this rank.
     * @param nNew Number of local particles to generate.
     * @param tStart Start time of the current tracker step.
     * @param dt Time step size.
     */
    void generateLocalParticles(
            size_type nlocalBefore, size_t globalBegin, size_t nNew, double tStart, double dt);

private:
    /**
     * @brief Raw row parsed from a legacy or explicit emitted distribution.
     */
    struct RawRecord {
        double x        = 0.0;    ///< Horizontal position from the file.
        double px       = 0.0;    ///< Horizontal momentum offset from the file.
        double y        = 0.0;    ///< Vertical position from the file.
        double py       = 0.0;    ///< Vertical momentum offset from the file.
        double z        = 0.0;    ///< Explicit longitudinal position; zero for legacy rows.
        double fileTime = 0.0;    ///< Legacy time or explicit birth-time offset.
        double pz       = 0.0;    ///< Longitudinal momentum offset from the file.
        size_t bin      = 0;      ///< Optional old-OPAL emission bin number.
        bool hasBin     = false;  ///< True if the optional bin number was present.
    };

    /**
     * @brief Selected file row plus its converted tracker birth time.
     */
    struct Record {
        RawRecord raw;           ///< Original parsed file row.
        double birthTime = 0.0;  ///< Birth time in the OPALX tracker convention.
    };

    /**
     * @brief Resolves relative file paths against the active input file directory.
     */
    void resolveFilenameFromInput();

    /**
     * @brief Reads and parses a legacy or explicit emitted particle file.
     *
     * Accepted files may contain a leading particle count and header, or a
     * comment header. Without an explicit spacetime header, data rows use the
     * legacy positional format `x px y py t pz [optional-bin]`.
     *
     * @param filename Path to the file to read.
     */
    void readFile(const std::string& filename);

    /**
     * @brief Selects records, converts file times to birth times, and sorts them.
     *
     * @param requested Number of records requested; zero selects all records.
     */
    void buildInventory(size_t requested);

    /**
     * @brief Computes the local contiguous subrange of a global emission batch.
     *
     * @param totalToEmit Number of particles emitted globally in the current step.
     * @return Pair of local offset into the global batch and local particle count.
     */
    std::pair<size_t, size_t> computeLocalEmitRange(size_t totalToEmit) const;

    std::string filename_m;                     ///< File path to read emitted particle data from.
    std::vector<RawRecord> rawRecords_m;        ///< Parsed file rows before selection and sorting.
    std::vector<Record> records_m;              ///< Selected records sorted by tracker birth time.
    bool explicitBirthFormat_m        = false;  ///< True for named x/y/z/.../birth_time files.
    size_t nextGlobalIndex_m          = 0;      ///< First global record index not emitted yet.
    bool inventoryBuilt_m             = false;  ///< True once records_m is ready for emission.
    Vector_t<double, 3> initialRefR_m = 0.0;    ///< Average initial reference position.
    Vector_t<double, 3> initialRefP_m = 0.0;    ///< Average initial reference momentum.
    double emissionTime_m             = 0.0;    ///< Total emission time spanned by records_m.
    double globalTimeShift_m          = 0.0;    ///< Shift needed to include the earliest birth.
    size_t emissionSteps_m            = 100;    ///< Number of steps used to derive emission dt.
};

#endif  // OPALX_EMITTED_FROM_FILE_H
