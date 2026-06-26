/**
 * @file EmittedFromFile.h
 * @brief Defines an emitted file distribution with old-OPAL time-column semantics.
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
 * @brief Reads old-OPAL emitted distribution dumps and emits particles by file time.
 *
 * The data rows are interpreted positionally as
 * x px y py t pz [optional-bin].  The file time is an old-OPAL pre-emission
 * time, so a row with t = -6e-13 is emitted at t = +6e-13 when the source t0 is
 * zero.
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
     * @brief Returns the global time shift needed to center old-OPAL file times.
     *
     * @return Non-negative shift from source time to the midpoint of the emission window.
     */
    double getGlobalTimeShift() const override {
        return std::max(0.0, 0.5 * emissionTime_m - t0_m);
    }

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
     * @brief Raw row parsed from an old-OPAL emitted distribution dump.
     */
    struct RawRecord {
        double x        = 0.0;  ///< Horizontal position from the file.
        double px       = 0.0;  ///< Horizontal momentum offset from the file.
        double y        = 0.0;  ///< Vertical position from the file.
        double py       = 0.0;  ///< Vertical momentum offset from the file.
        double fileTime = 0.0;  ///< Old-OPAL pre-emission time column.
        double pz       = 0.0;  ///< Longitudinal momentum offset from the file.
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
     * @brief Reads and parses an old-OPAL emitted particle file.
     *
     * Accepted files may contain a leading particle count and header, or a
     * comment header. Data rows are interpreted positionally as
     * x px y py t pz [optional-bin].
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
    size_t nextGlobalIndex_m          = 0;      ///< First global record index not emitted yet.
    bool inventoryBuilt_m             = false;  ///< True once records_m is ready for emission.
    Vector_t<double, 3> initialRefP_m = 0.0;    ///< Average initial reference momentum.
    double emissionTime_m             = 0.0;    ///< Total emission time spanned by records_m.
    size_t emissionSteps_m            = 100;    ///< Number of steps used to derive emission dt.
};

#endif  // OPALX_EMITTED_FROM_FILE_H
