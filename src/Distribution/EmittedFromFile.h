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

    EmittedFromFile(
            std::shared_ptr<ParticleContainer_t> pc, std::shared_ptr<FieldContainer_t> fc,
            Distribution_t* opalDist);

    EmittedFromFile(
            std::shared_ptr<ParticleContainer_t> pc, std::shared_ptr<FieldContainer_t> fc,
            const std::string& filename);

    void generateParticles(size_t& numberOfParticles, Vector_t<double, 3> nr) override;

    void emitParticles(double t, double dt) override;

    bool isEmissionDone(double /*t*/) const override {
        return inventoryBuilt_m && nextGlobalIndex_m >= records_m.size();
    }

    bool hasInitialReferenceMomentum() const override { return inventoryBuilt_m; }

    Vector_t<double, 3> getInitialReferenceMomentum() const override { return initialRefP_m; }

    double getGlobalTimeShift() const override {
        return std::max(0.0, 0.5 * emissionTime_m - t0_m);
    }

    double getEmissionTimeStep() const override {
        return emissionSteps_m > 0 && emissionTime_m > 0.0
                       ? emissionTime_m / static_cast<double>(emissionSteps_m)
                       : 0.0;
    }

    double getEmissionTime() const { return emissionTime_m; }

private:
    struct RawRecord {
        double x        = 0.0;
        double px       = 0.0;
        double y        = 0.0;
        double py       = 0.0;
        double fileTime = 0.0;
        double pz       = 0.0;
    };

    struct Record {
        RawRecord raw;
        double birthTime = 0.0;
    };

    void resolveFilenameFromInput();

    void readFile(const std::string& filename);

    void buildInventory(size_t requested);

    std::pair<size_t, size_t> computeLocalEmitRange(size_t totalToEmit) const;

    void generateLocalParticles(
            size_type nlocalBefore, size_t globalBegin, size_t nNew, double tStart, double dt);

    std::string filename_m;
    std::vector<RawRecord> rawRecords_m;
    std::vector<Record> records_m;
    size_t nextGlobalIndex_m  = 0;
    bool inventoryBuilt_m     = false;
    Vector_t<double, 3> initialRefP_m = 0.0;
    double emissionTime_m     = 0.0;
    size_t emissionSteps_m    = 100;
};

#endif  // OPALX_EMITTED_FROM_FILE_H
