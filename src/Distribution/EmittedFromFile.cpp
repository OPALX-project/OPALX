#include "EmittedFromFile.h"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>

#include <Kokkos_Core.hpp>

#include "AbstractObjects/OpalData.h"
#include "Physics/Physics.h"
#include "Utilities/OpalException.h"

namespace {
    std::string trim(const std::string& input) {
        const auto first = input.find_first_not_of(" \t\r\n");
        if (first == std::string::npos) {
            return "";
        }
        const auto last = input.find_last_not_of(" \t\r\n");
        return input.substr(first, last - first + 1);
    }

    std::string lowercase(std::string value) {
        std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
            return static_cast<char>(std::tolower(c));
        });
        return value;
    }

    bool parseSingleInteger(const std::string& line, long long& value) {
        std::istringstream stream(line);
        stream >> value;
        return stream && stream.eof();
    }

    bool looksLikeHeader(const std::string& line) {
        const std::string lower = lowercase(line);
        return lower.find("x") != std::string::npos && lower.find("px") != std::string::npos
               && lower.find("y") != std::string::npos
               && lower.find("py") != std::string::npos
               && lower.find("t") != std::string::npos
               && lower.find("pz") != std::string::npos;
    }

    bool parseRecordValues(const std::string& line, std::array<double, 6>& record) {
        std::istringstream stream(line);
        std::vector<double> values;
        double value = 0.0;
        while (stream >> value) {
            values.push_back(value);
        }
        if (values.size() < 6) {
            return false;
        }

        for (size_t i = 0; i < record.size(); ++i) {
            record[i] = values[i];
        }
        return true;
    }
}  // namespace

EmittedFromFile::EmittedFromFile(
        std::shared_ptr<ParticleContainer_t> pc, std::shared_ptr<FieldContainer_t> fc,
        Distribution_t* opalDist)
    : SamplingBase(pc, fc, opalDist) {
    filename_m = opalDist->getFilename();
    emissionSteps_m = opalDist->getEmissionSteps();
    if (filename_m.empty()) {
        throw OpalException(
                "EmittedFromFile::EmittedFromFile",
                "FNAME attribute must be set for EMITTEDFROMFILE distribution type.");
    }

    resolveFilenameFromInput();
    readFile(filename_m);
}

EmittedFromFile::EmittedFromFile(
        std::shared_ptr<ParticleContainer_t> pc, std::shared_ptr<FieldContainer_t> fc,
        const std::string& filename)
    : SamplingBase(pc, fc), filename_m(filename) {
    if (filename_m.empty()) {
        throw OpalException("EmittedFromFile::EmittedFromFile", "Filename must not be empty.");
    }

    readFile(filename_m);
}

void EmittedFromFile::resolveFilenameFromInput() {
    namespace fs = std::filesystem;
    if (fs::exists(filename_m)) {
        return;
    }

    const std::string inputFile = OpalData::getInstance()->getInputFn();
    if (!inputFile.empty()) {
        const fs::path filePath = fs::path(inputFile).parent_path() / filename_m;
        if (fs::exists(filePath)) {
            filename_m = filePath.string();
        }
    }
}

void EmittedFromFile::readFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw OpalException("EmittedFromFile::readFile", "Couldn't open file '" + filename + "'.");
    }

    rawRecords_m.clear();
    bool sawFirstPayload      = false;
    bool expectHeaderAfterN   = false;
    bool hasDeclaredCount     = false;
    size_t declaredCount      = 0;
    size_t lineNumber         = 0;

    std::string line;
    while (std::getline(file, line)) {
        ++lineNumber;
        std::string stripped = trim(line);
        if (stripped.empty()) {
            continue;
        }

        if (stripped[0] == '#') {
            // Old OPAL emitted dumps use a comment header.  The data order is
            // positional, so the header is accepted but not needed for mapping.
            continue;
        }

        if (!sawFirstPayload) {
            long long count = 0;
            if (parseSingleInteger(stripped, count)) {
                if (count < 0) {
                    throw OpalException(
                            "EmittedFromFile::readFile",
                            "Declared particle count in '" + filename + "' must be non-negative.");
                }
                declaredCount    = static_cast<size_t>(count);
                hasDeclaredCount = true;
                expectHeaderAfterN = true;
                sawFirstPayload  = true;
                continue;
            }
            sawFirstPayload = true;
        }

        if (expectHeaderAfterN) {
            if (!looksLikeHeader(stripped)) {
                throw OpalException(
                        "EmittedFromFile::readFile",
                        "Expected a header containing x px y py t pz after the count line in '"
                                + filename + "'.");
            }
            expectHeaderAfterN = false;
            continue;
        }

        if (looksLikeHeader(stripped)) {
            continue;
        }

        std::array<double, 6> values;
        if (!parseRecordValues(stripped, values)) {
            throw OpalException(
                    "EmittedFromFile::readFile",
                    "Line " + std::to_string(lineNumber) + " in '" + filename
                            + "' has fewer than six numeric columns.");
        }
        RawRecord record;
        record.x        = values[0];
        record.px       = values[1];
        record.y        = values[2];
        record.py       = values[3];
        record.fileTime = values[4];
        record.pz       = values[5];
        rawRecords_m.push_back(record);
    }

    if (expectHeaderAfterN) {
        throw OpalException(
                "EmittedFromFile::readFile",
                "Missing header line after declared particle count in '" + filename + "'.");
    }

    if (rawRecords_m.empty()) {
        throw OpalException(
                "EmittedFromFile::readFile",
                "No valid emitted particle data found in '" + filename + "'.");
    }

    if (hasDeclaredCount && rawRecords_m.size() != declaredCount) {
        throw OpalException(
                "EmittedFromFile::readFile",
                "Number of data lines (" + std::to_string(rawRecords_m.size())
                        + ") does not match declared count (" + std::to_string(declaredCount)
                        + ") in '" + filename + "'.");
    }
}

void EmittedFromFile::generateParticles(size_t& numberOfParticles, Vector_t<double, 3> /*nr*/) {
    if (emissionModel_m != "NONE") {
        throw OpalException(
                "EmittedFromFile::generateParticles",
                "EMISSIONMODEL '" + emissionModel_m
                        + "' is not supported for EMITTEDFROMFILE distributions");
    }

    size_t requested = numberOfParticles;
    if (opalDist_m && opalDist_m->getNumParticles() > 0) {
        requested = opalDist_m->getNumParticles();
    }

    buildInventory(requested);
    numberOfParticles = records_m.size();
}

void EmittedFromFile::buildInventory(size_t requested) {
    const size_t selected =
            requested > 0 ? std::min(requested, rawRecords_m.size()) : rawRecords_m.size();

    records_m.clear();
    records_m.reserve(selected);
    nextGlobalIndex_m = 0;
    inventoryBuilt_m  = false;
    initialRefP_m     = P0_m;
    emissionTime_m    = 0.0;

    if (selected == 0) {
        if (opalDist_m) {
            opalDist_m->setTEmission(emissionTime_m);
        }
        inventoryBuilt_m = true;
        return;
    }

    double minPulseTime = -rawRecords_m[0].fileTime;
    double maxPulseTime = minPulseTime;
    for (size_t i = 1; i < selected; ++i) {
        const double pulseTime = -rawRecords_m[i].fileTime;
        minPulseTime = std::min(minPulseTime, pulseTime);
        maxPulseTime = std::max(maxPulseTime, pulseTime);
    }
    emissionTime_m = std::max(0.0, maxPulseTime - minPulseTime);
    const double pulseCenter = 0.5 * (minPulseTime + maxPulseTime);
    if (opalDist_m) {
        opalDist_m->setTEmission(emissionTime_m);
    }

    double pxSum = 0.0;
    double pySum = 0.0;
    double pzSum = 0.0;
    for (size_t i = 0; i < selected; ++i) {
        Record record;
        record.raw       = rawRecords_m[i];
        record.birthTime = t0_m - record.raw.fileTime - pulseCenter;
        records_m.push_back(record);
        pxSum += record.raw.px + P0_m[0];
        pySum += record.raw.py + P0_m[1];
        pzSum += record.raw.pz + P0_m[2];
    }

    std::stable_sort(records_m.begin(), records_m.end(), [](const Record& lhs, const Record& rhs) {
        return lhs.birthTime < rhs.birthTime;
    });

    const double invSelected = 1.0 / static_cast<double>(selected);
    initialRefP_m[0]         = pxSum * invSelected;
    initialRefP_m[1]         = pySum * invSelected;
    initialRefP_m[2]         = pzSum * invSelected;
    inventoryBuilt_m = true;
}

std::pair<size_t, size_t> EmittedFromFile::computeLocalEmitRange(size_t totalToEmit) const {
    if (!pc_m || totalToEmit == 0) {
        return {0, 0};
    }

    const int nranks = ippl::Comm->size();
    const int rank   = ippl::Comm->rank();
    if (nranks <= 0) {
        return {0, totalToEmit};
    }

    const size_t capacity  = pc_m->R.size();
    const size_t localNum  = pc_m->getLocalNum();
    const size_t spaceLeft = (localNum <= capacity) ? (capacity - localNum) : size_t(0);

    std::vector<unsigned long> spaceLeftAll(static_cast<size_t>(nranks), 0);
    unsigned long mySpace = static_cast<unsigned long>(spaceLeft);
    MPI_Allgather(
            &mySpace, 1, MPI_UNSIGNED_LONG, spaceLeftAll.data(), 1, MPI_UNSIGNED_LONG,
            ippl::Comm->getCommunicator());

    const size_t nranksU = static_cast<size_t>(nranks);
    const size_t base    = totalToEmit / nranksU;
    const size_t rem     = totalToEmit % nranksU;

    std::vector<size_t> nlocal(nranksU, 0);
    size_t assigned = 0;
    for (int r = 0; r < nranks; ++r) {
        const size_t ideal = base + (static_cast<size_t>(r) < rem ? 1 : 0);
        const size_t cap   = static_cast<size_t>(spaceLeftAll[static_cast<size_t>(r)]);
        nlocal[static_cast<size_t>(r)] = std::min(ideal, cap);
        assigned += nlocal[static_cast<size_t>(r)];
    }

    size_t deficit = totalToEmit - assigned;
    while (deficit > 0) {
        int chosen = -1;
        for (int r = 0; r < nranks; ++r) {
            const size_t cap = static_cast<size_t>(spaceLeftAll[static_cast<size_t>(r)]);
            if (nlocal[static_cast<size_t>(r)] < cap) {
                chosen = r;
                break;
            }
        }
        if (chosen < 0) {
            throw OpalException(
                    "EmittedFromFile::computeLocalEmitRange",
                    "Particle container capacity is insufficient for EMITTEDFROMFILE emission. "
                    "Increase BEAM::NALLOC or reduce NPARTDIST.");
        }
        nlocal[static_cast<size_t>(chosen)] += 1;
        --deficit;
    }

    size_t offset = 0;
    for (int r = 0; r < rank; ++r) {
        offset += nlocal[static_cast<size_t>(r)];
    }
    return {offset, nlocal[static_cast<size_t>(rank)]};
}

void EmittedFromFile::emitParticles(double t, double dt) {
    if (!inventoryBuilt_m || records_m.empty()) {
        return;
    }
    if (nextGlobalIndex_m >= records_m.size()) {
        return;
    }

    const double tEnd = t + dt;
    auto emitEndIt    = std::upper_bound(
            records_m.begin() + nextGlobalIndex_m, records_m.end(), tEnd,
            [](double value, const Record& record) { return value < record.birthTime; });
    const size_t globalEnd = static_cast<size_t>(emitEndIt - records_m.begin());
    const size_t totalNew  = globalEnd - nextGlobalIndex_m;
    if (totalNew == 0) {
        return;
    }

    const double overdueTolerance = std::max(1.0e-18, std::abs(dt) * 1.0e-12);
    if (records_m[nextGlobalIndex_m].birthTime < t - overdueTolerance) {
        throw OpalException(
                "EmittedFromFile::emitParticles",
                "EMITTEDFROMFILE found an overdue birth time. This means the tracker "
                "started too late or skipped a timestep containing file-emitted particles.");
    }

    const auto [localOffset, nNew] = computeLocalEmitRange(totalNew);
    const size_type nlocalBefore   = pc_m->getLocalNum();
    pc_m->createParticles(static_cast<size_type>(nNew));

    generateLocalParticles(nlocalBefore, nextGlobalIndex_m + localOffset, nNew, t, dt);
    nextGlobalIndex_m = globalEnd;
}

void EmittedFromFile::generateLocalParticles(
        size_type nlocalBefore, size_t globalBegin, size_t nNew, double tStart, double dt) {
    if (nNew == 0) {
        return;
    }

    using HostView = Kokkos::View<double**, Kokkos::HostSpace>;
    HostView hostData("EmittedFromFile_hostData", nNew, 6);
    const double tEnd              = tStart + dt;
    const double overdueTolerance  = std::max(1.0e-18, std::abs(dt) * 1.0e-12);

    for (size_t i = 0; i < nNew; ++i) {
        const Record& record = records_m[globalBegin + i];
        if (record.birthTime < tStart - overdueTolerance) {
            throw OpalException(
                    "EmittedFromFile::generateLocalParticles",
                    "EMITTEDFROMFILE would need to pre-age a particle born before "
                    "the current step.");
        }
        const double effectiveBirthTime =
                record.birthTime < tStart ? tStart : record.birthTime;
        const double stepDt = tEnd - effectiveBirthTime;

        hostData(i, 0) = record.raw.x;
        hostData(i, 1) = record.raw.y;
        hostData(i, 2) = record.raw.px;
        hostData(i, 3) = record.raw.py;
        hostData(i, 4) = record.raw.pz;
        hostData(i, 5) = stepDt;
    }

    auto deviceData =
            Kokkos::create_mirror(Kokkos::DefaultExecutionSpace::memory_space(), hostData);
    Kokkos::deep_copy(deviceData, hostData);

    auto Rview      = pc_m->R.getView();
    auto Pview      = pc_m->P.getView();
    auto dtView     = pc_m->dt.getView();
    const Vector_t<double, 3> R0 = R0_m;
    const Vector_t<double, 3> P0 = P0_m;
    const size_type offset       = nlocalBefore;
    const double c               = Physics::c;

    Kokkos::parallel_for(
            "EmittedFromFile_generateLocalParticles", nNew, KOKKOS_LAMBDA(const size_t i) {
                const size_t j = offset + i;
                Vector_t<double, 3> p(
                        deviceData(i, 2) + P0[0], deviceData(i, 3) + P0[1],
                        deviceData(i, 4) + P0[2]);
                const double stepDt = deviceData(i, 5);

                Rview(j)[0] = deviceData(i, 0) + R0[0];
                Rview(j)[1] = deviceData(i, 1) + R0[1];
                Rview(j)[2] = R0[2];
                Pview(j)    = p;
                dtView(j)   = stepDt;

                const double gamma =
                        Kokkos::sqrt(1.0 + p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
                const double drift = 0.5 * c * stepDt / gamma;
                Rview(j)[0] += p[0] * drift;
                Rview(j)[1] += p[1] * drift;
                Rview(j)[2] += p[2] * drift;
            });
    Kokkos::fence();

    pc_m->markMomentsDirty();
}
