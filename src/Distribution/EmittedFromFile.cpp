#include "EmittedFromFile.h"

#include <mpi.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <optional>
#include <sstream>
#include <unordered_map>

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

    struct HeaderLayout {
        bool explicitBirth = false;
        std::unordered_map<std::string, size_t> column;
    };

    std::optional<HeaderLayout> parseHeader(const std::string& line) {
        std::istringstream stream(line);
        HeaderLayout layout;
        std::string token;
        size_t index = 0;
        while (stream >> token) {
            token = lowercase(token);
            if (!layout.column.emplace(token, index).second) {
                return std::nullopt;
            }
            ++index;
        }

        const auto has = [&layout](const std::string& name) {
            return layout.column.find(name) != layout.column.end();
        };
        const bool common = has("x") && has("y") && has("px") && has("py") && has("pz");
        if (common && has("z") && (has("birth_time") || has("t"))) {
            layout.explicitBirth = true;
            return layout;
        }
        if (common && !has("z") && has("t")) {
            return layout;
        }
        return std::nullopt;
    }

    double getConfiguredEmissionTime(const Distribution_t* opalDist) {
        if (!opalDist || !opalDist->emitting_m) {
            return 0.0;
        }

        const double sigmaTRise = opalDist->getSigmaTRise();
        const double sigmaTFall = opalDist->getSigmaTFall();
        const double pulseFwhm  = opalDist->getTPulseLengthFWHM();
        const double cutoffLong = opalDist->getCutoffR()[2];
        if (pulseFwhm <= 0.0 || cutoffLong <= 0.0 || (sigmaTRise <= 0.0 && sigmaTFall <= 0.0)) {
            return 0.0;
        }

        double flattopTime = pulseFwhm - std::sqrt(2.0 * std::log(2.0)) * (sigmaTRise + sigmaTFall);
        if (flattopTime < 0.0) {
            flattopTime = 0.0;
        }
        return cutoffLong * (sigmaTRise + sigmaTFall) + flattopTime;
    }

    bool parseNumericValues(const std::string& line, std::vector<double>& values) {
        std::istringstream stream(line);
        double value = 0.0;
        while (stream >> value) {
            values.push_back(value);
        }
        return stream.eof();
    }
}  // namespace

EmittedFromFile::EmittedFromFile(
        std::shared_ptr<ParticleContainer_t> pc, std::shared_ptr<FieldContainer_t> fc,
        Distribution_t* opalDist)
    : SamplingBase(pc, fc, opalDist) {
    filename_m      = opalDist->getFilename();
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
    explicitBirthFormat_m   = false;
    bool sawFirstPayload    = false;
    bool expectHeaderAfterN = false;
    bool hasDeclaredCount   = false;
    size_t declaredCount    = 0;
    size_t lineNumber       = 0;
    std::optional<HeaderLayout> header;

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
                declaredCount      = static_cast<size_t>(count);
                hasDeclaredCount   = true;
                expectHeaderAfterN = true;
                sawFirstPayload    = true;
                continue;
            }
            sawFirstPayload = true;
        }

        if (expectHeaderAfterN) {
            header = parseHeader(stripped);
            if (!header.has_value()) {
                throw OpalException(
                        "EmittedFromFile::readFile",
                        "Expected a legacy 'x px y py t pz' or explicit "
                        "'x y z px py pz birth_time' header after the count line in '"
                                + filename + "'.");
            }
            explicitBirthFormat_m = header->explicitBirth;
            expectHeaderAfterN    = false;
            continue;
        }

        if (const auto parsedHeader = parseHeader(stripped); parsedHeader.has_value()) {
            if (!rawRecords_m.empty()) {
                throw OpalException(
                        "EmittedFromFile::readFile", "A header appears after particle data on line "
                                                             + std::to_string(lineNumber) + " in '"
                                                             + filename + "'.");
            }
            header                = parsedHeader;
            explicitBirthFormat_m = header->explicitBirth;
            continue;
        }

        std::vector<double> values;
        if (!parseNumericValues(stripped, values)) {
            throw OpalException(
                    "EmittedFromFile::readFile", "Line " + std::to_string(lineNumber) + " in '"
                                                         + filename + "' is not fully numeric.");
        }
        RawRecord record;
        if (explicitBirthFormat_m) {
            const auto valueAt = [&](const std::string& name) -> double {
                const size_t column = header->column.at(name);
                if (column >= values.size()) {
                    throw OpalException(
                            "EmittedFromFile::readFile",
                            "Line " + std::to_string(lineNumber) + " in '" + filename
                                    + "' has fewer columns than its explicit header.");
                }
                return values[column];
            };
            record.x        = valueAt("x");
            record.y        = valueAt("y");
            record.z        = valueAt("z");
            record.px       = valueAt("px");
            record.py       = valueAt("py");
            record.pz       = valueAt("pz");
            record.fileTime = valueAt(header->column.count("birth_time") ? "birth_time" : "t");
        } else {
            if (values.size() < 6) {
                throw OpalException(
                        "EmittedFromFile::readFile",
                        "Line " + std::to_string(lineNumber) + " in '" + filename
                                + "' has fewer than six numeric columns.");
            }
            record.x        = values[0];
            record.px       = values[1];
            record.y        = values[2];
            record.py       = values[3];
            record.fileTime = values[4];
            record.pz       = values[5];
            if (values.size() >= 7) {
                const double rawBin     = values[6];
                const double roundedBin = std::round(rawBin);
                if (roundedBin <= 0.0 || std::fabs(rawBin - roundedBin) > 1.0e-9) {
                    throw OpalException(
                            "EmittedFromFile::readFile",
                            "Line " + std::to_string(lineNumber) + " in '" + filename
                                    + "' optional bin column must be a positive integer.");
                }
                record.bin    = static_cast<size_t>(roundedBin);
                record.hasBin = true;
            }
        }
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
    initialRefR_m     = R0_m;
    initialRefP_m     = P0_m;
    emissionTime_m    = 0.0;
    globalTimeShift_m = 0.0;

    if (selected == 0) {
        if (opalDist_m) {
            opalDist_m->setTEmission(emissionTime_m);
        }
        inventoryBuilt_m = true;
        return;
    }

    double minPulseTime =
            explicitBirthFormat_m ? rawRecords_m[0].fileTime : -rawRecords_m[0].fileTime;
    double maxPulseTime     = minPulseTime;
    bool allRecordsHaveBins = rawRecords_m[0].hasBin;
    size_t maxBin           = rawRecords_m[0].bin;
    for (size_t i = 1; i < selected; ++i) {
        const double pulseTime =
                explicitBirthFormat_m ? rawRecords_m[i].fileTime : -rawRecords_m[i].fileTime;
        minPulseTime       = std::min(minPulseTime, pulseTime);
        maxPulseTime       = std::max(maxPulseTime, pulseTime);
        allRecordsHaveBins = allRecordsHaveBins && rawRecords_m[i].hasBin;
        maxBin             = std::max(maxBin, rawRecords_m[i].bin);
    }

    double pulseCenter = 0.5 * (minPulseTime + maxPulseTime);
    if (explicitBirthFormat_m) {
        emissionTime_m = std::max(0.0, maxPulseTime - minPulseTime);
    } else if (allRecordsHaveBins && maxBin > 0) {
        double lowerEmissionTime = 0.0;
        double upperEmissionTime = std::numeric_limits<double>::infinity();
        for (size_t i = 0; i < selected; ++i) {
            const RawRecord& raw   = rawRecords_m[i];
            const double pulseTime = -raw.fileTime;
            lowerEmissionTime      = std::max(
                    lowerEmissionTime,
                    pulseTime * static_cast<double>(maxBin) / static_cast<double>(raw.bin));
            if (raw.bin > 1) {
                upperEmissionTime = std::min(
                        upperEmissionTime,
                        pulseTime * static_cast<double>(maxBin) / static_cast<double>(raw.bin - 1));
            }
        }

        emissionTime_m = std::isfinite(upperEmissionTime) && upperEmissionTime >= lowerEmissionTime
                                 ? 0.5 * (lowerEmissionTime + upperEmissionTime)
                                 : lowerEmissionTime;
        pulseCenter    = 0.5 * emissionTime_m;
    } else {
        emissionTime_m = std::max(0.0, maxPulseTime - minPulseTime);
    }

    const double configuredEmissionTime = getConfiguredEmissionTime(opalDist_m);
    if (!explicitBirthFormat_m && configuredEmissionTime > 0.0) {
        const double tolerance = std::max(1.0e-18, configuredEmissionTime * 1.0e-12);
        if (configuredEmissionTime + tolerance < maxPulseTime) {
            throw OpalException(
                    "EmittedFromFile::buildInventory",
                    "Configured emission window is shorter than the latest emitted file time.");
        }
        emissionTime_m = configuredEmissionTime;
        pulseCenter    = 0.5 * emissionTime_m;
    }
    if (opalDist_m) {
        opalDist_m->setTEmission(emissionTime_m);
    }
    globalTimeShift_m = explicitBirthFormat_m ? std::max(0.0, -(t0_m + minPulseTime))
                                              : std::max(0.0, 0.5 * emissionTime_m - t0_m);

    double pxSum = 0.0;
    double pySum = 0.0;
    double pzSum = 0.0;
    Vector_t<double, 3> rSum(0.0);
    for (size_t i = 0; i < selected; ++i) {
        Record record;
        record.raw       = rawRecords_m[i];
        record.birthTime = explicitBirthFormat_m ? t0_m + record.raw.fileTime
                                                 : t0_m - record.raw.fileTime - pulseCenter;
        records_m.push_back(record);
        rSum[0] += record.raw.x + R0_m[0];
        rSum[1] += record.raw.y + R0_m[1];
        rSum[2] += record.raw.z + R0_m[2];
        pxSum += record.raw.px + P0_m[0];
        pySum += record.raw.py + P0_m[1];
        pzSum += record.raw.pz + P0_m[2];
    }

    std::stable_sort(records_m.begin(), records_m.end(), [](const Record& lhs, const Record& rhs) {
        return lhs.birthTime < rhs.birthTime;
    });

    const double invSelected = 1.0 / static_cast<double>(selected);
    if (explicitBirthFormat_m) {
        initialRefR_m[0] = rSum[0] * invSelected;
        initialRefR_m[1] = rSum[1] * invSelected;
        initialRefR_m[2] = rSum[2] * invSelected;
    }
    initialRefP_m[0] = pxSum * invSelected;
    initialRefP_m[1] = pySum * invSelected;
    initialRefP_m[2] = pzSum * invSelected;
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
        const size_t ideal             = base + (static_cast<size_t>(r) < rem ? 1 : 0);
        const size_t cap               = static_cast<size_t>(spaceLeftAll[static_cast<size_t>(r)]);
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
            [](double value, const Record& record) {
                return value < record.birthTime;
            });
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
    HostView hostData("EmittedFromFile_hostData", nNew, 7);
    const double tEnd             = tStart + dt;
    const double overdueTolerance = std::max(1.0e-18, std::abs(dt) * 1.0e-12);

    for (size_t i = 0; i < nNew; ++i) {
        const Record& record = records_m[globalBegin + i];
        if (record.birthTime < tStart - overdueTolerance) {
            throw OpalException(
                    "EmittedFromFile::generateLocalParticles",
                    "EMITTEDFROMFILE would need to pre-age a particle born before "
                    "the current step.");
        }
        const double effectiveBirthTime = record.birthTime < tStart ? tStart : record.birthTime;
        const double stepDt             = tEnd - effectiveBirthTime;

        hostData(i, 0) = record.raw.x;
        hostData(i, 1) = record.raw.y;
        hostData(i, 2) = record.raw.z;
        hostData(i, 3) = record.raw.px;
        hostData(i, 4) = record.raw.py;
        hostData(i, 5) = record.raw.pz;
        hostData(i, 6) = stepDt;
    }

    auto deviceData =
            Kokkos::create_mirror(Kokkos::DefaultExecutionSpace::memory_space(), hostData);
    Kokkos::deep_copy(deviceData, hostData);

    auto Rview                   = pc_m->R.getView();
    auto Pview                   = pc_m->P.getView();
    auto dtView                  = pc_m->dt.getView();
    const Vector_t<double, 3> R0 = R0_m;
    const Vector_t<double, 3> P0 = P0_m;
    const size_type offset       = nlocalBefore;
    const double c               = Physics::c;

    Kokkos::parallel_for(
            "EmittedFromFile_generateLocalParticles", nNew, KOKKOS_LAMBDA(const size_t i) {
                const size_t j = offset + i;
                Vector_t<double, 3> p(
                        deviceData(i, 3) + P0[0], deviceData(i, 4) + P0[1],
                        deviceData(i, 5) + P0[2]);
                const double stepDt = deviceData(i, 6);

                Rview(j)[0] = deviceData(i, 0) + R0[0];
                Rview(j)[1] = deviceData(i, 1) + R0[1];
                Rview(j)[2] = deviceData(i, 2) + R0[2];
                Pview(j)    = p;
                dtView(j)   = stepDt;

                const double gamma = Kokkos::sqrt(1.0 + p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
                const double drift = 0.5 * c * stepDt / gamma;
                Rview(j)[0] += p[0] * drift;
                Rview(j)[1] += p[1] * drift;
                Rview(j)[2] += p[2] * drift;
            });
    Kokkos::fence();

    pc_m->markMomentsDirty();
}
