/**
 * @file CheckpointFile.cpp
 * @brief Atomic, MPI-rank-count-independent OPALX checkpoint I/O.
 */

#include "Structure/CheckpointFile.h"

// Define real MPI types before H5hut to avoid its serial MPI stubs.
#include <mpi.h>

extern "C" {
#include "H5hut.h"
}

#include "AbstractObjects/OpalData.h"
#include "Algorithms/CoordinateSystemTrafo.h"
#include "Algorithms/Quaternion.hpp"
#include "PartBunch/PartBunch.h"
#include "Utilities/OpalException.h"

#include "Utility/IpplInfo.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <string>
#include <system_error>
#include <vector>

namespace {
    constexpr h5_int64_t checkpointFormatVersion = 1;

    void requireSuccess(h5_int64_t rc, const std::string& operation) {
        if (rc != H5_SUCCESS) {
            throw OpalException("CheckpointFile", operation + " failed");
        }
    }

    h5_file_t openFile(const std::string& fileName, h5_int32_t flags) {
        h5_prop_t props = H5CreateFileProp();
        if (props == H5_ERR) {
            throw OpalException("CheckpointFile", "could not create HDF5 file properties");
        }

        MPI_Comm comm         = ippl::Comm->getCommunicator();
        const h5_err_t propRc = H5SetPropFileMPIOCollective(props, &comm);
        if (propRc == H5_ERR) {
            H5CloseProp(props);
            throw OpalException("CheckpointFile", "could not enable collective HDF5 I/O");
        }

        h5_file_t file = H5OpenFile(fileName.c_str(), flags, props);
        H5CloseProp(props);
        if (file == static_cast<h5_file_t>(H5_ERR)) {
            throw OpalException("CheckpointFile", "could not open '" + fileName + "'");
        }
        return file;
    }

    void closeFile(h5_file_t& file) {
        if (file != 0) {
            ippl::Comm->barrier();
            requireSuccess(H5CloseFile(file), "closing checkpoint file");
            file = 0;
        }
    }

    void requireFileAttribute(h5_file_t file, const char* name) {
        if (H5HasFileAttrib(file, name) <= 0) {
            throw OpalException(
                    "CheckpointFile::inspect",
                    "required file attribute '" + std::string(name) + "' is missing");
        }
    }

    void requireStepAttribute(h5_file_t file, const char* name) {
        if (H5HasStepAttrib(file, name) <= 0) {
            throw OpalException(
                    "CheckpointFile::read",
                    "required step attribute '" + std::string(name) + "' is missing");
        }
    }

    void requireDataset(h5_file_t file, const char* name) {
        if (H5PartHasDataset(file, name) <= 0) {
            throw OpalException(
                    "CheckpointFile::read",
                    "required particle dataset '" + std::string(name) + "' is missing");
        }
    }

    CheckpointFile::Metadata readMetadata(h5_file_t file) {
        requireFileAttribute(file, "OPALXCheckpointFormat");
        requireFileAttribute(file, "GlobalTrackStep");
        requireFileAttribute(file, "TIME");
        requireFileAttribute(file, "DT");
        requireFileAttribute(file, "NumContainers");
        requireFileAttribute(file, "StepSizeSegment");
        requireFileAttribute(file, "StepsCompletedInSegment");

        h5_int64_t format       = 0;
        h5_int64_t step         = 0;
        h5_int64_t count        = 0;
        h5_int64_t segment      = 0;
        h5_int64_t segmentSteps = 0;
        h5_float64_t time       = 0.0;
        h5_float64_t dt         = 0.0;
        requireSuccess(
                H5ReadFileAttribInt64(file, "OPALXCheckpointFormat", &format),
                "reading checkpoint format version");
        requireSuccess(
                H5ReadFileAttribInt64(file, "GlobalTrackStep", &step),
                "reading global tracking step");
        requireSuccess(H5ReadFileAttribFloat64(file, "TIME", &time), "reading checkpoint time");
        requireSuccess(H5ReadFileAttribFloat64(file, "DT", &dt), "reading checkpoint dt");
        requireSuccess(
                H5ReadFileAttribInt64(file, "NumContainers", &count),
                "reading particle-container count");
        requireSuccess(
                H5ReadFileAttribInt64(file, "StepSizeSegment", &segment),
                "reading step-size segment");
        requireSuccess(
                H5ReadFileAttribInt64(file, "StepsCompletedInSegment", &segmentSteps),
                "reading completed steps in segment");

        if (format != checkpointFormatVersion) {
            throw OpalException(
                    "CheckpointFile::inspect", "unsupported checkpoint format version "
                                                       + std::to_string(format)
                                                       + "; this OPALX build supports version "
                                                       + std::to_string(checkpointFormatVersion));
        }
        if (step < 0 || count <= 0 || dt == 0.0 || segment < 0 || segmentSteps < 0) {
            throw OpalException(
                    "CheckpointFile::inspect", "checkpoint metadata contains invalid values");
        }
        if (H5GetNumSteps(file) != count) {
            throw OpalException(
                    "CheckpointFile::inspect",
                    "checkpoint is incomplete: expected " + std::to_string(count)
                            + " container records, found " + std::to_string(H5GetNumSteps(file)));
        }

        return CheckpointFile::Metadata{
                static_cast<long long>(step),      static_cast<double>(time),
                static_cast<double>(dt),           static_cast<std::size_t>(count),
                static_cast<std::size_t>(segment), static_cast<unsigned long long>(segmentSteps)};
    }

    std::string readStringFileAttribute(h5_file_t file, const std::string& name) {
        h5_int64_t type = 0;
        h5_size_t size  = 0;
        requireSuccess(
                H5GetFileAttribInfoByName(file, name.c_str(), &type, &size),
                "querying file attribute '" + name + "'");
        std::vector<char> buffer(static_cast<std::size_t>(size) + 1, '\0');
        requireSuccess(
                H5ReadFileAttribString(file, name.c_str(), buffer.data()),
                "reading file attribute '" + name + "'");
        return std::string(buffer.data());
    }

    void restoreCavityPhases(h5_file_t file) {
        if (H5HasFileAttrib(file, "NumCavityPhases") <= 0) {
            return;
        }
        h5_int64_t count = 0;
        requireSuccess(
                H5ReadFileAttribInt64(file, "NumCavityPhases", &count),
                "reading cavity phase count");
        if (count < 0) {
            throw OpalException("CheckpointFile::read", "invalid cavity phase count");
        }
        for (h5_int64_t i = 0; i < count; ++i) {
            const std::string prefix = "CavityPhase" + std::to_string(i);
            requireFileAttribute(file, (prefix + "Name").c_str());
            requireFileAttribute(file, (prefix + "Value").c_str());
            h5_float64_t phase = 0.0;
            requireSuccess(
                    H5ReadFileAttribFloat64(file, (prefix + "Value").c_str(), &phase),
                    "reading cavity phase value");
            OpalData::getInstance()->setMaxPhase(
                    readStringFileAttribute(file, prefix + "Name"), phase);
        }
    }

    template <typename View>
    auto hostCopy(const View& view) {
        return Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), view);
    }

    bool compatibleScalar(double configured, double saved) {
        const double scale = std::max(
                {std::abs(configured), std::abs(saved), std::numeric_limits<double>::min()});
        return std::abs(configured - saved)
               <= 64.0 * std::numeric_limits<double>::epsilon() * scale;
    }

    template <typename Getter>
    void writeFloatDataset(
            h5_file_t file, const char* name, std::size_t localCount, Getter&& getter) {
        std::vector<h5_float64_t> values(std::max<std::size_t>(localCount, 1));
        for (std::size_t i = 0; i < localCount; ++i) {
            values[i] = static_cast<h5_float64_t>(getter(i));
        }
        requireSuccess(
                H5PartWriteDataFloat64(file, name, values.data()),
                "writing particle dataset '" + std::string(name) + "'");
    }

    template <typename Getter>
    void writeInt64Dataset(
            h5_file_t file, const char* name, std::size_t localCount, Getter&& getter) {
        std::vector<h5_int64_t> values(std::max<std::size_t>(localCount, 1));
        for (std::size_t i = 0; i < localCount; ++i) {
            values[i] = static_cast<h5_int64_t>(getter(i));
        }
        requireSuccess(
                H5PartWriteDataInt64(file, name, values.data()),
                "writing particle dataset '" + std::string(name) + "'");
    }

    template <typename Getter>
    void writeInt32Dataset(
            h5_file_t file, const char* name, std::size_t localCount, Getter&& getter) {
        std::vector<h5_int32_t> values(std::max<std::size_t>(localCount, 1));
        for (std::size_t i = 0; i < localCount; ++i) {
            values[i] = static_cast<h5_int32_t>(getter(i));
        }
        requireSuccess(
                H5PartWriteDataInt32(file, name, values.data()),
                "writing particle dataset '" + std::string(name) + "'");
    }

    std::vector<h5_float64_t> readFloatDataset(
            h5_file_t file, const char* name, std::size_t localCount) {
        requireDataset(file, name);
        std::vector<h5_float64_t> values(std::max<std::size_t>(localCount, 1));
        requireSuccess(
                H5PartReadDataFloat64(file, name, values.data()),
                "reading particle dataset '" + std::string(name) + "'");
        return values;
    }

    std::vector<h5_int64_t> readInt64Dataset(
            h5_file_t file, const char* name, std::size_t localCount) {
        requireDataset(file, name);
        std::vector<h5_int64_t> values(std::max<std::size_t>(localCount, 1));
        requireSuccess(
                H5PartReadDataInt64(file, name, values.data()),
                "reading particle dataset '" + std::string(name) + "'");
        return values;
    }

    std::vector<h5_int32_t> readInt32Dataset(
            h5_file_t file, const char* name, std::size_t localCount) {
        requireDataset(file, name);
        std::vector<h5_int32_t> values(std::max<std::size_t>(localCount, 1));
        requireSuccess(
                H5PartReadDataInt32(file, name, values.data()),
                "reading particle dataset '" + std::string(name) + "'");
        return values;
    }
}  // namespace

std::string CheckpointFile::defaultPath(const std::string& inputBasename) {
    return inputBasename + "_checkpoint.h5";
}

void CheckpointFile::write(
        const std::string& fileName, PartBunch_t& bunch, std::size_t stepSizeSegment,
        unsigned long long stepsCompletedInSegment) {
    namespace fs                    = std::filesystem;
    const std::string temporaryName = fileName + ".tmp";

    if (ippl::Comm->rank() == 0) {
        std::error_code ec;
        fs::remove(temporaryName, ec);
    }
    ippl::Comm->barrier();

    h5_file_t file = openFile(temporaryName, H5_O_WRONLY);
    try {
        h5_int64_t format       = checkpointFormatVersion;
        h5_int64_t step         = static_cast<h5_int64_t>(bunch.getGlobalTrackStep());
        h5_int64_t count        = static_cast<h5_int64_t>(bunch.getNumParticleContainers());
        h5_int64_t segment      = static_cast<h5_int64_t>(stepSizeSegment);
        h5_int64_t segmentSteps = static_cast<h5_int64_t>(stepsCompletedInSegment);
        h5_float64_t time       = bunch.getT();
        h5_float64_t dt         = bunch.getdT();
        requireSuccess(
                H5WriteFileAttribInt64(file, "OPALXCheckpointFormat", &format, 1),
                "writing checkpoint format version");
        requireSuccess(
                H5WriteFileAttribInt64(file, "GlobalTrackStep", &step, 1),
                "writing global tracking step");
        requireSuccess(H5WriteFileAttribFloat64(file, "TIME", &time, 1), "writing checkpoint time");
        requireSuccess(H5WriteFileAttribFloat64(file, "DT", &dt, 1), "writing checkpoint dt");
        requireSuccess(
                H5WriteFileAttribInt64(file, "NumContainers", &count, 1),
                "writing particle-container count");
        requireSuccess(
                H5WriteFileAttribInt64(file, "StepSizeSegment", &segment, 1),
                "writing step-size segment");
        requireSuccess(
                H5WriteFileAttribInt64(file, "StepsCompletedInSegment", &segmentSteps, 1),
                "writing completed steps in segment");

        h5_int64_t phaseCount = OpalData::getInstance()->getNumberOfMaxPhases();
        requireSuccess(
                H5WriteFileAttribInt64(file, "NumCavityPhases", &phaseCount, 1),
                "writing cavity phase count");
        std::size_t phaseIndex = 0;
        for (auto it = OpalData::getInstance()->getFirstMaxPhases();
             it != OpalData::getInstance()->getLastMaxPhases(); ++it, ++phaseIndex) {
            const std::string prefix = "CavityPhase" + std::to_string(phaseIndex);
            requireSuccess(
                    H5WriteFileAttribString(file, (prefix + "Name").c_str(), it->first.c_str()),
                    "writing cavity phase name");
            h5_float64_t phase = it->second;
            requireSuccess(
                    H5WriteFileAttribFloat64(file, (prefix + "Value").c_str(), &phase, 1),
                    "writing cavity phase value");
        }

        for (std::size_t ci = 0; ci < bunch.getNumParticleContainers(); ++ci) {
            const auto pc = bunch.getParticleContainer(ci);
            if (!pc) {
                throw OpalException(
                        "CheckpointFile::write",
                        "particle container " + std::to_string(ci) + " is null");
            }

            requireSuccess(
                    H5SetStep(file, static_cast<h5_int64_t>(ci)),
                    "selecting checkpoint container record");

            h5_int64_t containerIndex = static_cast<h5_int64_t>(ci);
            h5_int64_t species        = static_cast<h5_int64_t>(pc->Sp);
            h5_int64_t qmMode =
                    pc->getQMStorageMode() == decltype(pc->getQMStorageMode())::Attributes;
            h5_int64_t spinEnabled = pc->hasSpin() ? 1 : 0;
            h5_int64_t atSStop     = bunch.pcAtSStop(ci) ? 1 : 0;
            h5_float64_t spos      = pc->get_sPos();
            // Attribute-mode Q/M can vary per particle and some ranks may be empty. The particle
            // datasets are authoritative in that mode; use rank-independent placeholders here.
            h5_float64_t macroQ        = qmMode != 0 ? 0.0 : pc->getChargePerParticle();
            h5_float64_t macroM        = qmMode != 0 ? 0.0 : pc->getMassPerParticle();
            Vector_t<double, 3> refR   = pc->getRefPartR();
            Vector_t<double, 3> refP   = pc->getRefPartP();
            Vector_t<double, 3> origin = pc->getToLabTrafo().getOrigin();
            Quaternion rotation        = pc->getToLabTrafo().getRotation();
            h5_float64_t quaternion[4] = {rotation(0), rotation(1), rotation(2), rotation(3)};

            requireSuccess(
                    H5WriteStepAttribInt64(file, "ContainerIndex", &containerIndex, 1),
                    "writing container index");
            requireSuccess(
                    H5WriteStepAttribInt64(file, "Species", &species, 1),
                    "writing particle species");
            requireSuccess(
                    H5WriteStepAttribInt64(file, "QMStorageMode", &qmMode, 1),
                    "writing Q/M storage mode");
            requireSuccess(
                    H5WriteStepAttribInt64(file, "SpinEnabled", &spinEnabled, 1),
                    "writing spin flag");
            requireSuccess(
                    H5WriteStepAttribInt64(file, "AtSStop", &atSStop, 1),
                    "writing step-size stop flag");
            requireSuccess(
                    H5WriteStepAttribFloat64(file, "SPOS", &spos, 1), "writing path position");
            requireSuccess(
                    H5WriteStepAttribFloat64(file, "MacroCharge", &macroQ, 1),
                    "writing macroparticle charge");
            requireSuccess(
                    H5WriteStepAttribFloat64(file, "MacroMass", &macroM, 1),
                    "writing macroparticle mass");
            requireSuccess(
                    H5WriteStepAttribFloat64(file, "RefPartR", &refR[0], 3),
                    "writing reference position");
            requireSuccess(
                    H5WriteStepAttribFloat64(file, "RefPartP", &refP[0], 3),
                    "writing reference momentum");
            requireSuccess(
                    H5WriteStepAttribFloat64(file, "ToLabOrigin", &origin[0], 3),
                    "writing lab-frame origin");
            requireSuccess(
                    H5WriteStepAttribFloat64(file, "ToLabQuaternion", quaternion, 4),
                    "writing lab-frame rotation");

            const std::size_t localCount = pc->getLocalNum();
            requireSuccess(
                    H5PartSetNumParticles(file, static_cast<h5_size_t>(localCount)),
                    "setting local checkpoint particle count");

            auto rHost              = hostCopy(pc->R.getView());
            auto pHost              = hostCopy(pc->P.getView());
            auto dtHost             = hostCopy(pc->dt.getView());
            auto idHost             = hostCopy(pc->ID.getView());
            auto binHost            = hostCopy(pc->Bin.getView());
            auto qHost              = hostCopy(pc->getQView());
            auto mHost              = hostCopy(pc->getMView());
            const bool qmAttributes = qmMode != 0;

            writeFloatDataset(file, "x", localCount, [&](std::size_t i) {
                return rHost(i)(0);
            });
            writeFloatDataset(file, "y", localCount, [&](std::size_t i) {
                return rHost(i)(1);
            });
            writeFloatDataset(file, "z", localCount, [&](std::size_t i) {
                return rHost(i)(2);
            });
            writeFloatDataset(file, "px", localCount, [&](std::size_t i) {
                return pHost(i)(0);
            });
            writeFloatDataset(file, "py", localCount, [&](std::size_t i) {
                return pHost(i)(1);
            });
            writeFloatDataset(file, "pz", localCount, [&](std::size_t i) {
                return pHost(i)(2);
            });
            writeFloatDataset(file, "dt", localCount, [&](std::size_t i) {
                return dtHost(i);
            });
            writeFloatDataset(file, "q", localCount, [&](std::size_t i) {
                return qHost(qmAttributes ? i : 0);
            });
            writeFloatDataset(file, "m", localCount, [&](std::size_t i) {
                return mHost(qmAttributes ? i : 0);
            });
            writeInt64Dataset(file, "id", localCount, [&](std::size_t i) {
                return idHost(i);
            });
            writeInt32Dataset(file, "bin", localCount, [&](std::size_t i) {
                return binHost(i);
            });

            if (pc->hasSpin()) {
                auto polHost = hostCopy(pc->Pol.getView());
                writeFloatDataset(file, "polx", localCount, [&](std::size_t i) {
                    return static_cast<double>(polHost(i)(0));
                });
                writeFloatDataset(file, "poly", localCount, [&](std::size_t i) {
                    return static_cast<double>(polHost(i)(1));
                });
                writeFloatDataset(file, "polz", localCount, [&](std::size_t i) {
                    return static_cast<double>(polHost(i)(2));
                });
            }
        }
        closeFile(file);
    } catch (...) {
        if (file != 0) {
            H5CloseFile(file);
            file = 0;
        }
        throw;
    }

    int renameSucceeded = 1;
    std::string renameError;
    if (ippl::Comm->rank() == 0) {
        std::error_code ec;
        fs::rename(temporaryName, fileName, ec);
        if (ec) {
            renameSucceeded = 0;
            renameError     = ec.message();
        }
    }
    MPI_Bcast(&renameSucceeded, 1, MPI_INT, 0, ippl::Comm->getCommunicator());
    if (!renameSucceeded) {
        throw OpalException(
                "CheckpointFile::write",
                "could not atomically replace checkpoint '" + fileName + "': " + renameError);
    }
    ippl::Comm->barrier();
}

CheckpointFile::Metadata CheckpointFile::inspect(const std::string& fileName) {
    h5_file_t file = openFile(fileName, H5_O_RDONLY);
    try {
        Metadata metadata = readMetadata(file);
        closeFile(file);
        return metadata;
    } catch (...) {
        if (file != 0) {
            H5CloseFile(file);
        }
        throw;
    }
}

CheckpointFile::Metadata CheckpointFile::read(const std::string& fileName, PartBunch_t& bunch) {
    h5_file_t file = openFile(fileName, H5_O_RDONLY);
    try {
        const Metadata metadata = readMetadata(file);
        restoreCavityPhases(file);
        std::vector<bool> atSStop(metadata.numContainers, false);
        if (metadata.numContainers != bunch.getNumParticleContainers()) {
            throw OpalException(
                    "CheckpointFile::read",
                    "checkpoint contains " + std::to_string(metadata.numContainers)
                            + " particle containers, but the input defines "
                            + std::to_string(bunch.getNumParticleContainers()));
        }

        for (std::size_t ci = 0; ci < metadata.numContainers; ++ci) {
            auto pc = bunch.getParticleContainer(ci);
            if (!pc || pc->getLocalNum() != 0) {
                throw OpalException(
                        "CheckpointFile::read",
                        "target particle container " + std::to_string(ci) + " is not empty");
            }

            requireSuccess(
                    H5SetStep(file, static_cast<h5_int64_t>(ci)),
                    "selecting checkpoint container record");
            requireStepAttribute(file, "ContainerIndex");
            requireStepAttribute(file, "Species");
            requireStepAttribute(file, "QMStorageMode");
            requireStepAttribute(file, "SpinEnabled");
            requireStepAttribute(file, "AtSStop");
            requireStepAttribute(file, "SPOS");
            requireStepAttribute(file, "MacroCharge");
            requireStepAttribute(file, "MacroMass");
            requireStepAttribute(file, "RefPartR");
            requireStepAttribute(file, "RefPartP");
            requireStepAttribute(file, "ToLabOrigin");
            requireStepAttribute(file, "ToLabQuaternion");

            h5_int64_t containerIndex = -1;
            h5_int64_t species        = 0;
            h5_int64_t qmMode         = 0;
            h5_int64_t spinEnabled    = 0;
            h5_int64_t savedAtSStop   = 0;
            h5_float64_t spos         = 0.0;
            h5_float64_t macroQ       = 0.0;
            h5_float64_t macroM       = 0.0;
            Vector_t<double, 3> refR(0.0);
            Vector_t<double, 3> refP(0.0);
            Vector_t<double, 3> origin(0.0);
            h5_float64_t quaternion[4] = {1.0, 0.0, 0.0, 0.0};
            requireSuccess(
                    H5ReadStepAttribInt64(file, "ContainerIndex", &containerIndex),
                    "reading container index");
            requireSuccess(
                    H5ReadStepAttribInt64(file, "Species", &species), "reading particle species");
            requireSuccess(
                    H5ReadStepAttribInt64(file, "QMStorageMode", &qmMode),
                    "reading Q/M storage mode");
            requireSuccess(
                    H5ReadStepAttribInt64(file, "SpinEnabled", &spinEnabled), "reading spin flag");
            requireSuccess(
                    H5ReadStepAttribInt64(file, "AtSStop", &savedAtSStop),
                    "reading step-size stop flag");
            requireSuccess(H5ReadStepAttribFloat64(file, "SPOS", &spos), "reading path position");
            requireSuccess(
                    H5ReadStepAttribFloat64(file, "MacroCharge", &macroQ),
                    "reading macroparticle charge");
            requireSuccess(
                    H5ReadStepAttribFloat64(file, "MacroMass", &macroM),
                    "reading macroparticle mass");
            requireSuccess(
                    H5ReadStepAttribFloat64(file, "RefPartR", &refR[0]),
                    "reading reference position");
            requireSuccess(
                    H5ReadStepAttribFloat64(file, "RefPartP", &refP[0]),
                    "reading reference momentum");
            requireSuccess(
                    H5ReadStepAttribFloat64(file, "ToLabOrigin", &origin[0]),
                    "reading lab-frame origin");
            requireSuccess(
                    H5ReadStepAttribFloat64(file, "ToLabQuaternion", quaternion),
                    "reading lab-frame rotation");

            const bool targetQMAttributes =
                    pc->getQMStorageMode() == decltype(pc->getQMStorageMode())::Attributes;
            if (species != static_cast<h5_int64_t>(pc->Sp)) {
                throw OpalException(
                        "CheckpointFile::read",
                        "particle species in the input does not match container "
                                + std::to_string(ci) + " in the checkpoint");
            }
            if ((qmMode != 0) != targetQMAttributes) {
                throw OpalException(
                        "CheckpointFile::read",
                        "QM_MODE in the input does not match particle container "
                                + std::to_string(ci) + " in the checkpoint");
            }
            if ((spinEnabled != 0) != pc->hasSpin()) {
                throw OpalException(
                        "CheckpointFile::read",
                        "spin configuration in the input does not match particle container "
                                + std::to_string(ci) + " in the checkpoint");
            }
            if (!targetQMAttributes
                && (!compatibleScalar(pc->getChargePerParticle(), macroQ)
                    || !compatibleScalar(pc->getMassPerParticle(), macroM))) {
                throw OpalException(
                        "CheckpointFile::read",
                        "macroparticle charge or mass in the input does not match container "
                                + std::to_string(ci) + " in the checkpoint");
            }
            if (savedAtSStop != 0 && savedAtSStop != 1) {
                throw OpalException(
                        "CheckpointFile::read", "invalid step-size stop flag in checkpoint");
            }
            atSStop[ci] = savedAtSStop != 0;
            if (containerIndex != static_cast<h5_int64_t>(ci)) {
                throw OpalException(
                        "CheckpointFile::read", "checkpoint container records are out of order");
            }

            const h5_ssize_t globalCount = H5PartGetNumParticles(file);
            if (globalCount < 0) {
                throw OpalException("CheckpointFile::read", "could not read global particle count");
            }
            const std::size_t ranks      = static_cast<std::size_t>(ippl::Comm->size());
            const std::size_t rank       = static_cast<std::size_t>(ippl::Comm->rank());
            const std::size_t base       = static_cast<std::size_t>(globalCount) / ranks;
            const std::size_t extra      = static_cast<std::size_t>(globalCount) % ranks;
            const std::size_t localCount = base + (rank < extra ? 1 : 0);

            requireSuccess(
                    H5PartSetNumParticles(file, static_cast<h5_size_t>(localCount)),
                    "selecting this rank's checkpoint particle view");
            pc->createParticles(localCount);

            const auto x   = readFloatDataset(file, "x", localCount);
            const auto y   = readFloatDataset(file, "y", localCount);
            const auto z   = readFloatDataset(file, "z", localCount);
            const auto px  = readFloatDataset(file, "px", localCount);
            const auto py  = readFloatDataset(file, "py", localCount);
            const auto pz  = readFloatDataset(file, "pz", localCount);
            const auto dts = readFloatDataset(file, "dt", localCount);
            const auto q   = readFloatDataset(file, "q", localCount);
            const auto m   = readFloatDataset(file, "m", localCount);
            const auto id  = readInt64Dataset(file, "id", localCount);
            const auto bin = readInt32Dataset(file, "bin", localCount);

            auto rHost   = pc->R.getHostMirror();
            auto pHost   = pc->P.getHostMirror();
            auto dtHost  = pc->dt.getHostMirror();
            auto idHost  = pc->ID.getHostMirror();
            auto binHost = pc->Bin.getHostMirror();
            for (std::size_t i = 0; i < localCount; ++i) {
                rHost(i)   = Vector_t<double, 3>(x[i], y[i], z[i]);
                pHost(i)   = Vector_t<double, 3>(px[i], py[i], pz[i]);
                dtHost(i)  = dts[i];
                idHost(i)  = static_cast<std::int64_t>(id[i]);
                binHost(i) = static_cast<typename PartBunch_t::binIndex_t>(bin[i]);
            }
            Kokkos::deep_copy(pc->R.getView(), rHost);
            Kokkos::deep_copy(pc->P.getView(), pHost);
            Kokkos::deep_copy(pc->dt.getView(), dtHost);
            Kokkos::deep_copy(pc->ID.getView(), idHost);
            Kokkos::deep_copy(pc->Bin.getView(), binHost);

            if (targetQMAttributes) {
                auto qHost = Kokkos::create_mirror_view(pc->getQView());
                auto mHost = Kokkos::create_mirror_view(pc->getMView());
                for (std::size_t i = 0; i < localCount; ++i) {
                    qHost(i) = q[i];
                    mHost(i) = m[i];
                }
                Kokkos::deep_copy(pc->getQView(), qHost);
                Kokkos::deep_copy(pc->getMView(), mHost);
            } else {
                pc->setQ(macroQ);
                pc->setM(macroM);
            }

            if (pc->hasSpin()) {
                const auto polx = readFloatDataset(file, "polx", localCount);
                const auto poly = readFloatDataset(file, "poly", localCount);
                const auto polz = readFloatDataset(file, "polz", localCount);
                auto polHost    = pc->Pol.getHostMirror();
                for (std::size_t i = 0; i < localCount; ++i) {
                    polHost(i) = typename PartBunch_t::ParticleContainer_t::spin_vector_type{
                            static_cast<float>(polx[i]), static_cast<float>(poly[i]),
                            static_cast<float>(polz[i])};
                }
                Kokkos::deep_copy(pc->Pol.getView(), polHost);
            }

            pc->set_sPos(spos);
            pc->setRefPartR(refR);
            pc->setRefPartP(refP);
            pc->setToLabTrafo(CoordinateSystemTrafo(
                    origin,
                    Quaternion(quaternion[0], quaternion[1], quaternion[2], quaternion[3])));
            pc->markMomentsDirty();
            requireSuccess(H5PartResetView(file), "resetting checkpoint particle view");
        }

        bunch.setT(metadata.time);
        bunch.setdT(metadata.dt);
        bunch.setGlobalTrackStep(metadata.globalTrackStep);
        bunch.resetPcActive();
        for (std::size_t ci = 0; ci < atSStop.size(); ++ci) {
            if (atSStop[ci]) {
                bunch.setPcAtSStop(ci);
            }
        }
        Kokkos::fence();
        closeFile(file);
        return metadata;
    } catch (...) {
        if (file != 0) {
            H5CloseFile(file);
        }
        throw;
    }
}
