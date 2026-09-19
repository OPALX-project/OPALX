/**
 * @file CheckpointFile.h
 * @brief Atomic HDF5 checkpoint I/O for restarting OPALX parallel tracking.
 */

#ifndef OPALX_CHECKPOINTFILE_H
#define OPALX_CHECKPOINTFILE_H

#include <cstddef>
#include <string>

#include "OPALTypes.h"

class CheckpointFile {
public:
    /// Metadata shared by all particle containers in one checkpoint.
    struct Metadata {
        long long globalTrackStep   = 0;  ///< Next integration step to execute.
        double time                 = 0.0;
        double dt                   = 0.0;
        std::size_t numContainers   = 0;
        std::size_t stepSizeSegment = 0;  ///< Segment containing the next integration step.
        unsigned long long stepsCompletedInSegment = 0;
    };

    /// Return the default checkpoint path for an OPAL input basename.
    static std::string defaultPath(const std::string& inputBasename);

    /**
     * Write a complete checkpoint through a same-directory temporary file and atomically replace
     * @p fileName only after every MPI rank has closed the temporary HDF5 file.
     */
    static void write(
            const std::string& fileName, PartBunch_t& bunch, std::size_t stepSizeSegment,
            unsigned long long stepsCompletedInSegment);

    /// Read metadata without changing a bunch.
    static Metadata inspect(const std::string& fileName);

    /**
     * Restore @p bunch from @p fileName. The saved global particle table for every container is
     * repartitioned evenly over the MPI ranks in the restart run.
     */
    static Metadata read(const std::string& fileName, PartBunch_t& bunch);
};

#endif  // OPALX_CHECKPOINTFILE_H
