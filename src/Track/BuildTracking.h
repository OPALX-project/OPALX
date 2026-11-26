//
// buildTracking(): make the particle bunch, output, samplers, space-charge solver and step-size
// schedule of one RUN.
//
// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
// All rights reserved
//
#ifndef OPALX_BUILD_TRACKING_H
#define OPALX_BUILD_TRACKING_H

#include <memory>
#include <vector>

#include "Algorithms/TrackingSetup.h"
#include "PartBunch/PartBunch.h"
#include "SpaceCharge/SpaceChargeSolver.h"
#include "Structure/DataSink.h"
#include "Structure/H5PartWrapper.h"
#include "Track/RunDescription.h"

class ElementBase;

namespace opalx::run {

    /**
     * @brief The objects buildTracking() makes for one RUN; the tracker borrows them.
     *
     * Members are destroyed in reverse order: the space-charge solver before the bunch it
     * borrows particles from, then the data sink, which the solver borrows, and last the
     * phase-space files, which the data sink writes to.
     */
    struct RunObjects {
        /// Phase-space files, one per beam; the DataSink borrows them.
        std::vector<std::unique_ptr<H5PartWrapper>> phaseSpaceSinks;
        /// Statistics and phase-space output of this RUN.
        std::unique_ptr<DataSink> dataSink;
        std::unique_ptr<PartBunch_t> bunch;
        std::unique_ptr<spacecharge::SpaceChargeSolver> spaceChargeSolver;
    };

    /**
     * @brief Build everything the tracker needs for this RUN: beams, particle containers,
     *        output, samplers or the checkpoint, the space-charge solver, the step-size
     *        schedule and the stop rules.
     *
     * Reads only @p description, never the parsed input objects. The steps, numbered as in the
     * body:
     *  1. Print the field solver with its boundary conditions, the binning and the beams.
     *  2. For each beam, in TRACK BEAMS order: the beam whose reference particle it follows
     *     (its decay root), its emission sources (none for a decay daughter), its global
     *     processes, and the BunchBeam its particle container is made from, with the number of
     *     particles to reserve room for (BEAM NALLOC, or the sum of NPARTDIST).
     *  3. The space-charge configuration from the FIELDSOLVER, the BINNING, all emission
     *     sources and all distributions, and from it the domain the bunch is made on.
     *  4. The output: one .stat writer and, with ENABLEHDF5, one .h5 file per beam.
     *  5. The bunch: one particle container per beam. Each container gets the global processes
     *     of its beam, and decay daughters are connected to their parents.
     *  6. The start. A fresh run samples the distributions of every beam, starts at ZSTART (at
     *     0 if the first beam has no particles yet) and gives the bunch the DT of the segment
     *     that holds the start. A restart reads particles, reference particles, time, time step
     *     and step count from the checkpoint and rewinds the output files to it; it samples
     *     nothing.
     *  7. The space-charge solver, which allocates its fields and runs a warmup solve, also
     *     for TYPE = NONE.
     *  8. Print the start summary of each container.
     *  9. The TrackingSetup: the DT/ZSTOP/MAXSTEPS schedule, the start, the emitting samplers,
     *     the reference momenta from file distributions, the values a restart resumes with, the
     *     OPTION values, the DUMPEMFIELDS requests, RUN SCFIELDUPDATE and the stop rule of
     *     each beam.
     * 10. The beamline: every element gets the OPTION values and the run context, each field
     *     element is initialised and added, then all are placed, the lattice files written,
     *     the cavity phases of a checkpoint restored, the apertures checked against the reach
     *     and the element grid built.
     *
     * @param description What freeze() read for this RUN, checked by validate(); it names a
     *        FIELDSOLVER.
     * @param elements The elements buildElements() made for the tracked line.
     * @param objects Receives the bunch, the solver, the data sink and the phase-space files;
     *        it must outlive the tracker.
     * @return Everything the tracker reads; the pointers in it borrow from @p objects.
     * @throws OpalException for input that only this step can check: decay settings, sampler
     *         input, a checkpoint that cannot be read, and an aperture wider than the reach.
     * @throws LogicalError for a LASER in the tracked line.
     * @note The bunch time step must be set before step 10: an element may read it in
     *       initialise(). Step 6 sets it, and step 10 makes it positive first.
     */
    TrackingSetup buildTracking(
            const RunDescription& description, std::vector<std::shared_ptr<ElementBase>> elements,
            RunObjects& objects);

}  // namespace opalx::run

#endif  // OPALX_BUILD_TRACKING_H
