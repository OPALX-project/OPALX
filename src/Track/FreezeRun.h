//
// freeze(): copy the values one RUN command tracks with out of the parsed objects.
//
// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
// All rights reserved
//
#ifndef OPALX_FREEZE_RUN_H
#define OPALX_FREEZE_RUN_H

#include "Track/RunDescription.h"

class BeamSequence;
class Object;

namespace opalx::run {

    /**
     * @brief Copy everything the RUN command @p runCommand tracks with into a RunDescription.
     *
     * Takes the TRACK settings @p track and reads the RUN attributes, the beams, their emission
     * sources and distributions, the field solver, the tracked line and
     * the OPTION values. Every attribute is evaluated now, so ":=" expressions take the values
     * valid at this RUN. An object that cannot be found is left out; validate() reports that error
     * itself.
     *
     * @param runCommand The RUN command being executed.
     * @param track      The settings of the open TRACK (TrackCmd::openTrack()).
     * @param runIndex   1 for the first RUN of the input, 2 for the second, ...
     * @return The run description.
     */
    RunDescription freeze(const Object& runCommand, const TrackDescription& track, int runIndex);

    /**
     * @brief Describe one parsed object: its keyword, its name and every attribute, evaluated now.
     * @note An attribute that fails to evaluate is kept with its error message, so describing an
     *       object never throws for an attribute the run does not read.
     */
    ObjectDescription describeObject(const Object& object);

    /**
     * @brief Describe a LINE: its start position and direction, and every element occurrence
     *        in line order, nested lines entered.
     */
    LineDescription describeLine(const BeamSequence& sequence);

}  // namespace opalx::run

#endif  // OPALX_FREEZE_RUN_H
