//
// Class H5Writer
//   Interface for H5 writers.
//
// Copyright (c) 2019-2020, Matthias Frey, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
#ifndef OPAL_H5_WRITER_H
#define OPAL_H5_WRITER_H

#include "H5hut.h"
#include "PartBunch/PartBunch.h"
#include "Structure/H5PartWrapper.h"

class H5Writer {
public:
    explicit H5Writer(H5PartWrapper* h5wrapper);

    void rewindToTime(double checkpointTime, double timeStep);

    /** \brief Dumps Phase Space to H5 file.
     *
     * \param beam The beam.
     * \param FDext The external fields at the reference particle: FDext[0] = B, FDext[1] = E.
     */
    void writePhaseSpace(
            PartBunch_t& beam, const opalx::moments::BeamMoments& moments,
            Vector_t<double, 3> FDext[], size_t particleContainerIndex = 0);

private:
    H5PartWrapper* h5wrapper_m;
};

inline void H5Writer::rewindToTime(double checkpointTime, double timeStep) {
    h5wrapper_m->rewindToTime(checkpointTime, timeStep);
}

#endif
