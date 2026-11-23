//
// Class H5PartWrapperForPT
//   A class that manages all calls to H5Part for the Parallel-T tracker.
//
// Copyright (c) 200x-2021, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#ifndef OPAL_H5PARTWRAPPERFORPT_H
#define OPAL_H5PARTWRAPPERFORPT_H

#include "OPALTypes.h"
#include "Structure/DataSinkOptions.h"
#include "Structure/H5PartWrapper.h"

#include "H5hut.h"

class H5PartWrapperForPT : public H5PartWrapper {
public:
    H5PartWrapperForPT(
            const std::string& fileName, const DataSinkOptions& options,
            h5_int32_t flags = H5_O_WRONLY);
    H5PartWrapperForPT(
            const std::string& fileName, const DataSinkOptions& options, int restartStep,
            std::string sourceFile, h5_int32_t flags = H5_O_RDWR);
    virtual ~H5PartWrapperForPT();

    virtual void writeHeader();
    virtual void writeStep(
            PartBunch_t*, const opalx::moments::BeamMoments& moments,
            const std::map<std::string, double>& additionalStepAttributes,
            size_t particleContainerIndex = 0);

    /// Moment groups the step header reads.
    static constexpr opalx::moments::Request requiredMoments =
            opalx::moments::Moment::Bounds | opalx::moments::Moment::MeanPosition
            | opalx::moments::Moment::MeanMomentum | opalx::moments::Moment::Spread;

private:
    void writeStepHeader(
            PartBunch_t*, const opalx::moments::BeamMoments& moments,
            const std::map<std::string, double>&, size_t particleContainerIndex);
    void writeStepData(PartBunch_t*, size_t particleContainerIndex);

    DataSinkOptions options_m;

    IpplTimings::TimerRef h5TotalStepTimer_m;
    IpplTimings::TimerRef h5StepHeaderTimer_m;
    IpplTimings::TimerRef h5DeviceHostCopyTimer_m;
    IpplTimings::TimerRef h5DatasetWriteTimer_m;
};

#endif  // OPAL_H5PARTWRAPPERFORPT_H
