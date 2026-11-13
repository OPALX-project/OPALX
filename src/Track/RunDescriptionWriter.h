//
// Text dump of a RunDescription, to compare what two OPALX versions read from the same input.
//
// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
// All rights reserved
//
#ifndef OPALX_RUN_DESCRIPTION_WRITER_H
#define OPALX_RUN_DESCRIPTION_WRITER_H

#include <iosfwd>
#include <string>

#include "Track/RunDescription.h"

namespace opalx::run {

    /**
     * @brief Write @p description as text: one value per line, in a fixed order.
     * @note Doubles are written with 17 significant digits, so equal text means equal values.
     */
    void write(const RunDescription& description, std::ostream& os);

    /// File name of the dump of one RUN: "<input base>_run<n>.txt", n = context.runIndex.
    std::string dumpFileName(const RunDescription& description);

}  // namespace opalx::run

#endif  // OPALX_RUN_DESCRIPTION_WRITER_H
