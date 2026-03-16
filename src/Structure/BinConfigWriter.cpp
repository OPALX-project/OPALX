//
// Class BinConfigWriter
//   Writes binning configuration snapshots to a JSON file.
//
// See the corresponding header for copyright and license details.
//

#include "Structure/BinConfigWriter.h"

#include "Utilities/OpalException.h"
#include "Ippl.h"

#include <iomanip>

BinConfigWriter::BinConfigWriter(const std::string& fileName)
    : os_(fileName.c_str(), std::ios::out | std::ios::trunc)
    , first_(true) {
    Inform m("BinConfigWriter::BinConfigWriter");

    if (!os_) {
        m << level4 << "Failed to open bin configuration file \"" << fileName << "\"." << endl;
        throw OpalException("BinConfigWriter::BinConfigWriter",
                            "Failed to open bin configuration file \"" + fileName + "\"");
    }

    m << level4 << "Opened bin configuration JSON file \"" << fileName << "\" for writing." << endl;
    os_ << "[\n";
}

BinConfigWriter::~BinConfigWriter() {
    Inform m("BinConfigWriter::~BinConfigWriter");

    if (os_) {
        m << level5 << "Closing bin configuration JSON file." << endl;
        os_.close();
    }
}

void BinConfigWriter::writeEntry(long long step,
                                 double time,
                                 bool preMerge,
                                 const std::vector<std::size_t>& binCounts,
                                 const std::vector<double>& binWidths,
                                 double xMin) {
    if (!os_) {
        return;
    }

    Inform m("BinConfigWriter::writeEntry");
    m << level5
      << "Writing bin configuration entry: step=" << step
      << ", time=" << std::setprecision(17) << time
      << ", preMerge=" << (preMerge ? 1 : 0)
      << ", nBins=" << binCounts.size()
      << ", xMin=" << std::setprecision(17) << xMin << endl;
    
    // We maintain a *closed* top-level JSON array on disk after each call:
    //   [ <obj1>,\n  <obj2>,\n ... <objN>\n]
    //
    // For the first entry we just append after "[\n".
    // For subsequent entries we seek back over the trailing "\n]\n" and
    // insert ",\n" before appending the next object and the new "\n]\n".

    // Move write position to the end of the file.
    os_.seekp(0, std::ios::end);

    if (!first_) {
        // File currently ends with "\n]\n" (3 chars). Step back before that
        // closing sequence and append a comma+newline to separate entries.
        os_.seekp(-3, std::ios::end);
        os_ << ",\n";
    } else {
        first_ = false;
    }

    os_ << "  {\n";
    os_ << "    \"step\": " << step << ",\n";
    os_ << "    \"time\": " << std::setprecision(17) << time << ",\n";
    os_ << "    \"preMerge\": " << (preMerge ? "true" : "false") << ",\n";
    os_ << "    \"xMin\": " << std::setprecision(17) << xMin << ",\n";

    // binCounts array (pretty-printed, fixed number of entries per line).
    os_ << "    \"binCounts\": [";
    for (std::size_t i = 0; i < binCounts.size(); ++i) {
        if (i % 16 == 0) {
            os_ << "\n      ";
        }
        os_ << binCounts[i];
        if (i + 1 < binCounts.size()) {
            os_ << ", ";
        }
    }
    if (!binCounts.empty()) {
        os_ << "\n";
    }
    os_ << "    ],\n";

    // binWidths array (pretty-printed as well).
    os_ << "    \"binWidths\": [";
    for (std::size_t i = 0; i < binWidths.size(); ++i) {
        if (i % 16 == 0) {
            os_ << "\n      ";
        }
        os_ << std::setprecision(17) << binWidths[i];
        if (i + 1 < binWidths.size()) {
            os_ << ", ";
        }
    }
    if (!binWidths.empty()) {
        os_ << "\n";
    }
    os_ << "    ]\n";

    os_ << "  ";
    os_ << "}";

    // Close the top-level array so the file is valid JSON after each write.
    os_ << "\n]\n";

    // Flush each entry so that the JSON stays well-formed on disk even if the
    // simulation terminates unexpectedly between steps (between calls, the
    // file always ends with a complete, closed JSON array).
    os_.flush();
}

