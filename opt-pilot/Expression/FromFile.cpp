#include "Expression/FromFile.h"

/// reads a simple list of double values
void FromFile::readValues() {

    values_.clear();

    std::ifstream file;
    file.open(filename_.c_str(), std::ios::in);
    if(!file) {
        throw OptPilotException("FromFile::readValues()",
                "Error opening file " + filename_);
    }

    std::copy(std::istream_iterator<double>(file),
              std::istream_iterator<double>(),
              std::back_inserter(values_));

    file.close();
}

