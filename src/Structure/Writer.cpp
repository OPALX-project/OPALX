#include "Writer.h"

#include <boost/filesystem.hpp>

Writer::Writer(const std::string& fname)
    : fname_m(fname)
{ }


bool Writer::exists() {
    return boost::filesystem::exists(fname_m);
}


void Writer::rewindLines(size_t numberOfLines) {
    
}
