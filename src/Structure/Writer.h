#ifndef OPAL_OUTPUT_WRITER_H
#define OPAL_OUTPUT_WRITER_H

#include "Algorithms/PartBunchBase.h"

class Writer {
public:
    
    Writer(const std::string& fname);
    
    inline bool exists();
    
    virtual void write(PartBunchBase<double, 3> *beam) = 0;
    
    virtual void rewindLines(size_t numberOfLines);
    
    
    
protected:
    const std::string fname_m;
};

#endif
