#ifndef OPAL_LBAL_WRITER_H
#define OPAL_LBAL_WRITER_H

#include "SDDSWriter.h"

class LBalWriter : public SDDSWriter {

public:
    LBalWriter(const std::string& fname, bool restart);
    
    void write(PartBunchBase<double, 3> *beam);
    
private:
    void fillHeader_m(PartBunchBase<double, 3> *beam);
};


#endif
