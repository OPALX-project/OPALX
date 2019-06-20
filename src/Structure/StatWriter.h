#ifndef OPAL_STAT_WRITER_H
#define OPAL_STAT_WRITER_H

#include "SDDSWriter.h"

class StatWriter : public SDDSWriter {

private:
    void writeHeader();
    
    void writeData(PartBunchBase<double, 3> *beam);
};

#endif
