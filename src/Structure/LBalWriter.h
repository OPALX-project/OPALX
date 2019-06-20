#ifndef OPAL_LBAL_WRITER_H
#define OPAL_LBAL_WRITER_H

#include "SDDSWriter.h"

class LBalWriter : public SDDSWriter {
    
private:
    void writeHeader(PartBunchBase<double, 3> *beam);
    void writeData(PartBunchBase<double, 3> *beam);
}
#endif
