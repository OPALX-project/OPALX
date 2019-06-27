#ifndef OPAL_GRID_LBAL_WRITER_H
#define OPAL_GRID_LBAL_WRITER_H

#include "SDDSWriter.h"

class GridLBalWriter : public SDDSWriter {

private:
    void fillHeader(PartBunchBase<double, 3> *beam);
};

#endif