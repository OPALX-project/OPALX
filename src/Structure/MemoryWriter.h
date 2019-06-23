#ifndef OPAL_MEMORY_WRITER_H
#define OPAL_MEMORY_WRITER_H

#include "SDDSWriter.h"

class MemoryWriter : public SDDSWriter {
    
public:
    MemoryWriter(const std::string& fname, bool restart);
    
private:
    void fillHeader_m();
};

#endif
