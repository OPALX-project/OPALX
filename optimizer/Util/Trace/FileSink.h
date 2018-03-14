#ifndef __TRACE_FILE_SINK_H__
#define __TRACE_FILE_SINK_H__

#include <sstream>
#include <iostream>
#include <fstream>

#include "Util/Trace/TraceComponent.h"

class FileSink : public TraceComponent {

public:


    FileSink(std::string filename)
        : TraceComponent("FileSink")
        , filename_(filename)
    {}

    virtual ~FileSink()
    {}

    void execute(std::ostringstream &dump) {
        std::ofstream file;
        file.open(filename_.c_str(), std::ios::app);
        file << dump.str() << std::flush;
        file.close();
    }

private:

    std::string filename_;

};

#endif