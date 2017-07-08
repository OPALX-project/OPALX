#ifndef __TIMESTAMP_H__
#define __TIMESTAMP_H__

#include <sstream>
#include <boost/chrono.hpp>
#include <ctime>

#include "Util/Trace/TraceComponent.h"

#include "mpi.h"

class Timestamp : public TraceComponent {

public:


    Timestamp()
        : TraceComponent("Timestamp")
    {}

    void execute(std::ostringstream &dump) {

        boost::chrono::time_point<boost::chrono::system_clock> now;
        now = boost::chrono::system_clock::now();
        std::time_t now_time = boost::chrono::system_clock::to_time_t(now);

        std::ostringstream timestamp;
        timestamp << std::ctime(&now_time);

        prepend(dump, timestamp);
    }

    virtual ~Timestamp()
      {}

};

#endif