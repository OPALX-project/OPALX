#ifndef __TRACE_COUT_SINK_H__
#define __TRACE_COUT_SINK_H__

#include <iostream>
#include <string>

#include "Util/Trace/TraceComponent.h"

class CoutSink : public TraceComponent {

public:


    CoutSink(std::string prefix = "")
        : TraceComponent("CoutSink")
        , prefix_(prefix) {

        clear_color_ = "\e[0m";
    }

    ~CoutSink()
    {}


    void setColor(std::string color)      { color_ = color; }
    void setClearColor(std::string color) { clear_color_ = color; }


    void execute(std::ostringstream &dump) {
        std::cout << color_ << prefix_
                  << dump.str()
                  << clear_color_ << std::flush;
    }

private:

    std::string prefix_;
    std::string color_;
    std::string clear_color_;

};

#endif
