#ifndef __TRACE_COMPONENT_H__
#define __TRACE_COMPONENT_H__

#include <string>
#include <sstream>

class TraceComponent {

public:

    TraceComponent(std::string name) : name_(name)
    {}

    ~TraceComponent()
    {}

    virtual void execute(std::ostringstream &dump) = 0;

    void prepend(std::ostringstream &dump, std::ostringstream &prepender) {

        prepender << dump.str();
        dump.str("");
        dump.clear();
        dump << prepender.str();
    }

    void prepend(std::ostringstream &dump, std::string prepender) {

        std::ostringstream tmp;
        tmp << prepender << dump.str();
        dump.str("");
        dump.clear();
        dump << tmp.str();
    }

private:

    std::string name_;

};

#endif
