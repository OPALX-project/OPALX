#ifndef __OPTPILOTEXCEPTION_H__
#define __OPTPILOTEXCEPTION_H__

#include <string>

class OptPilotException {

public:

    OptPilotException(const std::string &meth, const std::string &descr) {
        descr_ = descr;
        meth_ = meth;
    }

    virtual const char* what() const throw() {
        return descr_.c_str();
    }

    virtual const char* where() const throw() {
        return meth_.c_str();
    }

private:

    std::string descr_;
    std::string meth_;

};

#endif