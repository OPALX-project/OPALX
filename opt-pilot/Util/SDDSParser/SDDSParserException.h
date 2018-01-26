#ifndef __SDDSPARSEREXCEPTION_H__
#define __SDDSPARSEREXCEPTION_H__

#include <string>

class SDDSParserException {

public:

    SDDSParserException(const std::string &meth, const std::string &descr) {
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