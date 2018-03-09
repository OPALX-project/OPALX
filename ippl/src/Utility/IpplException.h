#ifndef __IPPLEXCEPTION_H__
#define __IPPLEXCEPTION_H__

#include <string>

class IpplException {

public:

    IpplException(const std::string &meth, const std::string &descr) {
        descr_ = descr;
        meth_ = meth;
    }

    virtual const char* what() const throw() {
        return descr_.c_str();
    }

private:

    std::string descr_;
    std::string meth_;

};

#endif
