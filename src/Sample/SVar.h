#ifndef OPAL_SVAR_HH
#define OPAL_SVAR_HH

#include "AbstractObjects/Definition.h"

class SVar: public Definition {
public:
    SVar();
    ~SVar();

    virtual SVar *clone(const std::string &name);

    virtual void execute();
    
    std::string getVariable() const;
    
    double getLowerBound() const;
    
    double getUpperBound() const;
    
    std::string getType() const;

private:
    SVar(const std::string &name,
         SVar *parent);
};

inline
SVar* SVar::clone(const std::string &name) {
    return new SVar(name, this);
}

#endif