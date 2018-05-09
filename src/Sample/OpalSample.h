#ifndef OPAL_SAMPLE_H
#define OPAL_SAMPLE_H

#include <string>
#include <memory>
#include "AbstractObjects/Definition.h"

#include "Sample/SamplingMethod.h"

// Class OpalSample
// ------------------------------------------------------------------------
/// The SAMPLING definition.
//  A SAMPLING definition is used run the optimizer in sample mode.

class OpalSample: public Definition {

public:

    /// Exemplar constructor.
    OpalSample();

    /// Test if replacement is allowed.
    //  Can replace only by another OpalSample
    virtual bool canReplaceBy(Object *object);

    /// Make clone.
    virtual OpalSample *clone(const std::string &name);

    /// Check the OpalSample data.
    virtual void execute();

    /// Find named trim coil.
    static OpalSample *find(const std::string &name);

    /// Update the OpalSample data.
    virtual void update();

    /// Initialise implementation
    void initOpalSample();
    
    std::shared_ptr<SamplingMethod> sampleMethod_m;

private:

    ///@{ Not implemented.
    OpalSample  (const OpalSample &) = delete;
    void operator=(const OpalSample &) = delete;
    ///@}
    /// Private copy constructor, called by clone
    OpalSample(const std::string &name, OpalSample *parent);

};

#endif
