#include "Sample/OpalSample.h"

#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Utilities/OpalException.h"
#include "Utilities/Util.h"

#include "Sample/Uniform.h"
#include "Sample/Normal.h"
#include "Sample/SampleSequence.h"
#include "Sample/SampleGaussianSequence.h"
#include "Sample/FromFile.h"


// Class OpalSample
// ------------------------------------------------------------------------

// The attributes of class OpalSample.
namespace {
    enum {
        TYPE,       // The type of sampling
        VARIABLE,   // name of design variable
        SEED,       // for random sample methods
        FNAME,      // file to read from sampling points
        N,
        SIZE
    };
}

OpalSample::OpalSample():
    Definition(SIZE, "SAMPLING",
               "The \"SAMPLING\" statement defines methods used for the optimizer in sample mode.")
    , size_m(1)
{
    itsAttr[TYPE]       = Attributes::makeString
                          ("TYPE", "UNIFORM_INT, UNIFORM_REAL, SEQUENCE, FROMFILE");

    itsAttr[VARIABLE]   = Attributes::makeString
                          ("VARIABLE", "Name of design variable");

    itsAttr[SEED]       = Attributes::makeReal
                          ("SEED", "seed for random sampling (default: 42)", 42);

    itsAttr[FNAME]      = Attributes::makeString
                          ("FNAME", "File to read from the sampling points");

    itsAttr[N]          = Attributes::makeReal
                          ("N", "Number of sampling points", 1);

    registerOwnership(AttributeHandler::STATEMENT);
}


OpalSample::OpalSample(const std::string &name, OpalSample *parent):
    Definition(name, parent)
{}


OpalSample *OpalSample::clone(const std::string &name) {
    return new OpalSample(name, this);
}


void OpalSample::execute() {

}


OpalSample *OpalSample::find(const std::string &name) {
    OpalSample *sampling = dynamic_cast<OpalSample *>(OpalData::getInstance()->find(name));

    if (sampling == nullptr) {
        throw OpalException("OpalSample::find()",
                            "OpalSample \"" + name + "\" not found.");
    }
    return sampling;
}


void OpalSample::initialize(const std::string &dvarName,
                            double lower,
                            double upper,
                            int nSample,
                            size_t modulo,
                            bool sequence) {

    if ( lower >= upper )
        throw OpalException("OpalSample::initOpalSample()",
                                "Lower bound >= upper bound.");

    std::string type = Util::toUpper(Attributes::getString(itsAttr[TYPE]));

    int seed = Attributes::getReal(itsAttr[SEED]);
    size_m = Attributes::getReal(itsAttr[N]);

    if (sequence) {
        if (type == "UNIFORM_INT") {
            sampleMethod_m.reset( new SampleSequence<int>(lower, upper, modulo, size_m) );
        } else if (type == "UNIFORM") {
            sampleMethod_m.reset( new SampleSequence<double>(lower, upper, modulo, size_m) );
        } else if (type == "GAUSSIAN") {
            sampleMethod_m.reset( new SampleGaussianSequence(lower, upper, modulo, size_m) );
        } else if (type == "FROMFILE") {
            std::string fname = Attributes::getString(itsAttr[FNAME]);
            sampleMethod_m.reset( new FromFile(fname, dvarName, modulo) );
            size_m = static_cast<FromFile*>(sampleMethod_m.get())->getSize();
       } else {
            throw OpalException("OpalSample::initOpalSample()",
                                "Unkown sampling method: '" + type + "'.");
        }
    } else {
        if (type == "UNIFORM_INT") {
            sampleMethod_m.reset( new Uniform<int>(lower, upper, seed) );
        } else if (type == "UNIFORM") {
            sampleMethod_m.reset( new Uniform<double>(lower, upper, seed) );
        } else if (type == "GAUSSIAN") {
            sampleMethod_m.reset( new Normal(lower, upper, seed) );
        } else if (type == "FROMFILE") {
            std::string fname = Attributes::getString(itsAttr[FNAME]);
            sampleMethod_m.reset( new FromFile(fname, dvarName, modulo) );
            size_m = static_cast<FromFile*>(sampleMethod_m.get())->getSize();
        } else {
            throw OpalException("OpalSample::initOpalSample()",
                                "Unkown sampling method: '" + type + "'.");
        }
    }
}


std::string OpalSample::getVariable() const {
    return Attributes::getString(itsAttr[VARIABLE]);
}