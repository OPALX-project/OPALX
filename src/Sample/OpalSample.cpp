#include "Sample/OpalSample.h"

#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Utilities/OpalException.h"
#include "Utilities/Util.h"

#include "Sample/Uniform.h"


// Class OpalSample
// ------------------------------------------------------------------------

// The attributes of class OpalSample.
namespace {
    enum {
        TYPE,       // The type of sampling
        PARAMETERS, //
        SIZE
    };
}

OpalSample::OpalSample():
    Definition(SIZE, "SAMPLING",
               "The \"SAMPLING\" statement defines methods used for the optimizer in sample mode.")
{
    itsAttr[TYPE]       = Attributes::makeString
                          ("TYPE", "UNIFORM_INT, UNIFORM_REAL, SEQUENCE, FROMFILE");

    itsAttr[PARAMETERS] = Attributes::makeRealArray
                         ("PARAMETERS", "List of parameters (depends on TYPE)");


    registerOwnership(AttributeHandler::STATEMENT);

    OpalSample *defSampling = clone("UNNAMED_SAMPLING");
    defSampling->builtin = true;

    try {
        defSampling->update();
        OpalData::getInstance()->define(defSampling);
    } catch(...) {
        delete defSampling;
    }
}


OpalSample::OpalSample(const std::string &name, OpalSample *parent):
    Definition(name, parent)
{}


bool OpalSample::canReplaceBy(Object *object) {
    return dynamic_cast<OpalSample *>(object) != nullptr;
}


OpalSample *OpalSample::clone(const std::string &name) {
    return new OpalSample(name, this);
}


void OpalSample::execute() {
    update();
}


OpalSample *OpalSample::find(const std::string &name) {
    OpalSample *sampling = dynamic_cast<OpalSample *>(OpalData::getInstance()->find(name));

    if (sampling == nullptr) {
        throw OpalException("OpalSample::find()", "OpalSample \"" + name + "\" not found.");
    }
    return sampling;
}


void OpalSample::update() {
    // Set default name.
    if (getOpalName().empty()) setOpalName("UNNAMED_SAMPLING");
}


void OpalSample::initOpalSample() {
    
    std::string type = Util::toUpper(Attributes::getString(itsAttr[TYPE]));
        
    if (type == "UNIFORM_INT") {
        std::vector<double> params = Attributes::getRealArray(itsAttr[PARAMETERS]);
        
        if ( params.size() != 3 ) {
            throw OpalException("OpalSample::initOpalSample()",
                                type + " requires exactly 3 arguments.");
        }
        
        /* params[0]: lower bound
         * params[1]: upper bound
         * params[2]: seed
         */
        
        if ( params[0] >= params[1] ) {
            throw OpalException("OpalSample::initOpalSample()",
                                "Lower bound >= upper bound.");
        }
        
        sampleMethod_m.reset( new Uniform<int>(params[0], params[1], params[2]) );
            
            
    } else if (type == "UNIFORM_REAL") {
        std::vector<double> params = Attributes::getRealArray(itsAttr[PARAMETERS]);
    
        if ( params.size() != 3 ) {
            throw OpalException("OpalSample::initOpalSample()",
                                type + " requires exactly 3 arguments.");
        }
        
        if ( params[0] >= params[1] ) {
            throw OpalException("OpalSample::initOpalSample()",
                                "Lower bound >= upper bound.");
        }
            
        /* params[0]: lower bound
         * params[1]: upper bound
         * params[2]: seed
         */
        sampleMethod_m.reset( new Uniform<double>(params[0], params[1], params[2]) );
            
    } else {
        throw OpalException("OpalSample::initOpalSample()",
                            "Unkown sampling method: '" + type + "'.");
    }
}
