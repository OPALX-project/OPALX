// ------------------------------------------------------------------------
// $RCSfile: AttributeHandler.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: AttributeHandler
//   An abstract class used to parse and print attributes.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:34 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/AttributeHandler.h"
#include "Parser/Statement.h"
#include "Utilities/OpalException.h"


// Class AttributeHandler
// ------------------------------------------------------------------------

AttributeHandler::AttributeHandler
(const string &name, const string &help, AttributeBase *def):
    RCObject(), itsName(name), itsHelp(help), itsDefault(def),
    is_deferred(false), is_readonly(false)
{}


AttributeHandler::~AttributeHandler()
{}


AttributeHandler *AttributeHandler::clone() const {
    throw OpalException("AttributeHandler::clone()",
                        "Internal error: should not call this method.");
}


AttributeBase *AttributeHandler::getDefault() const {
    if(itsDefault.isValid()) {
        return &*itsDefault;
    } else {
        throw OpalException("AttributeHandler::getDefault()",
                            "Attribute \"" + itsName + "\" has no default value.");
    }
}


const string &AttributeHandler::getHelp() const {
    return itsHelp;
}


const string &AttributeHandler::getName() const {
    return itsName;
}


void AttributeHandler::parseComponent
(Attribute &, Statement &, bool, int) const {
    // Default behaviour.
    throw OpalException("AttributeHandler::parseComponent()",
                        "You cannot assign to a component of \"" + itsName +
                        "\" which is not a vector value.");
}


bool AttributeHandler::isDeferred() const {
    return is_deferred;
}


void AttributeHandler::setDeferred(bool flag) {
    is_deferred = flag;
}


bool AttributeHandler::isReadOnly() const {
    return is_readonly;
}


void AttributeHandler::setReadOnly(bool flag) {
    is_readonly = flag;
}

