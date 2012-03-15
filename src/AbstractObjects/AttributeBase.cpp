// ------------------------------------------------------------------------
// $RCSfile: AttributeBase.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: AttributeBase
//   The abstract base class for all attribute types.
//
// ------------------------------------------------------------------------
//
// $Date: 2002/12/09 15:06:07 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/AttributeBase.h"
#if defined(__GNUC__) && __GNUC__ < 3
#include <strstream>
#else
#include <sstream>
#endif

// Class AttributeBase
// ------------------------------------------------------------------------

AttributeBase::AttributeBase():
    RCObject()
{}


AttributeBase::~AttributeBase()
{}


string AttributeBase::getImage() const {
#if defined(__GNUC__) && __GNUC__ < 3
    std::ostrstream os;
#else
    std::ostringstream os;
#endif
    print(os);
    os << std::ends;
#if defined(__GNUC__) && __GNUC__ < 3
    string img(os.str());
    os.freeze(0);
    return img;
#else
    return os.str();
#endif
}


bool AttributeBase::isExpression() const {
    return false;
}
