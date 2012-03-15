// ------------------------------------------------------------------------
// $RCSfile: OpalMarker.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.4 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalMarker
//   The class of OPAL markers.
//
// ------------------------------------------------------------------------
//
// $Date: 2002/01/22 15:16:02 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "Elements/OpalMarker.h"
#include "AbstractObjects/OpalData.h"
#include "Utilities/Options.h"
#include "BeamlineCore/MarkerRep.h"

#include <iostream>

using namespace std;

// Class OpalMarker
// ------------------------------------------------------------------------

OpalMarker::OpalMarker():
  OpalElement(COMMON, "MARKER",
	     "The \"MARKER\" element defines a marker.")
{
  setElement(new MarkerRep("MARKER"));

  // Construct the begin marker for beam lines.
  OpalMarker *S = new OpalMarker("#S", this);
  S->getElement()->makeSharable();
  S->builtin = true;
  OPAL.create(S);

  // Construct the end marker for beam lines.
  OpalMarker *E = new OpalMarker("#E", this);
  E->getElement()->makeSharable();
  E->builtin = true;
  OPAL.create(E);
}


OpalMarker::OpalMarker(const string &name, OpalMarker *parent):
  OpalElement(name, parent)
{
  setElement(new MarkerRep(name));
}


OpalMarker::~OpalMarker()
{}


OpalMarker *OpalMarker::clone(const string &name)
{
  return new OpalMarker(name, this);
}


void OpalMarker::print(std::ostream &os) const
{
  if (Options::opal8) {
    string head = getOpalName();
    Object *parent = getParent();
    if (parent != 0  &&  ! parent->getOpalName().empty()) {
      if (! getOpalName().empty()) head += ':';
      head += parent->getOpalName();
      parent->print();
    }

    //    os << head;
    int len = head.length();

    // Only TYPE attribute permitted in OPAL-8.
    itsAttr[TYPE].print(len);

    //JMJ 30/11/2000 Semi-colons are OK in OPAL8 so suppress following conditional
    // JMJ    if (! Options::opal8)
    //    os << ';';
    //os << std::endl;
  } else {
    OpalElement::print(os);
  }
}


void OpalMarker::update()
{
  // Transmit "unknown" attributes.
  MarkerRep *mark = dynamic_cast<MarkerRep *>(getElement());
  OpalElement::updateUnknown(mark);
}
