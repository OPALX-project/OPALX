// ------------------------------------------------------------------------
// $RCSfile: ElementBase.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ElementBase
//   The very base class for beamline representation objects.  A beamline
//   is modelled as a composite structure having a single root object
//   (the top level beamline), which contains both "single" leaf-type
//   elements (Components), as well as sub-lines (composites).
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/12/16 16:26:43 $
// $Author: mad $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/ElementBase.h"
#include "AbsBeamline/AlignWrapper.h"
#include "AbsBeamline/ElementImage.h"
#include "BeamlineGeometry/Euclid3D.h"
#include "BeamlineGeometry/Geometry.h"
#include "Channels/Channel.h"
#include <string>
using namespace std;



// Class ElementBase
// ------------------------------------------------------------------------

ElementBase::ElementBase():
  RCObject(),
  shareFlag(true),
  elementID(""),
  userAttribs()
{}


ElementBase::ElementBase(const ElementBase &right):
  RCObject(),
  shareFlag(true),
  elementID(right.elementID),
  userAttribs(right.userAttribs)
{}


ElementBase::ElementBase(const string &name):
  RCObject(),
  shareFlag(true),
  elementID(name),
  userAttribs()
{}


ElementBase::~ElementBase()

{}


const string &ElementBase::getName() const

{
  return elementID;
}


void ElementBase::setName(const string &name)
{
  elementID = name;
}


double ElementBase::getAttribute(const string &aKey) const
{
  const ConstChannel *aChannel = getConstChannel(aKey);

  if (aChannel != NULL) {
    double val = *aChannel;
    delete aChannel;
    return val;
  } else {
    return 0.0;
  }
}


bool ElementBase::hasAttribute(const string &aKey) const
{
  const ConstChannel *aChannel = getConstChannel(aKey);

  if (aChannel != NULL) {
    delete aChannel;
    return true;
  } else {
    return false;
  }
}


void ElementBase::removeAttribute(const string &aKey)
{
  userAttribs.removeAttribute(aKey);
}


void ElementBase::setAttribute(const string &aKey, double val)
{
  Channel *aChannel = getChannel(aKey, true);

  if (aChannel != NULL  &&  aChannel->isSettable()) {
    *aChannel = val;
    delete aChannel;
  }
  else
    cout << "Channel NULL or not Settable" << endl;
}


Channel *ElementBase::getChannel(const string &aKey, bool create)
{
  return userAttribs.getChannel(aKey, create);
}


const ConstChannel *ElementBase::getConstChannel(const string &aKey) const
{
  // Use const_cast to allow calling the non-const method GetChannel().
  // The const return value of this method will nevertheless inhibit set().
  return const_cast<ElementBase*>(this)->getChannel(aKey);
}


ElementImage *ElementBase::getImage() const
{
  return new ElementImage(getName(), getType(), userAttribs);
}


ElementBase *ElementBase::copyStructure()
{
  if (isSharable()) {
    return this;
  } else {
    return clone();
  }
}


void ElementBase::makeSharable()
{
  shareFlag = true;
}


ElementBase *ElementBase::makeAlignWrapper()
{
  ElementBase *wrap = new AlignWrapper(this);
  wrap->setName(getName());
  return wrap;
}


ElementBase *ElementBase::makeFieldWrapper()
{
  return this;
}


ElementBase *ElementBase::makeWrappers()
{
  return makeFieldWrapper()->makeAlignWrapper();
}


ElementBase *ElementBase::removeAlignWrapper()
{
  return this;
}


const ElementBase *ElementBase::removeAlignWrapper() const
{
  return this;
}


ElementBase *ElementBase::removeFieldWrapper()
{
  return this;
}


const ElementBase *ElementBase::removeFieldWrapper() const
{
  return this;
}


ElementBase *ElementBase::removeWrappers()
{
  return this;
}


const ElementBase *ElementBase::removeWrappers() const
{
  return this;
}


bool ElementBase::update(const AttributeSet &set)
{
  for (AttributeSet::const_iterator i = set.begin(); i != set.end(); ++i) {
    setAttribute(i->first, i->second);
  }

  return true;
}
