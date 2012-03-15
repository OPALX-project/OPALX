// ------------------------------------------------------------------------
// $RCSfile: AttributeSet.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: AttributeSet
//   A map of name (string) versus value (double) intended to store
//   user-defined attributes.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/12/16 16:26:43 $
// $Author: mad $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/AttributeSet.h"
#include "AbsBeamline/ElementBase.h"
#include "Channels/DirectChannel.h"


// Class AttributeSet
// ------------------------------------------------------------------------

AttributeSet::AttributeSet():
  itsMap()
{}


AttributeSet::AttributeSet(const AttributeSet &rhs):
  itsMap(rhs.itsMap)
{}


AttributeSet::~AttributeSet()
{}


const AttributeSet &AttributeSet::operator=(const AttributeSet &rhs)
{
  itsMap = rhs.itsMap;
  return *this;
}


double AttributeSet::getAttribute(const string &aKey) const
{
  const_iterator index = itsMap.find(aKey);

  if (index == itsMap.end()) {
    return 0.0;
  } else {
    return index->second;
  }
}


bool AttributeSet::hasAttribute(const string &aKey) const
{
  return (itsMap.find(aKey) != itsMap.end());
}


void AttributeSet::removeAttribute(const string &aKey)
{
  itsMap.erase(aKey);
}


void AttributeSet::setAttribute(const string &aKey, double value)
{
  itsMap[aKey] = value;
}


// This method is inlined so its const version can wrap it.
// ada 3-7-2000 remove inline because KCC does not like it.

Channel *AttributeSet::getChannel(const string &aKey, bool create)
{
  NameMap::iterator index = itsMap.find(aKey);

  if (index == itsMap.end()) {
    if (create){
      itsMap[aKey] = 0.0;
      return new DirectChannel(itsMap[aKey]);
    }
    //    for (NameMap::iterator index = itsMap.begin(); index != itsMap.end(); index++)
    return NULL;
  } else {
    return new DirectChannel((*index).second);
  }
}


const ConstChannel *AttributeSet::getConstChannel(const string &aKey) const
{
  // Use const_cast to allow calling the non-const GetChannel().
  // The const return value will nevertheless inhibit set().
  return const_cast<AttributeSet*>(this)->getChannel(aKey);
}

