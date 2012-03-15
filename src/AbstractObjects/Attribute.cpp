// ------------------------------------------------------------------------
// $RCSfile: Attribute.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.3.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class Attribute:
//   Interface a polymorphic attribute, including its parser.
//
// ------------------------------------------------------------------------
//
// $Date: 2002/12/09 15:06:07 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Attribute.h"
#include "AbstractObjects/AttributeBase.h"
#include "Utilities/Options.h"
#include <set>
#include <iostream>

extern Inform *gmsg;


// Class Attribute
// ------------------------------------------------------------------------

Attribute::Attribute():
  base(), handler()
{}


Attribute::Attribute(const Attribute &rhs):
  base(rhs.base), handler(rhs.handler)
{}


Attribute::Attribute(AttributeHandler *h, AttributeBase *b):
  base(b), handler(h)
{}


Attribute::~Attribute()
{}


const Attribute &Attribute::operator=(const Attribute &rhs)
{
  if (&rhs != this) {
    base = rhs.base;
    handler = rhs.handler;
  }

  return *this;
}


void Attribute::doomGet(const DoomReader &reader, int index)
{
  handler->doomGet(*this, reader, index);
}


void Attribute::doomPut(DoomWriter &writer, int index) const
{
  handler->doomPut(*this, writer, index);
}


AttributeBase &Attribute::getBase() const
{
  return *base;
}


AttributeHandler &Attribute::getHandler() const
{
  return *handler;
}


const string &Attribute::getHelp() const
{
  return handler->getHelp();
}


string Attribute::getImage() const
{
  return base->getImage();
}


const string &Attribute::getName() const
{
  return handler->getName();
}


const string &Attribute::getType() const
{
  return handler->getType();
}


bool Attribute::isDeferred() const
{
  return handler->isDeferred();
}


void Attribute::setDeferred(bool flag)
{
  handler->setDeferred(flag);
}


bool Attribute::isExpression() const
{
  return base->isExpression();
}


bool Attribute::isReadOnly() const
{
  return handler->isReadOnly();
}


void Attribute::setReadOnly(bool flag)
{
  handler->setReadOnly(flag);
}


void Attribute::parse(Statement &stat, bool eval)
{
  handler->parse(*this, stat, eval);
}


void Attribute::parseComponent(Statement &stat, bool eval, int index)
{
  handler->parseComponent(*this, stat, eval, index);
}


void Attribute::set(AttributeBase *newBase)
{
  base = newBase;
}


void Attribute::setDefault()
{
  base = handler->getDefault();
}


void Attribute::print(int &pos) const
{
  if (*this) {
    string name  = getName();
    string image = getImage();
    int step = name.length() + image.length() + 2;
    pos += step;

    if (pos > 74) {
      if (Options::opal8) {
        *gmsg << ",&" << endl << "  ";
      } else {
        *gmsg << ',' << endl << "  ";
      }
      pos = step + 2;
    } else {
      *gmsg << ',';
    }

    // JMJ 18/12/2000 adding this branch when opal8 option is on
    //     Before there was only the second version.
    //     This block and above could probably be merged to tidy up, do it one day.
    if (Options::opal8) {
       *gmsg << name << "=" << image;
    } else {
       *gmsg << name << (isExpression() ? ":=" : "=") << image;
    }
  }
}


std::ostream &operator<<(std::ostream &os, const Attribute &attr)
{
  if (attr) {
    attr.getBase().print(os);
  } else {
    *gmsg << "*undefined*";
  }

  return os;
}
