// ------------------------------------------------------------------------
// $RCSfile: TokenListArray.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class TokenListArray:
//   A class used to parse a "token list array" attribute.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:36 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Attributes/TokenListArray.h"
#include "Attributes/Attributes.h"
#include "AbstractObjects/Expressions.h"
#include "Expressions/AValue.h"
#include "Utilities/OpalException.h"
#include <vector>

using namespace Expressions;


// Class TokenListArray
// ------------------------------------------------------------------------

namespace Attributes {

  TokenListArray::TokenListArray(const string &name, const string &help):
    AttributeHandler(name, help, 0)
    {}


  TokenListArray::~TokenListArray()
    {}
 

  void TokenListArray::doomGet
    (Attribute &attr, const DoomReader &, int index) const
    {
      // Cannot read token list from DOOM data base.
    }


  void TokenListArray::doomPut
    (const Attribute &attr, DoomWriter &, int index) const
    {
      // Cannot write token list to DOOM data base.
    }


  const string &TokenListArray::getType() const
    {
      static const string type("token list array");
      return type;
    }


  void TokenListArray::parse(Attribute &attr, Statement &stat, bool) const
    {
      Attributes::setTokenListArray(attr, parseTokenListArray(stat));
    }


  void TokenListArray::parseComponent
    (Attribute &attr, Statement &statement, bool, int index) const
    {
      std::vector<std::list<Token> > array;

      if (AttributeBase *base = &attr.getBase()) {
	array = dynamic_cast<AValue<std::list<Token> > *>(base)->evaluate();
      }

      while (int(array.size()) < index) {
	array.push_back(std::list<Token>());
      }

      array[index-1] = Expressions::parseTokenList(statement);
      Attributes::setTokenListArray(attr, array);
    }

};
