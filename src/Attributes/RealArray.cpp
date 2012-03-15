// ------------------------------------------------------------------------
// $RCSfile: RealArray.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class RealArray:
//   A class used to parse real array attributes.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:36 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Attributes/RealArray.h"
#include "Attributes/Attributes.h"
#include "AbstractObjects/DoomReader.h"
#include "AbstractObjects/DoomWriter.h"
#include "AbstractObjects/Expressions.h"
#include "Expressions/AAutomatic.h"
#include "Expressions/ADeferred.h"
#include "Expressions/AList.h"
#include "Expressions/SConstant.h"
#include "Utilities/OpalException.h"
#include <vector>

using namespace Expressions;


// Class RealArray
// ------------------------------------------------------------------------

namespace Attributes {

  RealArray::RealArray(const string &name, const string &help):
    AttributeHandler(name, help, 0)
  {}


  RealArray::~RealArray()
  {}


  void RealArray::doomGet
  (Attribute &attr, const DoomReader &reader, int index) const
  {
    int size = reader.getInt(index);
    std::vector<double> array(size, 0.0);

    for (int i = 0; i < size; i++) {
      array[i] = reader.getReal(index + i + 1);
    }

    Attributes::setRealArray(attr, array);
  }


  void RealArray::doomPut
  (const Attribute &attr, DoomWriter &writer, int index) const
  {
    std::vector<double> array = Attributes::getRealArray(attr);
    int size = array.size();
    writer.putInt(index, size);

    for (int i = 0; i < size; i++) {
      writer.putReal(index + i + 1, array[i] ? 1.0 : 0.0);
    }
  }


  const string &RealArray::getType() const
  {
    static string type = "real array";
    return type;
  }


  void RealArray::parse(Attribute &attr, Statement &statement, bool eval) const
  {
    PtrToArray<double> expr = parseRealArray(statement);

    if (eval) {
      // Use ADeferred here, since a component may be overridden later
      // by an expression.
      attr.set(new ADeferred<double>(expr->evaluate()));
    } else if (is_deferred) {
      attr.set(new ADeferred<double>(expr));
    } else {
      attr.set(new AAutomatic<double>(expr));
    }
  }


  void RealArray::parseComponent
  (Attribute &attr, Statement &stat, bool eval, int index) const
  {
    ADeferred<double> *array = 0;

    if (AttributeBase *base = &(attr.getBase())) {
      array = dynamic_cast<ADeferred<double>*>(base);
    } else {
      attr.set(array = new ADeferred<double>());
    }

    PtrToScalar<double> expr = parseReal(stat);
    if (eval) expr = new SConstant<double>(expr->evaluate());
    array->setComponent(index, expr);
  }

};
