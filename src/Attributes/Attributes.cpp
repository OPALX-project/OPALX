// ------------------------------------------------------------------------
// $RCSfile: Attributes.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Namespace Attributes
//
// ------------------------------------------------------------------------
//
// $Date: 2002/01/17 22:18:36 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Attribute.h"
#include "Attributes/Bool.h"
#include "Attributes/BoolArray.h"
#include "Attributes/Place.h"
#include "Attributes/Range.h"
#include "Attributes/Real.h"
#include "Attributes/RealArray.h"
#include "Attributes/Reference.h"
#include "Attributes/opalstr.h"
#include "Attributes/StringArray.h"
#include "Attributes/TableRow.h"
#include "Attributes/TokenList.h"
#include "Attributes/TokenListArray.h"
#include "Expressions/AValue.h"
#include "Expressions/SRefAttr.h"
#include "Expressions/SValue.h"
#include "Utilities/OpalException.h"

using std::string;
using namespace Expressions;


// Namespace Attributes.
// ------------------------------------------------------------------------

namespace Attributes {

    // ----------------------------------------------------------------------
    // Boolean value.

    Attribute makeBool(const std::string &name, const std::string &help) {
        return Attribute(new Bool(name, help), 0);
    }


    Attribute makeBool(const std::string &name, const std::string &help, bool ini) {
        return Attribute(new Bool(name, help), new SValue<bool>(ini));
    }


    bool getBool(const Attribute &attr) {
        if(AttributeBase *base = &attr.getBase()) {
            if(dynamic_cast<Bool *>(&attr.getHandler())) {
                return dynamic_cast<SValue<bool> *>(base)->evaluate();
            } else if(SValue<SRefAttr<bool> > *ref =
                          dynamic_cast<SValue<SRefAttr<bool> > *>(base)) {
                const SRefAttr<bool> &value = ref->evaluate();
                return value.evaluate();
            } else {
                throw OpalException("Attributes::get()", "Attribute \"" +
                                    attr.getName() + "\" is not logical.");
            }
        } else {
            return false;
        }
    }


    void setBool(Attribute &attr, bool val) {
        if(dynamic_cast<const Bool *>(&attr.getHandler())) {
            attr.set(new SValue<bool>(val));
        } else if(SValue<SRefAttr<bool> > *ref =
                      dynamic_cast<SValue<SRefAttr<bool> >*>(&attr.getBase())) {
            const SRefAttr<bool> &value = ref->evaluate();
            value.set(val);
        } else {
            throw OpalException("Attributes::setBool()", "Attribute \"" +
                                attr.getName() + "\" is not logical.");
        }
    }


    // ----------------------------------------------------------------------
    // Boolean array value.
    Attribute makeBoolArray(const std::string &name, const std::string &help) {
        return Attribute(new BoolArray(name, help), 0);
    }


    std::vector<bool> getBoolArray(const Attribute &attr) {
        if(AttributeBase *base = &attr.getBase()) {
            if(AValue<bool> *value =
                   dynamic_cast<AValue<bool>*>(base)) {
                return value->evaluate();
            } else {
                throw OpalException("Attributes::getBoolArray()", "Attribute \"" +
                                    attr.getName() + "\" is not a logical array.");
            }
        } else {
            return std::vector<bool>();
        }
    }


    void setBoolArray(Attribute &attr, const std::vector<bool> &value) {
        if(dynamic_cast<const BoolArray *>(&attr.getHandler())) {
            // Use ADeferred here, since a component may be overridden later
            // by an expression.
            attr.set(new ADeferred<bool>(value));
        } else {
            throw OpalException("Attributes::setBoolArray()", "Attribute \"" +
                                attr.getName() + "\" is not a logical array");
        }
    }


    // ----------------------------------------------------------------------
    // Place value.

    Attribute makePlace(const std::string &name, const std::string &help) {
        return Attribute(new Place(name, help),
                         new SValue<PlaceRep>(PlaceRep("#S")));
    }


    PlaceRep getPlace(const Attribute &attr) {
        if(&attr.getBase()) {
            if(SValue<PlaceRep> *place =
                   dynamic_cast<SValue<PlaceRep> *>(&attr.getBase())) {
                return place->evaluate();
            } else {
                throw OpalException("Attributes::getPlace()", "Attribute \"" +
                                    attr.getName() + "\" is not a place reference.");
            }
        } else {
            return PlaceRep();
        }
    }


    void setPlace(Attribute &attr, const PlaceRep &rep) {
        if(dynamic_cast<const Place *>(&attr.getHandler())) {
            attr.set(new SValue<PlaceRep>(rep));
        } else {
            throw OpalException("Attributes::getPlace()", "Attribute \"" +
                                attr.getName() + "\" is not a place reference.");
        }
    }


    // ----------------------------------------------------------------------
    // Range value.

    Attribute makeRange(const std::string &name, const std::string &help) {
        return Attribute(new Range(name, help),
                         new SValue<RangeRep>(RangeRep()));
    }


    RangeRep getRange(const Attribute &attr) {
        if(&attr.getBase()) {
            if(SValue<RangeRep> *range =
                   dynamic_cast<SValue<RangeRep> *>(&attr.getBase())) {
                return range->evaluate();
            } else {
                throw OpalException("Attributes::get()", "Attribute \"" +
                                    attr.getName() + "\" is not a range reference.");
            }
        } else {
            return RangeRep();
        }
    }


    void setRange(Attribute &attr, const RangeRep &rep) {
        if(dynamic_cast<const Range *>(&attr.getHandler())) {
            attr.set(new SValue<RangeRep>(rep));
        } else {
            throw OpalException("Attributes::get()", "Attribute \"" +
                                attr.getName() + "\" is not a range reference.");
        }
    }


    // ----------------------------------------------------------------------
    // Real value.

    Attribute makeReal(const std::string &name, const std::string &help) {
        return Attribute(new Real(name, help), 0);
    }


    Attribute
    makeReal(const std::string &name, const std::string &help, double initial) {
        return Attribute(new Real(name, help),
                         new SValue<double>(initial));
    }


    double getReal(const Attribute &attr) {
        if(AttributeBase *base = &attr.getBase()) {
            if(dynamic_cast<Real *>(&attr.getHandler())) {
                return dynamic_cast<SValue<double> *>(base)->evaluate();
            } else if(SValue<SRefAttr<double> > *ref =
                          dynamic_cast<SValue<SRefAttr<double> > *>(base)) {
                const SRefAttr<double> &value = ref->evaluate();
                return value.evaluate();
            } else {
                throw OpalException("Attributes::getReal()", "Attribute \"" +
                                    attr.getName() + "\" is not real.");
            }
        } else {
            return 0.0;
        }
    }


    void setReal(Attribute &attr, double val) {
        if(dynamic_cast<const Real *>(&attr.getHandler())) {
            attr.set(new SValue<double>(val));
        } else if(SValue<SRefAttr<double> > *ref =
                      dynamic_cast<SValue<SRefAttr<double> >*>(&attr.getBase())) {
            const SRefAttr<double> &value = ref->evaluate();
            value.set(val);
        } else {
            throw OpalException("Attributes::setReal()", "Attribute \"" +
                                attr.getName() + "\" is not logical.");
        }
    }


    // ----------------------------------------------------------------------
    // Real array value.

    Attribute makeRealArray(const std::string &name, const std::string &help) {
        return Attribute(new RealArray(name, help), 0);
    }


    std::vector<double> getRealArray(const Attribute &attr) {
        if(AttributeBase *base = &attr.getBase()) {
            if(dynamic_cast<RealArray *>(&attr.getHandler())) {
                return dynamic_cast<AValue<double>*>(base)->evaluate();
            } else {
                throw OpalException("Attributes::getRealArray()", "Attribute \"" +
                                    attr.getName() + "\" is not a real array.");
            }
        } else {
            return std::vector<double>();
        }
    }


    void setRealArray(Attribute &attr, const std::vector<double> &value) {
        if(dynamic_cast<const RealArray *>(&attr.getHandler())) {
            // Use ADeferred here, since a component may be overridden later
            // by an expression.
            attr.set(new ADeferred<double>(value));
        } else {
            throw OpalException("Attributes::setRealArray()", "Attribute \"" +
                                attr.getName() + "\" is not a real array.");
        }
    }


    // ----------------------------------------------------------------------
    // Reference value.

    Attribute makeReference(const std::string &name, const std::string &help) {
        return Attribute(new Reference(name, help), 0);
    }


    // ----------------------------------------------------------------------
    // String value.

    Attribute makeString(const std::string &name, const std::string &help) {
        return Attribute(new String(name, help), 0);
    }


    Attribute
    makeString(const std::string &name, const std::string &help, const std::string &initial) {
        return Attribute(new String(name, help), new SValue<std::string>(initial));
    }


    std::string getString(const Attribute &attr) {
        if(AttributeBase *base = &attr.getBase()) {
            if(dynamic_cast<String *>(&attr.getHandler())) {
                return dynamic_cast<SValue<std::string> *>(base)->evaluate();
            } else if(SValue<SRefAttr<std::string> > *ref =
                          dynamic_cast<SValue<SRefAttr<std::string> > *>(base)) {
                const SRefAttr<std::string> &value = ref->evaluate();
                return value.evaluate();
            } else {
                throw OpalException("Attributes::getString()", "Attribute \"" +
                                    attr.getName() + "\" is not string.");
            }
        } else {
            return std::string();
        }
    }


    void setString(Attribute &attr, const std::string &val) {
        if(dynamic_cast<const String *>(&attr.getHandler())) {
            attr.set(new SValue<std::string>(val));
        } else if(SValue<SRefAttr<std::string> > *ref =
                      dynamic_cast<SValue<SRefAttr<std::string> >*>(&attr.getBase())) {
            const SRefAttr<std::string> &value = ref->evaluate();
            value.set(val);
        } else {
            throw OpalException("Attributes::setString()", "Attribute \"" +
                                attr.getName() + "\" is not a string.");
        }
    }


    // ----------------------------------------------------------------------
    // String array value.

    Attribute makeStringArray(const std::string &name, const std::string &help) {
        return Attribute(new StringArray(name, help), 0);
    }


    std::vector<std::string> getStringArray(const Attribute &attr) {
        if(AttributeBase *base = &attr.getBase()) {
            if(dynamic_cast<StringArray *>(&attr.getHandler())) {
                return dynamic_cast<AValue<std::string>*>(base)->evaluate();
            } else {
                throw OpalException("Attributes::getStringArray()", "Attribute \"" +
                                    attr.getName() + "\" is not a string array.");
            }
        } else {
            return std::vector<std::string>();
        }
    }


    void setStringArray(Attribute &attr, const std::vector<std::string> &value) {
        if(dynamic_cast<const StringArray *>(&attr.getHandler())) {
            // Strings are never expressions, so AValue will do here.
            attr.set(new AValue<std::string>(value));
        } else {
            throw OpalException("Attributes::setStringArray()", "Attribute \"" +
                                attr.getName() + "\" is not a string array.");
        }
    }


    // ----------------------------------------------------------------------
    // Table row reference value.

    Attribute makeTableRow(const std::string &name, const std::string &help) {
        return Attribute(new TableRow(name, help), 0);
    }


    TableRowRep getTableRow(const Attribute &attr) {
        if(&attr.getBase()) {
            if(SValue<TableRowRep> *row =
                   dynamic_cast<SValue<TableRowRep> *>(&attr.getBase())) {
                return row->evaluate();
            } else {
                throw OpalException("Attributes::get()", "Attribute \"" +
                                    attr.getName() +
                                    "\" is not a table row reference.");
            }
        } else {
            return TableRowRep();
        }
    }


    void setTableRow(Attribute &attr, const TableRowRep &rep) {
        if(dynamic_cast<const TableRow *>(&attr.getHandler())) {
            attr.set(new SValue<TableRowRep>(rep));
        } else {
            throw OpalException("Attributes::get()", "Attribute \"" +
                                attr.getName() +
                                "\" is not a table row reference.");
        }
    }


    // ----------------------------------------------------------------------
    // Token list value.

    Attribute makeTokenList(const std::string &name, const std::string &help) {
        return Attribute(new TokenList(name, help), 0);
    }


    std::list<Token> getTokenList(const Attribute &attr) {
        if(AttributeBase *base = &attr.getBase()) {
            if(dynamic_cast<TokenList *>(&attr.getHandler())) {
                return dynamic_cast<SValue<std::list<Token> > *>(base)->evaluate();
            } else {
                throw OpalException("Attributes::getTokenList()", "Attribute \"" +
                                    attr.getName() + "\" is not a token list.");
            }
        } else {
            return std::list<Token>();
        }
    }


    void setTokenList(Attribute &attr, const std::list<Token> &val) {
        if(dynamic_cast<const TokenList *>(&attr.getHandler())) {
            attr.set(new SValue<std::list<Token> >(val));
        } else {
            throw OpalException("Attributes::set()", "Attribute \"" + attr.getName() +
                                "\" is not a token list.");
        }
    }


    // ----------------------------------------------------------------------
    // Token list array value.

    Attribute makeTokenListArray(const std::string &name, const std::string &help) {
        return Attribute(new TokenListArray(name, help), 0);
    }


    std::vector<std::list<Token> > getTokenListArray(const Attribute &attr) {
        if(AttributeBase *base = &attr.getBase()) {
            if(dynamic_cast<TokenListArray *>(&attr.getHandler())) {
                return dynamic_cast<AValue<std::list<Token> > *>(base)->evaluate();
            } else {
                throw OpalException("Attributes::getTokenListArray()", "Attribute \"" +
                                    attr.getName() + "\" is not a token list array.");
            }
        } else {
            return std::vector<std::list<Token> >();
        }
    }


    void
    setTokenListArray(Attribute &attr,
                      const std::vector<std::list<Token> > &value) {
        // Token lists are never expressions, so AValue will do here.
        attr.set(new AValue<std::list<Token> >(value));
    }

}
