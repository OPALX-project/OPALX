// ------------------------------------------------------------------------
// $RCSfile: Bool.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class Attributes::Bool
//   A class used to parse boolean attributes.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:36 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Attributes/Bool.h"
#include "AbstractObjects/DoomReader.h"
#include "AbstractObjects/DoomWriter.h"
#include "AbstractObjects/Expressions.h"
#include "Expressions/SAutomatic.h"
#include "Expressions/SDeferred.h"
#include "Expressions/SValue.h"
#include "Parser/SimpleStatement.h"
#include "Parser/StringStream.h"
#include "Parser/Token.h"
#include "Utilities/ClassicException.h"
#include "Utilities/ParseError.h"

using namespace Expressions;
using std::cerr;
using std::endl;


// Class Attributes::Bool
// ------------------------------------------------------------------------

namespace Attributes {

    Bool::Bool(const string &name, const string &help):
        AttributeHandler(name, help, new SValue<bool>(true))
    {}


    Bool::~Bool()
    {}


    void Bool::
    doomGet(Attribute &attr, const DoomReader &reader, int index) const {
        attr.set(new SValue<bool>(false));

        switch(reader.getInt(index)) {
            case 0:
                break;
            case 1:
                attr.set(new SValue<bool>(reader.getReal(index) != 0.0));
                break;
            case 2: {
                SimpleStatement stat(getName(), 1);
                try {
                    string exprString = reader.getString(index);
                    StringStream is(exprString);
                    Token token = is.readToken();

                    while(! token.isEOF()) {
                        stat.append(token);
                        token = is.readToken();
                    }

                    stat.start();
                    PtrToScalar<bool> expr = parseBool(stat);
                    attr.set(new SAutomatic<bool>(expr));
                } catch(ParseError &ex) {
                    cerr << endl << "*** Parse error detected by function \""
                         << ex.where() << "\"" << endl;
                    stat.printWhere(true);
                    cerr << "    ";
                    stat.print();
                    cerr << "    " << ex.what() << endl << endl;
                    throw;
                } catch(ClassicException &ex) {
                    cerr << endl << "*** User error detected by function \""
                         << ex.where() << "\"" << endl;
                    stat.printWhere(false);
                    cerr << "    ";
                    stat.print();
                    cerr << "    " << ex.what() << endl << endl;
                    throw;
                }
            }
            break;
        }
    }


    void Bool::
    doomPut(const Attribute &attr, DoomWriter &writer, int index) const {
        if(attr) {
            SValue<bool> &value = dynamic_cast<SValue<bool> &>(attr.getBase());
            writer.putReal(index, value.evaluate() ? 1.0 : 0.0);
            if(attr.isExpression()) {
                writer.putInt(index, 2);
                writer.putString(index, attr.getImage());
            } else {
                writer.putInt(index, 1);
            }
        } else {
            writer.putInt(index, 0);
        }
    }


    const string &Bool::getType() const {
        static const string type("logical");
        return type;
    }


    void Bool::parse(Attribute &attr, Statement &stat, bool eval) const {
        PtrToScalar<bool> expr = parseBool(stat);

        if(eval || expr->isConstant()) {
            attr.set(new SValue<bool>(expr->evaluate()));
        } else if(is_deferred) {
            attr.set(new SDeferred<bool>(expr));
        } else {
            attr.set(new SAutomatic<bool>(expr));
        }
    }

};
