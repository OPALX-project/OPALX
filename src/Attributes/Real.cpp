// ------------------------------------------------------------------------
// $RCSfile: Real.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class Real:
//   A class used to parse real attributes.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/12/15 10:08:28 $
// $Author: opal $
//
// ------------------------------------------------------------------------

#include "Attributes/Real.h"
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


// Class Real
// ------------------------------------------------------------------------

namespace Attributes {

    Real::Real(const string &name, const string &help):
        AttributeHandler(name, help, 0)
    {}


    Real::~Real()
    {}


    void Real::doomGet(Attribute &attr, const DoomReader &reader, int index) const {
        attr.set(0);

        switch(reader.getInt(index)) {
            case 0:
                break;
            case 1:
                attr.set(new SValue<double>(reader.getReal(index)));
                break;
            case 2: {
                SimpleStatement statement(getName(), 1);

                try {
                    string exprString = reader.getString(index);
                    StringStream is(exprString);
                    Token token = is.readToken();

                    while(! token.isEOF()) {
                        statement.append(token);
                        token = is.readToken();
                    }

                    statement.start();
                    PtrToScalar<double> expr = parseReal(statement);
                    attr.set(new SAutomatic<double>(expr));
                } catch(ParseError &ex) {
                    cerr << endl << "*** Parse error detected by function \""
                         << ex.where() << "\"" << endl;
                    statement.printWhere(true);
                    cerr << "    ";
                    statement.print();
                    cerr << "    " << ex.what() << endl << endl;
                    throw;
                } catch(ClassicException &ex) {
                    cerr << endl << "*** User error detected by function \""
                         << ex.where() << "\"" << endl;
                    statement.printWhere(false);
                    cerr << "    ";
                    statement.print();
                    cerr << "    " << ex.what() << endl << endl;
                    throw;
                }
            }
            break;
        }
    }


    void Real::doomPut(const Attribute &attr, DoomWriter &writer, int index) const {
        if(attr) {
            SValue<double> &value = dynamic_cast<SValue<double> &>(attr.getBase());
            writer.putReal(index, value.evaluate());
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


    const string &Real::getType() const {
        static const string type("real");
        return type;
    }


    void Real::parse(Attribute &attr, Statement &statement, bool eval) const {
        PtrToScalar<double> expr = parseReal(statement);

        if(eval || expr->isConstant()) {
            attr.set(new SValue<double>(expr->evaluate()));
        } else if(isDeferred()) {
            attr.set(new SDeferred<double>(expr));
        } else {
            attr.set(new SAutomatic<double>(expr));
        }
    }

};
