// ------------------------------------------------------------------------
// $RCSfile: Line.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.3.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Line
//   The class for OPAL beam lines.
//
// ------------------------------------------------------------------------
//
// $Date: 2002/12/09 15:06:08 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "Lines/Line.h"
#include "AbsBeamline/ElementBase.h"
#include "AbstractObjects/BeamSequence.h"
#include "AbstractObjects/DoomDB.h"
#include "AbstractObjects/DoomWriter.h"
#include "AbstractObjects/Expressions.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Beamlines/Beamline.h"
#include "Beamlines/FlaggedElmPtr.h"
#include "Errors/AlignHandler.h"
#include "Errors/AlignWriter.h"
#include "Errors/MPHandler.h"
#include "Errors/MPWriter.h"
#include "Expressions/SAutomatic.h"
#include "Expressions/SBinary.h"
#include "Expressions/SRefExpr.h"
#include "Expressions/TFunction2.h"
#include "Lines/LineTemplate.h"
#include "Lines/LineWriter.h"
#include "Lines/Replacer.h"
#include "Parser/Statement.h"
#include "Utilities/Options.h"
#include "Utilities/ParseError.h"
#include <algorithm>
#include <iostream>

using namespace Expressions;


// Class Line
// ------------------------------------------------------------------------

// The attributes of class Line.
namespace {
    enum {
        TYPE,        // The type attribute.
        LENGTH,      // The line length.
        SIZE
    };

    // Auxiliary function for adding two element lengths.
    double AttAdd(double a, double b)
    { return a + b; }
}


Line::Line():
    BeamSequence(SIZE, "LINE",
                 "The \"LINE\" statement defines a beamline list.\n"
                 "\t<name> : line = (<list>)") {
    itsAttr[TYPE] = Attributes::makeString
                    ("TYPE", "Design type");

    itsAttr[LENGTH] = Attributes::makeReal
                      ("L", "Total length of line in m");
    itsAttr[LENGTH].setReadOnly(true);

    setElement((new FlaggedBeamline("LINE"))->makeAlignWrapper());
}


Line::Line(const string &name, Line *parent):
    // Do not copy list members; they will be filled at parse time.
    BeamSequence(name, parent) {
    setElement((new FlaggedBeamline(name))->makeAlignWrapper());
}


Line::~Line()
{}


Line *Line::clone(const string &name) {
    return new Line(name, this);
}


Line *Line::copy(const string &name) {
    Line *clone = new Line(name, this);
    FlaggedBeamline *oldLine = fetchLine();
    FlaggedBeamline *newLine = clone->fetchLine();
    std::copy(oldLine->begin(), oldLine->end(), std::back_inserter(*newLine));
    return &*clone;
}


void Line::doomGet(const DoomReader &) {
    // Cannot be called, since lines are stored in DOOM as sequence objects.
}


void Line::doomPut(DoomWriter &writer) const {
    // Write sequence object.
    writer.setParentName("SEQUENCE");
    const FlaggedBeamline *line = fetchLine();
    Object::doomPut(writer);
    const string &lineName = getOpalName();
    LineWriter lineWriter(*line, lineName, lineName);
    lineWriter.execute();

    // Write misalignments.
    {
        AlignWriter alignWriter(lineName);
        AlignHandler extractor(*line, alignWriter, true);
        extractor.execute();
    }

    // Write field errors.
    {
        MPWriter fieldWriter(lineName);
        MPHandler extractor(*line, fieldWriter, true);
        extractor.execute();
    }
}


double Line::getLength() const {
    return Attributes::getReal(itsAttr[LENGTH]);
}


Object *Line::makeTemplate
(const string &name, TokenStream &is, Statement &statement) {
    LineTemplate *macro = new LineTemplate(name, this);

    try {
        macro->parseTemplate(is, statement);
        return macro;
    } catch(...) {
        delete macro;
        macro = 0;
        throw;
    }
}


void Line::parse(Statement &stat) {
    static const TFunction2<double, double> plus = { "+", 4, AttAdd };

    // Check for delimiters.
    parseDelimiter(stat, '=');
    parseDelimiter(stat, '(');

    // Clear all reference counts.
    OPAL.apply(OpalData::ClearReference());

    // Parse the line list.
    parseList(stat);

    // Insert the begin and end markers.
    FlaggedBeamline *line = fetchLine();
    ElementBase *markS = Element::find("#S")->getElement();
    FlaggedElmPtr first(markS, false);
    line->push_front(first);

    ElementBase *markE = Element::find("#E")->getElement();
    FlaggedElmPtr last(markE, false);
    line->push_back(last);

    // Construct expression for length, and fill in occurrence counters.
    PtrToScalar<double> expr;
    for(FlaggedBeamline::iterator i = line->begin(); i != line->end(); ++i) {
        // Accumulate length.
        PtrToScalar<double> temp =
            new SRefExpr<double>(i->getElement()->getName(), "L");
        if(expr.isValid() &&  temp.isValid()) {
            expr = SBinary<double, double>::make(plus, expr, temp);
        } else {
            expr = temp;
        }

        // Do the required instantiations.
        ElementBase *base = i->getElement();
        i->setElement(base->copyStructure());
        Element *elem = Element::find(base->getName());
        i->setCounter(elem->increment());
    }
    if(expr.isValid()) itsAttr[LENGTH].set(new SAutomatic<double>(expr));
}


void Line::print(std::ostream &os) const {
    if(! Options::opal8 && isShared()) os << "SHARED ";
    os << getOpalName() << ":LINE=(";
    const FlaggedBeamline *line = fetchLine();
    bool seen = false;

    for(FlaggedBeamline::const_iterator i = line->begin();
        i != line->end(); ++i) {
        const string name = i->getElement()->getName();
        if(name[0] != '#') {
            if(seen) os << ',';
            if(i->getReflectionFlag()) os << '-';
            os << name;
            seen = true;
        }
    }

    os << ')';
    //JMJ 30/11/2000 Comment out conditional since semi-colons are harmless
    // in OPAL8 anyway
    //if (! Options::opal8)
    os << ';';
    os << std::endl;
}


FlaggedBeamline *Line::fetchLine() const {
    return dynamic_cast<FlaggedBeamline *>(getElement()->removeWrappers());
}


void Line::parseList(Statement &stat) {
    FlaggedBeamline *line = fetchLine();

    do {
        // Reversed member ?
        bool rev = stat.delimiter('-');

        // Repetition count.
        int repeat = 1;
        if(stat.integer(repeat)) parseDelimiter(stat, '*');

        // List member.
        if(stat.delimiter('(')) {
            // Anonymous sub-line is expanded immediately.
            Line nestedLine;
            nestedLine.parseList(stat);
            FlaggedBeamline *subLine = nestedLine.fetchLine();

            while(repeat-- > 0) {
                if(rev) {
                    for(FlaggedBeamline::reverse_iterator i = subLine->rbegin();
                        i != subLine->rend(); ++i) {
                        FlaggedElmPtr fep(*i);
                        fep.setReflectionFlag(false);
                        line->push_back(fep);
                    }
                } else {
                    for(FlaggedBeamline::iterator i = subLine->begin();
                        i != subLine->end(); ++i) {
                        FlaggedElmPtr fep(*i);
                        line->push_back(fep);
                    }
                }
            }
        } else {
            // Identifier.
            string name = parseString(stat, "Line member expected.");
            Pointer<Object> obj = OPAL.find(name);

            if(! obj.isValid()) {
                // Name not found; create marker instead.
                obj = OPAL.find("MARKER")->clone(name);
                OPAL.define(&*obj);

                if(Options::verify) {
                    std::cerr << "\nElement \"" << name << "\" is undefined; "
                              << "\"MARKER\" created with this name.\n" << std::endl;
                }
            }

            if(stat.delimiter('(')) {
                // Line or sequence macro.
                // This instance will always be an anonymous object.
                obj = obj->makeInstance(name, stat, 0);
            }

            if(Element *elem = dynamic_cast<Element *>(&*obj)) {
                while(repeat-- > 0) {
                    ElementBase *base = elem->getElement();
                    FlaggedElmPtr member(base, rev);
                    line->push_back(member);
                }
            } else {
                throw ParseError("Line::parseList()",
                                 "Object \"" + name + "\" cannot be a line member.");
            }
        }
    } while(stat.delimiter(','));

    parseDelimiter(stat, ')');
}


void Line::replace(Object *oldObject, Object *newObject) {
    Element *oldElement = dynamic_cast<Element *>(oldObject);
    Element *newElement = dynamic_cast<Element *>(newObject);

    if(oldElement != 0  &&  newElement != 0) {
        Replacer rep(*fetchLine(),
                     oldElement->getOpalName(),
                     newElement->getElement());
        rep.execute();
    }
}
