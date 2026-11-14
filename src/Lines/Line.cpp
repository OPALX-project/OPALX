//
// Class Line
//   The LINE definition.
//   A Line keeps its members by name, in order. Nested anonymous lines and
//   line macro expansions are copied in; a nested named LINE stays one member
//   and is looked up by name when the line is described.
//
// Copyright (c) 200x - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
#include "Lines/Line.h"
#include <algorithm>
#include <iostream>
#include <memory>
#include "AbstractObjects/BeamSequence.h"
#include "AbstractObjects/Expressions.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Expressions/SAutomatic.h"
#include "Expressions/SBinary.h"
#include "Expressions/SRefExpr.h"
#include "Expressions/TFunction2.h"
#include "Lines/LineTemplate.h"
#include "Lines/Ring.h"
#include "OpalParser/Statement.h"
#include "Utilities/Options.h"
#include "Utilities/ParseError.h"

using namespace Expressions;

// The attributes of class Line.
namespace {
    enum {
        TYPE,         // The type attribute.
        LENGTH,       // The line length.
        ORIGIN,       // The location of the particle source
        ORIENTATION,  // The orientation of the particle source
        X,
        Y,
        Z,
        THETA,
        PHI,
        PSI,
        SIZE
    };

    // Auxiliary function for adding two element lengths.
    double AttAdd(double a, double b) { return a + b; }
}  // namespace

Line::Line()
    : Line("LINE",
           "The \"LINE\" statement defines a beamline list.\n"
           "\t<name> : line = (<list>)") {}

Line::Line(const char* keyword, const char* help) : BeamSequence(SIZE, keyword, help) {
    itsAttr[TYPE] = Attributes::makeString("TYPE", "Design type");

    itsAttr[LENGTH] = Attributes::makeReal("L", "Total length of line in m");
    itsAttr[LENGTH].setReadOnly(true);

    itsAttr[ORIGIN] = Attributes::makeRealArray("ORIGIN", "The location of the particle source");

    itsAttr[ORIENTATION] = Attributes::makeRealArray(
            "ORIENTATION", "The Tait-Bryan angles for the orientation of the particle source");

    itsAttr[X] =
            Attributes::makeReal("X", "The x-coordinate of the location of the particle source", 0);

    itsAttr[Y] =
            Attributes::makeReal("Y", "The y-coordinate of the location of the particle source", 0);

    itsAttr[Z] =
            Attributes::makeReal("Z", "The z-coordinate of the location of the particle source", 0);

    itsAttr[THETA] = Attributes::makeReal(
            "THETA", "The rotation about the y-axis of the particle source", 0);

    itsAttr[PHI] =
            Attributes::makeReal("PHI", "The rotation about the x-axis of the particle source", 0);

    itsAttr[PSI] =
            Attributes::makeReal("PSI", "The rotation about the z-axis of the particle source", 0);

    registerOwnership(AttributeHandler::STATEMENT);
}

Line::Line(const std::string& name, Line* parent)
    :  // Do not copy list members; they will be filled at parse time.
      BeamSequence(name, parent) {}

Line::~Line() {}

Line* Line::clone(const std::string& name) { return new Line(name, this); }

Line* Line::copy(const std::string& name) {
    Line* clone      = new Line(name, this);
    clone->members_m = members_m;
    return clone;
}

double Line::getLength() const { return Attributes::getReal(itsAttr[LENGTH]); }

Object* Line::makeTemplate(const std::string& name, TokenStream& is, Statement& statement) {
    LineTemplate* macro = new LineTemplate(name, this);

    try {
        macro->parseTemplate(is, statement);
        return macro;
    } catch (...) {
        delete macro;
        macro = 0;
        throw;
    }
}

void Line::parse(Statement& stat) { parseDefinition(stat, true); }

void Line::parseDefinition(Statement& stat, bool allowModifiers) {
    static const TFunction2<double, double> plus = {"+", 4, AttAdd};

    // Check for delimiters.
    parseDelimiter(stat, '=');
    parseDelimiter(stat, '(');

    // Parse the line list.
    parseList(stat, allowModifiers);

    // Insert the begin and end markers.
    members_m.insert(members_m.begin(), Member{"#S", false});
    members_m.push_back(Member{"#E", false});

    // Construct expression for length.
    PtrToScalar<double> expr;
    for (const Member& member : members_m) {
        PtrToScalar<double> temp = new SRefExpr<double>(member.name, "L");
        if (expr.isValid() && temp.isValid()) {
            expr = SBinary<double, double>::make(plus, expr, temp);
        } else {
            expr = temp;
        }
    }
    if (expr.isValid()) itsAttr[LENGTH].set(new SAutomatic<double>(expr));

    while (stat.delimiter(',')) {
        std::string name = Expressions::parseString(stat, "Attribute name expected.");
        Attribute* attr  = findAttribute(name);

        if (attr != 0) {
            if (stat.delimiter('=')) {
                attr->parse(stat, true);
            } else if (stat.delimiter(":=")) {
                attr->parse(stat, false);
            } else {
                attr->setDefault();
            }
        }
    }

    if (itsAttr[ORIGIN] || itsAttr[ORIENTATION]) {
        std::vector<double> origin = Attributes::getRealArray(itsAttr[ORIGIN]);
        if (origin.size() == 3) {
            origin_m = Vector_t<double, 3>(origin[0], origin[1], origin[2]);
        } else {
            origin_m = Vector_t<double, 3>(0.0);
            if (itsAttr[ORIGIN]) {
                throw OpalException(
                        "Line::parse", "Parameter origin is array of 3 values (x, y, z);\n"
                                               + std::to_string(origin.size())
                                               + " values provided");
            }
        }

        std::vector<double> direction = Attributes::getRealArray(itsAttr[ORIENTATION]);
        if (direction.size() == 3) {
            const double& theta = direction[0];
            const double& phi   = direction[1];
            const double& psi   = direction[2];

            Quaternion rotTheta(cos(0.5 * theta), 0, -sin(0.5 * theta), 0);
            Quaternion rotPhi(cos(0.5 * phi), -sin(0.5 * phi), 0, 0);
            Quaternion rotPsi(cos(0.5 * psi), 0, 0, -sin(0.5 * psi));
            direction_m = rotPsi * rotPhi * rotTheta;
        } else {
            direction_m = Quaternion(1, 0, 0, 0);
            if (itsAttr[ORIENTATION]) {
                throw OpalException(
                        "Line::parse",
                        "Parameter orientation is array of 3 values (theta, phi, psi);\n"
                                + std::to_string(direction.size()) + " values provided");
            }
        }
    } else {
        const Vector_t<double, 3> origin(
                Attributes::getReal(itsAttr[X]), Attributes::getReal(itsAttr[Y]),
                Attributes::getReal(itsAttr[Z]));

        const double theta = Attributes::getReal(itsAttr[THETA]);
        const double phi   = Attributes::getReal(itsAttr[PHI]);
        const double psi   = Attributes::getReal(itsAttr[PSI]);

        Quaternion rotTheta(cos(0.5 * theta), 0, -sin(0.5 * theta), 0);
        Quaternion rotPhi(cos(0.5 * phi), -sin(0.5 * phi), 0, 0);
        Quaternion rotPsi(cos(0.5 * psi), 0, 0, -sin(0.5 * psi));

        origin_m    = origin;
        direction_m = rotPsi * rotPhi * rotTheta;
    }
}

void Line::print(std::ostream& os) const {
    os << getOpalName() << ':' << getSequenceKeyword() << "=(";
    bool seen = false;

    for (const Member& member : members_m) {
        if (member.name[0] != '#') {
            if (seen) os << ',';
            if (member.reflected) os << '-';
            os << member.name;
            seen = true;
        }
    }

    os << ')';
    os << ';';
    os << std::endl;
}

const char* Line::getSequenceKeyword() const { return "LINE"; }

void Line::parseList(Statement& stat, bool allowModifiers) {
    do {
        // Reversed member ?
        bool rev = stat.delimiter('-');
        if (rev && !allowModifiers) {
            throw ParseError("Ring::parseList()", "RING members cannot be reflected.");
        }

        // Repetition count.
        int repeat = 1;
        if (stat.integer(repeat)) {
            if (!allowModifiers) {
                throw ParseError("Ring::parseList()", "RING members cannot be repeated.");
            }
            parseDelimiter(stat, '*');
        }

        // List member.
        if (stat.delimiter('(')) {
            // Anonymous sub-line is expanded immediately. Reversing it reverses the order of
            // its members and drops their own reflection.
            Line nestedLine;
            nestedLine.parseList(stat, allowModifiers);
            const std::vector<Member>& subLine = nestedLine.members_m;

            while (repeat-- > 0) {
                if (rev) {
                    for (auto i = subLine.rbegin(); i != subLine.rend(); ++i) {
                        members_m.push_back(Member{i->name, false});
                    }
                } else {
                    members_m.insert(members_m.end(), subLine.begin(), subLine.end());
                }
            }
        } else {
            // Identifier.
            std::string name = parseString(stat, "Line member expected.");
            Object* obj      = OpalData::getInstance()->find(name);

            if (!obj && OpalData::getInstance()->isKeyword(name)) {
                throw ParseError(
                        "Line::parseList()",
                        "\"" + name + "\" is a keyword, not a defined element.");
            } else if (!obj) {
                throw ParseError("Line::parseList()", "Element \"" + name + "\" is undefined.");
            }

            if (stat.delimiter('(')) {
                // Line macro: its expansion is an anonymous LINE, whose members (its #S and #E
                // included) are copied in, in order; -NAME reflects each of them.
                std::unique_ptr<Object> instance(obj->makeInstance(name, stat, 0));
                const auto* expanded = dynamic_cast<const Line*>(instance.get());
                if (expanded == nullptr) {
                    throw ParseError(
                            "Line::parseList()",
                            "Object \"" + name + "\" cannot be a line member.");
                }
                while (repeat-- > 0) {
                    for (const Member& member : expanded->members_m) {
                        members_m.push_back(Member{member.name, member.reflected != rev});
                    }
                }
            } else if (dynamic_cast<Ring*>(obj)) {
                throw ParseError(
                        "Line::parseList()", "A RING cannot be nested in a LINE or another RING.");
            } else if (dynamic_cast<Element*>(obj)) {
                while (repeat-- > 0) {
                    members_m.push_back(Member{name, rev});
                }
            } else {
                throw ParseError(
                        "Line::parseList()", "Object \"" + name + "\" cannot be a line member.");
            }
        }
    } while (stat.delimiter(','));

    parseDelimiter(stat, ')');
}
