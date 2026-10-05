//
// Class Line
//   The LINE definition.
//   A Line contains an OPALX TBeamline<FlaggedElmPtr> which represents the
//   sequence of elements in the line.  The line is always flat in the sense
//   that nested anonymous lines are flattened.
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
#ifndef OPAL_Line_HH
#define OPAL_Line_HH

#include "AbsBeamline/ElementBase.h"
#include "AbstractObjects/BeamSequence.h"
#include "Beamlines/FlaggedBeamline.h"

class TokenStream;

class Line : public BeamSequence {
    friend class LineTemplate;

public:
    /// Exemplar constructor.
    Line();

    virtual ~Line();

    /// Make clone.
    //  The new object is an empty line, it will be filled by the parser.
    virtual Line* clone(const std::string& name);

    /// Make complete copy.
    //  Copy also the line list.
    virtual Line* copy(const std::string& name);

    /// Return line length.
    virtual double getLength() const;

    /// Make a line template.
    //  The template gets the name [b]name[/b], [b]is[/b] is ignored,
    //  and the formals and the line list are read from [b]stat[/b].
    virtual Object* makeTemplate(const std::string& name, TokenStream& is, Statement& stat);

    /// Parse the line object.
    //  Read the definition from [b]stat[/b].
    virtual void parse(Statement& stat);

    /// Print the line.
    virtual void print(std::ostream& stream) const;

protected:
    /// Exemplar constructor for LINE-compatible sequence statements.
    Line(const char* keyword, const char* help);

    // Clone constructor for LINE-compatible sequence statements.
    Line(const std::string& name, Line* parent);

    /// Parse a LINE-compatible definition.
    void parseDefinition(Statement& stat, bool allowModifiers);

    /// Return the keyword used when printing this sequence.
    virtual const char* getSequenceKeyword() const;

    // Return the embedded OPALX beam line.
    // The result is the ideal line.
    virtual FlaggedBeamline* fetchLine() const;

private:
    // Not implemented.
    Line(const Line&);
    void operator=(const Line&);

    // Parse sub-list.
    void parseList(Statement&, bool allowModifiers);

    // Replace references to elements.
    virtual void replace(Object* oldObject, Object* newObject);
};

#endif  // OPAL_Line_HH
