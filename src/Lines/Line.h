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
#ifndef OPAL_Line_HH
#define OPAL_Line_HH

#include <string>
#include <vector>

#include "AbstractObjects/BeamSequence.h"
#include "Algorithms/Quaternion.hpp"
#include "OPALTypes.h"

class TokenStream;

class Line : public BeamSequence {
    friend class LineTemplate;

public:
    /// One entry of the member list: an element or another LINE, by name.
    struct Member {
        std::string name;
        bool reflected = false;  ///< Written as -NAME.
    };

    /// Exemplar constructor.
    Line();

    virtual ~Line();

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

    /// The members in line order, the #S and #E markers included; N*NAME gives N entries.
    const std::vector<Member>& getMembers() const { return members_m; }

    /// Start position of the line in the global frame (m).
    const Vector_t<double, 3>& getOrigin() const { return origin_m; }

    /// Start direction of the line in the global frame.
    const Quaternion& getDirection() const { return direction_m; }

protected:
    /// Exemplar constructor for LINE-compatible sequence statements.
    Line(const char* keyword, const char* help);

    /// Parse a LINE-compatible definition.
    void parseDefinition(Statement& stat, bool allowModifiers);

    /// Return the keyword used when printing this sequence.
    virtual const char* getSequenceKeyword() const;

    /// The members in line order, the #S and #E markers included.
    std::vector<Member> members_m;

private:
    // Not implemented.
    Line(const Line&);
    void operator=(const Line&);

    // Parse sub-list.
    void parseList(Statement&, bool allowModifiers);

    Vector_t<double, 3> origin_m{0.0, 0.0, 0.0};
    Quaternion direction_m{1.0, 0.0, 0.0, 0.0};
};

#endif  // OPAL_Line_HH
