#ifndef OPAL_Line_HH
#define OPAL_Line_HH

// ------------------------------------------------------------------------
// $RCSfile Line.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Line
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:42 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/BeamSequence.h"
#include "AbsBeamline/ElementBase.h"
#include "Beamlines/FlaggedBeamline.h"
#include "MemoryManagement/Pointer.h"

class TokenStream;


// Class Line
// ------------------------------------------------------------------------
/// The LINE definition.
//  A Line contains a CLASSIC TBeamline<FlaggedElmPtr> which represents the
//  sequence of elements in the line.  The line is always flat in the sense
//  that nested anonymous lines are flattened.

class Line: public BeamSequence {

    friend class LineTemplate;

public:

    /// Exemplar constructor.
    Line();

    virtual ~Line();

    /// Make clone.
    //  The new object is an empty line, it will be filled by the parser.
    virtual Line *clone(const std::string &name);

    /// Make complete copy.
    //  Copy also the line list.
    virtual Line *copy(const std::string &name);

    /// Return line length.
    virtual double getLength() const;

    /// Make a line template.
    //  The template gets the name [b]name[/b], [b]is[/b] is ignored,
    //  and the formals and the line list are read from [b]stat[/b].
    virtual Object *makeTemplate
    (const std::string &name, TokenStream &is, Statement &stat);

    /// Parse the line object.
    //  Read the definition from [b]stat[/b].
    virtual void parse(Statement &stat);

    /// Print the line.
    virtual void print(std::ostream &stream) const;

private:

    // Not implemented.
    Line(const Line &);
    void operator=(const Line &);

    // Clone constructor.
    Line(const std::string &name, Line *parent);

    // Return the embedded CLASSIC beam line.
    // The result it the ideal line.
    virtual FlaggedBeamline *fetchLine() const;

    // Parse sub-list.
    void parseList(Statement &);

    // Replace references to elements.
    virtual void replace(Object *oldObject, Object *newObject);
};

#endif // OPAL_Line_HH
