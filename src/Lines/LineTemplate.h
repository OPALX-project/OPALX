#ifndef OPAL_LineTemplate_HH
#define OPAL_LineTemplate_HH

// ------------------------------------------------------------------------
// $RCSfile: LineTemplate.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: LineTemplate
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/29 10:41:40 $
// $Author: opal $
//
// ------------------------------------------------------------------------

#include "OpalParser/Macro.h"
#include "OpalParser/MacroStream.h"

class Line;
class Statement;
class TokenStream;


// Class LineTemplate
// ------------------------------------------------------------------------
/// An ``archetype'' for a OPAL beam line with arguments.
//  The model is stored in form of a MacroStream.  A call to the macro line
//  is expanded by first replacing the arguments, and then parsing the
//  resulting stream as a LINE definition.

class LineTemplate: public Macro {

    friend class Line;

public:

    LineTemplate();
    virtual ~LineTemplate();

    /// Make clone.
    //  Throw OpalException, since the template cannot be cloned.
    virtual LineTemplate *clone(const string &name);

    /// Make line instance.
    //  The instance gets the name [b]name[/b], and its actual arguments
    //  are read from [b]stat[/b]. The parser is ignored.
    virtual Object *makeInstance
    (const string &name, Statement &stat, const Parser *);

    /// Make a line template.
    //  Return NULL, since one cannot make a template from a template.
    virtual Object *makeTemplate(const string &, TokenStream &, Statement &);

    /// Parse the line template.
    //  Read the actual arguments from [b]stat[/b]. [b]is[/b] is not used.
    void parseTemplate(TokenStream &is, Statement &stat);

private:

    // Not implemented.
    LineTemplate(const LineTemplate &);
    void operator=(const LineTemplate &);

    // Clone constructor.
    LineTemplate(const string &name, Object *parent);

    // The contained beam line element list.
    MacroStream body;
};

#endif // OPAL_LineTemplate_HH
