#ifndef OPAL_Edit_HH
#define OPAL_Edit_HH

// ------------------------------------------------------------------------
// $RCSfile: Edit.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Edit
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:38 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Editor/EditParser.h"
#include "MemoryManagement/Pointer.h"
#include "Lines/Sequence.h"
#include <string>

using std::string;

class Attribute;
class EditParser;
class ElementBase;
class PlaceRep;
class RangeRep;


// Class Edit
// ------------------------------------------------------------------------
/// This class contains all data for the sequence editor.
//  It acts as a communication area between the sequence editor commands.

class Edit {

public:

    /// Constructor.
    //  Prepares the given sequence for editing.
    //  Makes a copy with all drifts removed.
    Edit(Sequence *);

    ~Edit();

    /// The type of line contained in a sequence.
    typedef Sequence::TLine TLine;

    /// The line iterator for a sequence.
    typedef Sequence::TLine::iterator iterator;

    /// Cycle the edit sequence.
    //  The new start point is at [b]start[/b].
    bool cycle(const PlaceRep &start);

    /// Finish editing.
    //  Reconstruct the modified sequence as required.
    //  If a new name is given, make a copy, otherwise if the sequence is
    //  modified, overwrite the original.
    void finish(const string &newName);

    /// Flatten the edit sequence.
    void flatten();

    /// Install multiple elements.
    //  New element [b]elem[/b] at position [b]at[/b] from all selected
    //  elements.
    int installMultiple(ElementBase *, double);

    /// Install element relative to place.
    //  New element [b]elem[/b] at position [b]at[/b] from [b]from[/b]
    //  (if given) or from origin.
    int installSingle(const PlaceRep &, ElementBase *, double);

    /// Move multiple elements.
    //  Move all selected elements by [b]by[/b].
    int moveMultiple(double by);

    /// Move single element.
    //  Move element at [b]pos[/b] to absolute position.
    int moveSingleAbs(const PlaceRep &, double to);

    /// Move single element.
    //  Move element at [b]pos[/b] by given amount.
    int moveSingleRel(const PlaceRep &, const PlaceRep &, double to);

    /// Reflect the edit sequence.
    void reflect();

    /// Remove multiple elements.
    //  Remove all selected elements.
    int removeMultiple();

    /// Remove single element.
    //  Remove element at [b]pos[/b]
    int removeSingle(const PlaceRep &);

    /// Replace multiple elements.
    //  Replace all selected elements by [b]elem[/b].
    int replaceMultiple(ElementBase *elem);

    /// Replace single element.
    //  Replace element at [b]pos[/b] by [b]elem[/b].
    int replaceSingle(const PlaceRep &, ElementBase *elem);

    /// Select elements in the edit sequence.
    //  Use range, class and regular expression.
    int select(const RangeRep &rng, const string &cls,
               const string &typ, const string &patt);

    /// Clear all selection flags.
    void selectClear();

    /// Set all selection flags.
    void selectFull();

    /// The original sequence.
    Pointer<Sequence> itsSequence;

    /// The edit sequence.
    Pointer<TLine> itsLine;

    /// Modify flag.
    //  If true, the edit sequence is different from the original.
    bool isModified;

    /// The parser used during a sequence edit.
    EditParser parser;

    /// Pointer to the edit data.
    static Edit *block;

private:

    // Not implemented.
    Edit();
    Edit(const Edit &);
    void operator=(const Edit &);

    // Add an element to the install list.
    void install(TLine &, ElementBase *, double);

    // Install multiple elements.
    int installMultiple(bool, TLine &, ElementBase *, double);

    // Install single element.
    int installSingle(bool, TLine &, PlaceRep &, ElementBase *, double);

    // Move one element.
    void merge(TLine &, TLine &);

    // Move multiple elements.
    int moveMultiple(bool, TLine &, double by);

    // Move single element.
    int moveSingleAbs(bool, TLine &, PlaceRep &, double to);
    int moveSingleRel(bool, TLine &, PlaceRep &, PlaceRep &, double to);

    // Reflect the edit sequence.
    TLine *reflect(TLine &);

    // Remove multiple elements.
    int removeMultiple(bool, TLine &);

    // Remove single element.
    int removeSingle(bool, TLine &, PlaceRep &);

    // Replace multiple elements.
    int replaceMultiple(bool, TLine &, ElementBase *elem);

    // Replace single element.
    int replaceSingle(bool, TLine &, PlaceRep &, ElementBase *elem);

    // Warning message.
    void invalidLine(const char msg[]);
    void invalidShare(const char msg[]);
};

#endif // OPAL_Edit_HH
