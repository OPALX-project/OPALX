#ifndef OPAL_Selector_HH
#define OPAL_Selector_HH

// ------------------------------------------------------------------------
// $RCSfile: Selector.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Selector
//
// ------------------------------------------------------------------------
//
// $Date: 2001/08/13 15:25:22 $
// $Author: jowett $
//
// ------------------------------------------------------------------------

#include "Tables/RangeSelector.h"
#include <string>

class Element;
class RangeRep;
class RegularExpression;


// Class Selector
// ------------------------------------------------------------------------
/// Set selection flags for a given range in a beam line.

class Selector: public RangeSelector {

public:

    /// Constructor.
    //  Attach visitor to [b]bl[/b].  Remember range [b]range[/b],
    //  class name [b]cName[/b], type name [b]tName[/b], and pattern
    //  string [b]pString[/b].
    Selector(const Beamline &, const RangeRep &range, const std::string &cName,
             const std::string &tName, const std::string &pString);

    virtual ~Selector();

    /// Execute the selection.
    virtual void execute();

    /// Return the count of selected elements.
    int getCount() const;

protected:

    /// The operation to be done for elements.
    virtual void handleElement(const FlaggedElmPtr &);

private:

    // Not implemented.
    Selector();
    Selector(const Selector &);
    void operator=(const Selector &);

    // Class of "SELECT", or zero.
    const Element *itsClass;

    // Type name of "SELECT".
    const std::string itsType;

    // Pattern of "SELECT", or zero.
    const RegularExpression *itsPattern;

    // The count of selected elements.
    int itsCount;
};

#endif // OPAL_Selector_HH
