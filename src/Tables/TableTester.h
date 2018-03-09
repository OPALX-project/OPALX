// ------------------------------------------------------------------------
// $RCSfile: TableTester.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: TableTester
//
// ------------------------------------------------------------------------
//
// $Date: 2001/08/13 15:25:22 $
// $Author: jowett $
//
// ------------------------------------------------------------------------

#include "Algorithms/DefaultVisitor.h"
#include <string>

class ElementBase;


// Class TableTester
// ------------------------------------------------------------------------
/// Test dependency of a table.
//  Visitor for testing for dependency of a table on a name.
//  It walks through a beam line and test for occurrence of a given
//  element name.

class TableTester: public DefaultVisitor {

public:

    /// Constructor.
    //  Attach visitor to [b]bl[/b], remember the name [b]name[/b].
    TableTester(const Beamline &bl, const std::string &name);

    ~TableTester();

    /// Apply default operation.
    //  Throw an exception, if the contained element has the given name.
    virtual void applyDefault(const ElementBase &);

private:

    // Not implemented.
    TableTester();
    TableTester(const TableTester &);
    void operator=(const TableTester &);

    // The name to be tested.
    const std::string itsName;
};
