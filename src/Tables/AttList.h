#ifndef OPAL_AttList_HH
#define OPAL_AttList_HH

// ------------------------------------------------------------------------
// $RCSfile: AttList.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: AttList
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:46 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Action.h"
#include <iosfwd>

class Beamline;


// Class AttList
// ------------------------------------------------------------------------
/// The ATTLIST command.

class AttList: public Action {

public:

    /// Exemplar constructor.
    AttList();

    virtual ~AttList();

    /// Make clone.
    virtual AttList *clone(const std::string &name);

    /// Execute the command.
    virtual void execute();

private:

    // Not implemented.
    AttList(const AttList &);
    void operator=(const AttList &);

    // Clone constructor.
    AttList(const std::string &name, AttList *parent);

    // The working routine.
    void writeTable(const Beamline &line, std::ostream &os);
};

#endif // OPAL_AttList_HH
