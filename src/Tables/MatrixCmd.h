#ifndef OPAL_MatrixCmd_HH
#define OPAL_MatrixCmd_HH

// ------------------------------------------------------------------------
// $RCSfile: MatrixCmd.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: MatrixCmd
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:45 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/Action.h"
#include "AbstractObjects/Attribute.h"
#include <iosfwd>

class Twiss;


// Class MatrixCmd
// ------------------------------------------------------------------------
/// The MATRIX command.

class MatrixCmd: public Action {

public:

    /// Exemplar constructor.
    MatrixCmd();

    virtual ~MatrixCmd();

    /// Make clone.
    virtual MatrixCmd *clone(const std::string &name);

    /// Execute the command.
    virtual void execute();

private:

    // Not implemented.
    MatrixCmd(const MatrixCmd &);
    void operator=(const MatrixCmd &);

    // Clone constructor.
    MatrixCmd(const std::string &name, MatrixCmd *parent);

    // Do the listing.
    void format(std::ostream &, const Twiss *);

    /// Print Twiss table in accumulated map representation.
    void formatPrint(std::ostream &, const Twiss *) const;

};

#endif // OPAL_MatrixCmd_HH
