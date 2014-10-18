#ifndef OPAL_Eigen_HH
#define OPAL_Eigen_HH

// ------------------------------------------------------------------------
// $RCSfile: Eigen.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Eigen
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


// Class Eigen
// ------------------------------------------------------------------------
/// The EIGEN command.

class Eigen: public Action {

public:

    /// Exemplar constructor.
    Eigen();

    virtual ~Eigen();

    /// Make clone.
    virtual Eigen *clone(const std::string &name);

    /// Execute the command.
    virtual void execute();

private:

    // Not implemented.
    Eigen(const Eigen &);
    void operator=(const Eigen &);

    // Clone constructor.
    Eigen(const std::string &name, Eigen *parent);

    // Do the listing.
    void format(std::ostream &, const Twiss *);

    /// Print Twiss table in eigenvector representation.
    void formatPrint(std::ostream &, const Twiss *) const;

};

#endif // OPAL_Eigen_HH
