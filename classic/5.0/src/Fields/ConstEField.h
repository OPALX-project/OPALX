#ifndef CLASSIC_ConstEField_HH
#define CLASSIC_ConstEField_HH

// ------------------------------------------------------------------------
// $RCSfile: ConstEField.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ConstEField
//
// ------------------------------------------------------------------------
// Class category: Fields
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:35 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Fields/StaticElectricField.h"


// Class ConstEField
// ------------------------------------------------------------------------
/// A homogeneous electricstatic field.

class ConstEField: public StaticElectricField {

public:

    /// Default constructor.
    //  Constructs null field.
    ConstEField();

    virtual ~ConstEField();

    /// Get component.
    //  Return the x-component of the electric field in A/m.
    virtual double getEx() const;

    /// Get component.
    //  Return the y-component of the electric field in A/m.
    virtual double getEy() const;

    /// Get component.
    //  Return the z-component of the electric field in A/m.
    virtual double getEz() const;

    /// Set component.
    //  Assign the x-component of the electric field in A/m.
    virtual void setEx(double);

    /// Set component.
    //  Assign the y-component of the electric field in A/m.
    virtual void setEy(double);

    /// Set component.
    //  Assign the z-component of the electric field in A/m.
    virtual void setEz(double);
};

#endif // CLASSIC_ConstEField_HH
