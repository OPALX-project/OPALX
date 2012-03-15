#ifndef CLASSIC_Separator_HH
#define CLASSIC_Separator_HH

// ------------------------------------------------------------------------
// $RCSfile: Separator.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Separator
//   Defines the abstract interface for an  separator.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Component.h"


// Class Separator
// ------------------------------------------------------------------------
/// Interface for electrostatic separator.
//  Class Separator defines the abstract interface for electrostatic
//  separators.

class Separator: public Component {

public:

    /// Constructor with given name.
    explicit Separator(const string &name);

    Separator();
    Separator(const Separator &);
    virtual ~Separator();

    /// Apply visitor to Separator.
    virtual void accept(BeamlineVisitor &) const;

    /// Get horizontal component Ex of field in V/m.
    virtual double getEx() const = 0;

    /// Get vertical component Ey of field in V/m.
    virtual double getEy() const = 0;

    virtual bool apply(const int &i, const double &t, double E[], double B[]);

    virtual bool apply(const int &i, const double &t, Vector_t &E, Vector_t &B);
  
    virtual bool apply(const Vector_t &R, const Vector_t &centroid, const double &t, Vector_t &E, Vector_t &B);
  
    virtual void initialise(PartBunch *bunch, double &startField, double &endField, const double &scaleFactor);

    virtual void finalise();

    virtual bool bends() const;

    virtual const string& getType() const;

    virtual void getDimensions(double &zBegin, double &zEnd) const;

private:

    // Not implemented.
    void operator=(const Separator &);
};

#endif // CLASSIC_Separator_HH
