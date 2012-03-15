#ifndef OPAL_OpalParallelPlate_HH
#define OPAL_OpalParallelPlate_HH

// ------------------------------------------------------------------------
// $RCSfile: OpalParallelPlate.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalParallelPlate
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:39 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalElement.h"
#include "Structure/BoundaryGeometry.h"

// Class OpalParallelPlate
// ------------------------------------------------------------------------
/// The ParallelPlate element.


class OpalParallelPlate: public OpalElement {

public:

    /// The attributes of class OpalParallelPlate.
    enum {
        VOLT = COMMON,  // The peak voltage.
        GEOMETRY,       // geometry of boundary
        FREQ,           // The RF frequency.
        LAG,            // The phase lag.
       
        PLENGTH,           //distance between two plates or length in s direction 
        DX,             // Misalignment: translation in x direction
        DY,             // Misalignment: translation in y direction
        DZ,             // Misalignment: translation in z direction
     
        SIZE
    };

    /// Exemplar constructor.
    OpalParallelPlate();

    virtual ~OpalParallelPlate();

    /// Make clone.
    virtual OpalParallelPlate *clone(const string &name);

    /// Fill in all registered attributes.
    virtual void fillRegisteredAttributes(const ElementBase &, ValueFlag);

    /// Update the embedded CLASSIC cavity.
    virtual void update();

private:

    // Not implemented.
    OpalParallelPlate(const OpalParallelPlate &);
    void operator=(const OpalParallelPlate &);

    // Clone constructor.
    OpalParallelPlate(const string &name, OpalParallelPlate *parent);

   

    BoundaryGeometry *obgeo_m;


};

#endif // OPAL_OpalParallelPlate_HH
