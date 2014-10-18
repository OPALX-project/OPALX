#ifndef CLASSIC_CyclotronValley_HH
#define CLASSIC_CyclotronValley_HH

// ------------------------------------------------------------------------
// $RCSfile: CyclotronValley.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: CyclotronValley
//   Defines the abstract interface for an element modeling the valley of 
//   a cyclotron.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2010/12/8 09:51:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------


#include "AbsBeamline/Component.h"

class Fieldmap;

// Class CyclotronValley
// ------------------------------------------------------------------------
/// Interface for cyclotron valley.
//  Class CyclotronValley defines the abstract interface for the magnetic field of cyclotron valley.


class CyclotronValley: public Component {

public:

   
    /// Constructor with given name.
    explicit CyclotronValley(const std::string &name);

    CyclotronValley();
    CyclotronValley(const CyclotronValley &);
    virtual ~CyclotronValley();

    /// Apply visitor to CyclotronValley.
    virtual void accept(BeamlineVisitor &) const;

   
    /// Set the name of the field map
    void setFieldMapFN(std::string fmapfn);

    std::string getFieldMapFN() const;

    void setFast(bool fast);

    bool getFast() const;
    const std::string &getType() const;
    virtual bool apply(const size_t &i, const double &t, double E[], double B[]);
   
    virtual bool apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B);

    virtual bool apply(const Vector_t &R, const Vector_t &centroid, const double &t, Vector_t &E, Vector_t &B);

    virtual void initialise(PartBunch *bunch, double &startField, double &endField, const double &scaleFactor);

    //virtual void initialise(PartBunch *bunch, const double &scaleFactor);

    virtual void finalise();

    virtual bool bends() const;

    virtual void goOnline();

    virtual void goOffline();

    virtual void getDimensions(double &zBegin, double &zEnd) const;

private:
    std::string filename_m;             /**< The name of the inputfile*/
    Fieldmap *fieldmap_m;
    double scale_m;              /**< scale multiplier*/
   
    double ElementEdge_m;
    double startField_m;         /**< starting point of field(m)*/
    double endField_m;
    bool fast_m;
   

    // Not implemented.
    void operator=(const CyclotronValley &);
};

#endif // CLASSIC_CyclotronValley_HH
