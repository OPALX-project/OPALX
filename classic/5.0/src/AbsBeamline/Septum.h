#ifndef CLASSIC_Septum_HH
#define CLASSIC_Septum_HH

// ------------------------------------------------------------------------
// $RCSfile: Septum.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Septum
//   *** MISSING *** Septum interface is still incomplete.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:32 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Component.h"


// Class Septum
// ------------------------------------------------------------------------
/// Interface for septum magnet.
//  Class Septum defines the abstract interface for a septum magnet.

class Septum: public Component {

public:

    /// Constructor with given name.
    explicit Septum(const std::string &name);

    Septum();
    Septum(const Septum &);
    virtual ~Septum();

    /// Apply visitor to Septum.
    virtual void accept(BeamlineVisitor &) const;

    virtual bool apply(const size_t &i, const double &t, double E[], double B[]);

    virtual bool apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B);

    virtual bool apply(const Vector_t &R, const Vector_t &centroid, const double &t, Vector_t &E, Vector_t &B);

    virtual void initialise(PartBunch *bunch, double &startField, double &endField, const double &scaleFactor);
    virtual void initialise(PartBunch *bunch, const double &scaleFactor);

    virtual void finalise();

    virtual bool bends() const;

    virtual void goOffline();

    void setXstart(double xstart);

    void setXend(double xend);

    void setYstart(double ystart);
    void setYend(double yend);


    void setWidth(double width);

    virtual double getXstart() const;

    virtual double getXend() const;

    virtual double getYstart() const;
    virtual double getYend() const;


    virtual double getWidth() const;
    bool  checkSeptum(PartBunch &bunch);
    double calculateAngle(double x, double y);

    virtual ElementBase::ElementType getType() const;

    virtual void getDimensions(double &zBegin, double &zEnd) const;

private:
    std::string filename_m;             /**< The name of the inputfile*/
    double position_m;
    double xstart_m;
    double xend_m;
    double ystart_m;
    double yend_m;
    double width_m;
    // Not implemented.
    void operator=(const Septum &);
};

#endif // CLASSIC_Septum_HH