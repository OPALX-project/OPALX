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

    virtual void initialise(PartBunchBase<double, 3> *bunch, double &startField, double &endField);
    virtual void initialise(PartBunchBase<double, 3> *bunch);

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
    bool  checkSeptum(PartBunchBase<double, 3> &bunch);
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