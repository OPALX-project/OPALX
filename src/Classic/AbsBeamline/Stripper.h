#ifndef CLASSIC_Stripper_HH
#define CLASSIC_Stripper_HH

// ------------------------------------------------------------------------
// $RCSfile: Stripper.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Stripper
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2011/07/08 11:14:04 $
// $Author: Jianjun Yang $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Component.h"

class LossDataSink;

// Class Stripper
// ------------------------------------------------------------------------
//  Class Stripper defines the abstract interface for a striping foil.

class Stripper: public Component {

public:

    /// Constructor with given name.
    explicit Stripper(const std::string &name);

    Stripper();
    Stripper(const Stripper &);
    virtual ~Stripper();

    /// Apply visitor to Stripper.
    virtual void accept(BeamlineVisitor &) const;

    virtual void initialise(PartBunchBase<double, 3> *bunch, double &startField, double &endField);
    virtual void initialise(PartBunchBase<double, 3> *bunch);

    virtual void finalise();

    virtual bool bends() const;

    virtual void goOffline();

    void setXstart(double xstart);
    virtual double getXstart() const;

    void setXend(double xend);
    virtual double getXend() const;

    void setYstart(double ystart);
    virtual double getYstart() const;

    void setYend(double yend);
    virtual double getYend() const;

    void setOPCharge(double charge);
    virtual double getOPCharge() const;

    void setOPMass(double mass);
    virtual double getOPMass() const;

    void setOPYield(double yield);
    virtual double getOPYield() const;

    void setWidth(double width);
    virtual double getWidth() const;

    void setStop(bool stopflag);
    virtual bool getStop() const;

    bool  checkStripper(PartBunchBase<double, 3> &bunch, const int turnnumber, const double t, const double tsetp);

    virtual void getDimensions(double &zBegin, double &zEnd) const;

    virtual ElementBase::ElementType getType() const;

private:
    std::string filename_m;             /**< The name of the inputfile*/
    double position_m;
    double xstart_m;
    double xend_m;
    double ystart_m;
    double yend_m;
    double width_m;
    double opcharge_m;
    double opmass_m;
    double opyield_m;
    Point  geom_m[5];
    bool   stop_m;
    std::vector<int> idrec_m;
    int step_m;

    double A_m, B_m,R_m, C_m;
    void setGeom(const double dist);
    int  checkPoint( const double & x, const double & y );

    std::unique_ptr<LossDataSink> lossDs_m;

    // Not implemented.
    void operator=(const Stripper &);
};

#endif // CLASSIC_Stripper_HH