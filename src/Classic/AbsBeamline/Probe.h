#ifndef CLASSIC_Probe_HH
#define CLASSIC_Probe_HH

// ------------------------------------------------------------------------
// $RCSfile: Probe.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Probe
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2009/10/07 09:32:32 $
// $Author: Bi, Yang $
// 2012/03/01: fix bugs and change the algorithm in the checkProbe()
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Component.h"
#include <string>
#include <utility>

template <class T, unsigned Dim>
class PartBunchBase;

class LossDataSink;
class PeakFinder;

// Class Probe
// ------------------------------------------------------------------------
/// Interface for probe.
//  Class Probe defines the abstract interface for a probe.

class Probe: public Component {

public:

    /// Constructor with given name.
    explicit Probe(const std::string &name);

    Probe();
    Probe(const Probe &);
    virtual ~Probe();

    /// Apply visitor to Probe.
    virtual void accept(BeamlineVisitor &) const;

    virtual void initialise(PartBunchBase<double, 3> *bunch, double &startField, double &endField);
    virtual void initialise(PartBunchBase<double, 3> *bunch);

    virtual void finalise();

    virtual bool bends() const;

    virtual void goOffline();

    /// Set dimensions and consistency checks
    void setDimensions(double xstart, double xend, double ystart, double yend, double width);

    void setStep(double step);

    virtual double getXstart() const;
    virtual double getXend() const;

    virtual double getYstart() const;
    virtual double getYend() const;

    virtual double getWidth() const;
    virtual double getStep() const;

    bool  checkProbe(PartBunchBase<double, 3> *bunch, const int turnnumber, const double t, const double tstep);
    virtual ElementBase::ElementType getType() const;

    virtual void getDimensions(double &zBegin, double &zEnd) const;

private:
    std::string filename_m;             /**< The name of the inputfile*/
    double position_m;
    ///@{ input geometry positions
    double xstart_m;
    double xend_m;
    double ystart_m;
    double yend_m;
    double rstart_m;
    double rend_m;
    ///@}
    double width_m;   ///< bin width, not used
    Point  geom_m[5]; ///< actual geometry positions with adaptive width such that each particle hits probe once per turn
    int step_m; ///< Step size of the probe (bin width in histogram file)

    double A_m, B_m,R_m, C_m;
    void setGeom(const double dist);
    int  checkPoint( const double & x, const double & y );

    std::unique_ptr<PeakFinder> peakfinder_m;

    std::unique_ptr<LossDataSink> lossDs_m;

    // Not implemented.
    void operator=(const Probe &);
};

#endif // CLASSIC_Probe_HH
