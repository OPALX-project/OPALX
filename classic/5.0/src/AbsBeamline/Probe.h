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
// $Author: bi $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Component.h"


// Class Probe
// ------------------------------------------------------------------------
/// Interface for septum magnet.
//  Class Probe defines the abstract interface for a septum magnet.

class Probe: public Component {

public:

    /// Constructor with given name.
    explicit Probe(const string &name);

    Probe();
    Probe(const Probe &);
    virtual ~Probe();

    /// Apply visitor to Probe.
    virtual void accept(BeamlineVisitor &) const;

    virtual bool apply(const int &i, const double &t, double E[], double B[]);

    virtual bool apply(const int &i, const double &t, Vector_t &E, Vector_t &B);
  
    virtual bool apply(const Vector_t &R, const Vector_t &centroid, const double &t, Vector_t &E, Vector_t &B);
  
    virtual void initialise(PartBunch *bunch, double &startField, double &endField, const double &scaleFactor);
    virtual void initialise(PartBunch *bunch,const double &scaleFactor);

    virtual void finalise();

    virtual bool bends() const;
    virtual void goOnline();

    virtual void goOffline();

    void setXstart(double xstart );

    void setXend(double xend);

    void setYstart(double ystart);
    void setYend(double yend);


    void setWidth(double width);

    virtual double getXstart() const;

    virtual double getXend() const;

    virtual double getYstart() const;
    virtual double getYend() const;


    virtual double getWidth() const;
    bool  checkProbe(PartBunch &bunch,int turnnumber,double deltaTheta);
    double calculateAngle(double x, double y);

    virtual const string& getType() const;

    virtual void getDimensions(double &zBegin, double &zEnd) const;

private:
    string filename_m;             /**< The name of the inputfile*/
    double position_m;
    double xstart_m;
    double xend_m;
    double ystart_m;
    double yend_m;
    double width_m;
    vector<int> idrec_m;
    int step_m;

    LossDataSink *lossDs_m;

    // Not implemented.
    void operator=(const Probe &);
};

#endif // CLASSIC_Probe_HH
