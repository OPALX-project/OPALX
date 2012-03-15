#ifndef CLASSIC_Monitor_HH
#define CLASSIC_Monitor_HH

// ------------------------------------------------------------------------
// $RCSfile: Monitor.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Monitor
//   Defines the abstract interface for a beam position monitor.
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
#include "AbsBeamline/BeamlineVisitor.h"
#include "BeamlineGeometry/StraightGeometry.h"

#include <hdf5.h>
#include "H5Part.h"
#include <list>

// Class Monitor
// ------------------------------------------------------------------------
/// Interface for beam position monitors.
//  Class Monitor defines the abstract interface for general beam position
//  monitors.

class Monitor: public Component {

public:

    /// Plane selection.
    enum Plane {
        /// Monitor is off (inactive).
        OFF,
        /// Monitor acts on x-plane.
        X,
        /// Monitor acts on y-plane.
        Y,
        /// Monitor acts on both planes.
        XY
    };

    /// Constructor with given name.
    explicit Monitor(const string &name);

    Monitor();
    Monitor(const Monitor &);
    virtual ~Monitor();

    /// Apply visitor to Monitor.
    virtual void accept(BeamlineVisitor &) const;

    /// Get geometry.
    virtual StraightGeometry &getGeometry() = 0;

    /// Get geometry. Version for const object.
    virtual const StraightGeometry &getGeometry() const = 0;

    /// Get plane on which monitor observes.
    virtual Plane getPlane() const = 0;

    virtual bool apply(const int &i, const double &t, double E[], double B[]);

    virtual bool apply(const int &i, const double &t, Vector_t &E, Vector_t &B);

    virtual bool apply(const Vector_t &R, const Vector_t &centroid, const double &t, Vector_t &E, Vector_t &B);

    virtual void initialise(PartBunch *bunch, double &startField, double &endField, const double &scaleFactor);

    virtual void finalise();

    virtual bool bends() const;

    virtual void goOnline();

    virtual void goOffline();

    virtual const string &getType() const;

    virtual void getDimensions(double &zBegin, double &zEnd) const;

    void setOutputFN(string fn);

    void moveBy(const double & dz);
private:

    // Not implemented.
    void operator=(const Monitor &);
    string filename_m;               /**< The name of the outputfile*/
    Plane plane_m;
    double position_m;
    list<double> PosX_m;
    list<double> PosY_m;
    list<double> MomentumX_m;
    list<double> MomentumY_m;
    list<double> MomentumZ_m;
    list<double> time_m;
    list<int> id_m;
    bool informed_m;
    unsigned int step_m;

    static map<string, unsigned int> h5pfiles_s;
};

#endif // CLASSIC_Monitor_HH
