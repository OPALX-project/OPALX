#ifndef CLASSIC_Collimator_HH
#define CLASSIC_Collimator_HH

// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Collimator
//   Defines the abstract interface for a beam Collimator.
//   *** MISSING *** Collimator interface is still incomplete.
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

//#include "H5hut.h"
#include <vector>

typedef struct h5_file h5_file_t;
class LossDataSink;

// Class Collimator
// ------------------------------------------------------------------------
/// Abstract collimator.
//  Class Collimator defines the abstract interface for a collimator.

class Collimator: public Component {

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
    explicit Collimator(const string &name);

    Collimator();
    Collimator(const Collimator &rhs);
    virtual ~Collimator();

    /// Apply visitor to Collimator.
    virtual void accept(BeamlineVisitor &) const;

    /// Return the horizontal half-aperture.
    virtual double getXsize() const {return a_m;}

    /// Return the vertical half-aperture.
    virtual double getYsize() const {return b_m;}

    virtual bool apply(const int &i, const double &t, double E[], double B[]);

    virtual bool apply(const int &i, const double &t, Vector_t &E, Vector_t &B);

    virtual bool apply(const Vector_t &R, const Vector_t &centroid, const double &t, Vector_t &E, Vector_t &B);

    virtual bool checkCollimator(PartBunch &bunch, const int turnnumber, const double t, const double tstep);

    virtual void initialise(PartBunch *bunch, double &startField, double &endField, const double &scaleFactor);

    virtual void initialise(PartBunch *bunch, const double &scaleFactor);

    virtual void finalise();

    virtual bool bends() const;

    virtual void goOnline();

    virtual void goOffline();

    virtual const string &getType() const;

    virtual void getDimensions(double &zBegin, double &zEnd) const;

    string  getCollimatorShape();
    void setOutputFN(string fn);
    string getOutputFN();

    void setXsize(double a) ;

    void setYsize(double b) ;

    void setXpos(double x0) ;

    void setYpos(double y0) ;

    double getXsize(double a) ;

    double getYsize(double b) ;

    double getXpos() ;

    double getYpos() ;

    // --------Cyclotron collimator

    void setXStart(double xstart) ;

    void setYStart(double ystart) ;

    void setXEnd(double xend) ;

    void setYEnd(double yend) ;

    void setWidth(double width) ;

    double getXStart() ;

    double getYStart() ;

    double getXEnd() ;

    double getYEnd() ;

    double getWidth() ;

    //-----------------

    void setRHole(double r) ;

    void setNHoles(unsigned int nx, unsigned int ny) ;

    void setPitch(double p) ;

    void setPepperPot() ;
    void setSlit() ;
    void setRColl() ;
    void setCColl() ;
    void setWire() ;

private:

    // Not implemented.
    void operator=(const Collimator &);
    h5_file_t *H5file_m;
    string filename_m;               /**< The name of the outputfile*/
    Plane plane_m;
    double position_m;
    std::vector<double> PosX_m;
    std::vector<double> PosY_m;
    std::vector<double> PosZ_m;
    std::vector<double> MomentumX_m;
    std::vector<double> MomentumY_m;
    std::vector<double> MomentumZ_m;
    std::vector<double> time_m;
    std::vector<int> id_m;
    bool informed_m;
    double a_m;
    double b_m;
    double x0_m;
    double y0_m;

    //parameters for CCollimator
    double xstart_m;
    double xend_m;
    double ystart_m;
    double yend_m;
    double width_m;


    /** This defines a pepperpot */
    bool isAPepperPot_m;
    bool isASlit_m;
    bool isARColl_m;
    bool isACColl_m;
    bool isAWire_m;
    double rHole_m;
    unsigned int nHolesX_m;
    unsigned int nHolesY_m;
    double pitch_m;

    Point  geom_m[5];
    void setGeom();
    int  checkPoint( const double & x, const double & y );
    LossDataSink *lossDs_m;
};

#endif // CLASSIC_Collimator_HH
