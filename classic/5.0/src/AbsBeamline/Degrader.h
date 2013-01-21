#ifndef CLASSIC_Degrader_HH
#define CLASSIC_Degrader_HH

// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Degrader
//   Defines the abstract interface for a beam Degrader.
//   *** MISSING *** Degrader interface is still incomplete.
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
#include <vector>

class LossDataSink;

// Class Degrader
// ------------------------------------------------------------------------
/// Abstract collimator.
//  Class Degrader defines the abstract interface for a collimator.

class Degrader: public Component {

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
    explicit Degrader(const string &name);

    Degrader();
    Degrader(const Degrader &rhs);
    virtual ~Degrader();

    /// Apply visitor to Degrader.
    virtual void accept(BeamlineVisitor &) const;

    virtual bool apply(const size_t &i, const double &t, double E[], double B[]);

    virtual bool apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B);

    virtual bool apply(const Vector_t &R, const Vector_t &centroid, const double &t, Vector_t &E, Vector_t &B);

    virtual void initialise(PartBunch *bunch, double &startField, double &endField, const double &scaleFactor);

    virtual void initialise(PartBunch *bunch, const double &scaleFactor);

    virtual void finalise();

    virtual bool bends() const;

    virtual void goOnline();

    virtual void goOffline();

    virtual const string &getType() const;

    virtual void getDimensions(double &zBegin, double &zEnd) const;

    string  getDegraderShape(); // AAA
 
    void setOutputFN(string fn);
    string getOutputFN();

    void setZSize( double z) ;
   
    void setZStart(double zstart) ; 
    void setZEnd(double zend) ; 
   
    double getZStart() ; 
    double getZEnd() ; 
    double getZSize();

    virtual bool isInMaterial(double z);

private:

    // Not implemented.
    void operator=(const Degrader &);

    string filename_m;               /**< The name of the outputfile*/

    double position_m;
    double deg_width_m;

    std::vector<double> PosX_m;
    std::vector<double> PosY_m;
    std::vector<double> PosZ_m;
    std::vector<double> MomentumX_m;
    std::vector<double> MomentumY_m;
    std::vector<double> MomentumZ_m;
    std::vector<double> time_m;
    std::vector<int> id_m;
    bool informed_m;
    
    double zstart_m;
    double zend_m;

    std::unique_ptr<LossDataSink> lossDs_m;
};

#endif // CLASSIC_Degrader_HH
