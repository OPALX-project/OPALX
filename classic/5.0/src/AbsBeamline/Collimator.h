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

#include <hdf5.h>
#include "H5Part.h"
#include <vector>


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
  
    virtual bool apply(const Vector_t &R, const double &t, Vector_t &E, Vector_t &B);
  
    virtual void initialise(const PartBunch *bunch, double &startField, double &endField, const double &scaleFactor);

    virtual void finalise();

    virtual void rescaleFieldMap(const double &scaleFactor);

    virtual bool bends() const;

    virtual void goOnline();

    virtual void goOffline();

    virtual const string& getType() const;

    virtual void getDimensions(double &zBegin, double &zEnd) const;

    void setOutputFN(string fn);

    void setXsize(double a) ;

    void setYsize(double b) ;
  
    double getXsize(double a) ;

    double getYsize(double b) ;

    void setRHole (double r) ;

    void setNHoles (unsigned int nx, unsigned int ny) ;
  
    void setPitch(double p) ;

    void setPepperPot() ;

private:

    // Not implemented.
    void operator=(const Collimator &);
    H5PartFile *H5file_m;
    string filename_m;               /**< The name of the outputfile*/
    Plane plane_m;
    double position_m;
    vector<double> PosX_m;
    vector<double> PosY_m;
    vector<double> PosZ_m;
    vector<double> MomentumX_m;
    vector<double> MomentumY_m;
    vector<double> MomentumZ_m;
    vector<double> time_m;
    vector<int> id_m;
    bool informed_m;
    double a_m; 
    double b_m;

    /** This defines a pepperpot */
    bool isAPepperPot_m;
    double rHole_m;
    unsigned int nHolesX_m;
    unsigned int nHolesY_m;
    double pitch_m;
};

#endif // CLASSIC_Collimator_HH
