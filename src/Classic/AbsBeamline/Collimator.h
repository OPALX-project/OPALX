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

#include <vector>

class BeamlineVisitor;
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
    explicit Collimator(const std::string &name);

    Collimator();
    Collimator(const Collimator &rhs);
    virtual ~Collimator();

    /// Apply visitor to Collimator.
    virtual void accept(BeamlineVisitor &) const;

    /// Return the horizontal half-aperture.
    virtual double getXsize() const {return a_m;}

    /// Return the vertical half-aperture.
    virtual double getYsize() const {return b_m;}

    virtual bool apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B);

    virtual bool applyToReferenceParticle(const Vector_t &R, const Vector_t &P, const double &t, Vector_t &E, Vector_t &B);

    virtual bool checkCollimator(PartBunchBase<double, 3> *bunch, const int turnnumber, const double t, const double tstep);

    virtual bool checkCollimator(Vector_t r, Vector_t rmin, Vector_t rmax);

    virtual void initialise(PartBunchBase<double, 3> *bunch, double &startField, double &endField);

    virtual void initialise(PartBunchBase<double, 3> *bunch);

    virtual void finalise();

    virtual bool bends() const;

    virtual void goOnline(const double &kineticEnergy);

    virtual void goOffline();

    virtual ElementBase::ElementType getType() const;

    virtual void getDimensions(double &zBegin, double &zEnd) const;

    void print();

    std::string  getCollimatorShape();
    void setOutputFN(std::string fn);
    std::string getOutputFN();

    unsigned int getLosses() const;

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
    void setZStart(double zstart) ;
    void setXEnd(double xend) ;
    void setYEnd(double yend) ;
    void setZEnd(double zend) ;
    void setWidth(double width) ;

    double getXStart() ;
    double getYStart() ;
    double getZStart() ;
    double getXEnd() ;
    double getYEnd() ;
    double getZEnd() ;
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

    virtual bool isInColl(Vector_t R, Vector_t P, double recpgamma);

    int  checkPoint(const double & x, const double & y );
private:

    // Not implemented.
    void operator=(const Collimator &);

    std::string filename_m;               /**< The name of the outputfile*/
    Plane plane_m;

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
    double zstart_m;
    double zend_m;
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

    unsigned int losses_m;

    std::unique_ptr<LossDataSink> lossDs_m;

    ParticleMaterInteractionHandler *parmatint_m;


};

inline
unsigned int Collimator::getLosses() const {
    return losses_m;
}

inline
void Collimator::setXsize(double a) {
    a_m = a;
}

inline
void Collimator::setYsize(double b) {
    b_m = b;
}

inline
void Collimator::setXpos(double x0) {
    x0_m = x0;
}

inline
void Collimator::setYpos(double y0) {
    y0_m = y0;
}

inline
double Collimator::getXsize(double a) {
    return a_m;
}

inline
double Collimator::getYsize(double b) {
    return b_m;
}

inline
double Collimator::getXpos() {
    return x0_m;
}

inline
double Collimator::getYpos() {
    return y0_m;
}

inline
void Collimator::setXStart(double xstart) {
    xstart_m = xstart;
}

inline
void Collimator::setXEnd(double xend) {
    xend_m = xend;
}

inline
void Collimator::setYStart(double ystart) {
    ystart_m = ystart;
}

inline
void Collimator::setYEnd(double yend) {
    yend_m = yend;
}

inline
void Collimator::setZStart(double zstart) {
    zstart_m = zstart;
}

inline
void Collimator::setZEnd(double zend) {
    zend_m = zend;
}

inline
void Collimator::setWidth(double width) {
    width_m = width;
}

inline
double Collimator::getXStart() {
    return xstart_m;
}

inline
double Collimator::getXEnd() {
    return xend_m;
}

inline
double Collimator::getYStart() {
    return ystart_m;
}

inline
double Collimator::getYEnd() {
    return yend_m;
}

inline
double Collimator::getZStart() {
    return zstart_m;
}

inline
double Collimator::getZEnd() {
    return zend_m;
}

inline
double Collimator::getWidth() {
    return width_m;
}

inline
void Collimator::setRHole(double r) {
    rHole_m = r;
}

inline
void Collimator::setNHoles(unsigned int nx, unsigned int ny) {
    nHolesX_m = nx;
    nHolesY_m = ny;
}

inline
void Collimator::setPitch(double p) {
    pitch_m = p;
}

inline
void Collimator::setPepperPot() {
    isAPepperPot_m = true;
}

inline
void Collimator::setSlit() {
    isASlit_m = true;
}

inline
void Collimator::setRColl() {
    isARColl_m = true;
}

inline
void Collimator::setCColl() {
    isACColl_m = true;
}

inline
void Collimator::setWire() {
    isAWire_m = true;
}

#endif // CLASSIC_Collimator_HH
