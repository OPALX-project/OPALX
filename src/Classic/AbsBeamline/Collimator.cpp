// ------------------------------------------------------------------------
// $RCSfile: Collimator.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Collimator
//   Defines the abstract interface for a beam Collimator.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Collimator.h"
#include "Physics/Physics.h"
#include "Algorithms/PartBunch.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "Fields/Fieldmap.h"
#include "Structure/LossDataSink.h"
#include "Utilities/Options.h"
#include "Solvers/SurfacePhysicsHandler.hh"
#include <memory>

extern Inform *gmsg;

using namespace std;

// Class Collimator
// ------------------------------------------------------------------------

Collimator::Collimator():
    Component(),
    filename_m(""),
    plane_m(OFF),
    position_m(0.0),
    PosX_m(0),
    PosY_m(0),
    PosZ_m(0),
    MomentumX_m(0),
    MomentumY_m(0),
    MomentumZ_m(0),
    time_m(0),
    id_m(0),
    informed_m(false),
    a_m(0.0),
    b_m(0.0),
    x0_m(0.0),
    y0_m(0.0),
    xstart_m(0.0),
    xend_m(0.0),
    ystart_m(0.0),
    yend_m(0.0),
    width_m(0.0),
    isAPepperPot_m(false),
    isASlit_m(false),
    isARColl_m(false),
    isACColl_m(false),
    isAWire_m(false),
    rHole_m(0.0),
    nHolesX_m(0),
    nHolesY_m(0),
    pitch_m(0.0),
    losses_m(0),
    lossDs_m(nullptr),
    sphys_m(NULL)
{}


Collimator::Collimator(const Collimator &right):
    Component(right),
    filename_m(right.filename_m),
    plane_m(right.plane_m),
    position_m(right.position_m),
    PosX_m(right.PosX_m),
    PosY_m(right.PosY_m),
    PosZ_m(right.PosZ_m),
    MomentumX_m(right.MomentumX_m),
    MomentumY_m(right.MomentumY_m),
    MomentumZ_m(right.MomentumZ_m),
    time_m(right.time_m),
    id_m(right.id_m),
    informed_m(right.informed_m),
    a_m(right.a_m),
    b_m(right.b_m),
    x0_m(right.x0_m),
    y0_m(right.y0_m),
    xstart_m(right.xstart_m),
    xend_m(right.xend_m),
    ystart_m(right.ystart_m),
    yend_m(right.yend_m),
    zstart_m(right.zstart_m),
    zend_m(right.zend_m),
    width_m(right.width_m),
    isAPepperPot_m(right.isAPepperPot_m),
    isASlit_m(right.isASlit_m),
    isARColl_m(right.isARColl_m),
    isACColl_m(right.isACColl_m),
    isAWire_m(right.isAWire_m),
    rHole_m(right.rHole_m),
    nHolesX_m(right.nHolesX_m),
    nHolesY_m(right.nHolesY_m),
    pitch_m(right.pitch_m),
    losses_m(0),
    lossDs_m(nullptr),
    sphys_m(NULL)
{
    setGeom();
}


Collimator::Collimator(const std::string &name):
    Component(name),
    filename_m(""),
    plane_m(OFF),
    position_m(0.0),
    PosX_m(0),
    PosY_m(0),
    PosZ_m(0),
    MomentumX_m(0),
    MomentumY_m(0),
    MomentumZ_m(0),
    time_m(0),
    id_m(0),
    informed_m(false),
    a_m(0.0),
    b_m(0.0),
    x0_m(0.0),
    y0_m(0.0),
    xstart_m(0.0),
    xend_m(0.0),
    ystart_m(0.0),
    yend_m(0.0),
    zstart_m(0.0),
    zend_m(0.0),
    width_m(0.0),
    isAPepperPot_m(false),
    isASlit_m(false),
    isARColl_m(false),
    isACColl_m(false),
    isAWire_m(false),
    rHole_m(0.0),
    nHolesX_m(0),
    nHolesY_m(0),
    pitch_m(0.0),
    losses_m(0),
    lossDs_m(nullptr),
    sphys_m(NULL)
{}


Collimator::~Collimator() {
    if(online_m)
        goOffline();
}


void Collimator::accept(BeamlineVisitor &visitor) const {
    visitor.visitCollimator(*this);
}


inline bool Collimator::isInColl(Vector_t R, Vector_t P, double recpgamma) {
    /**
       check if we are in the longitudinal
       range of the collimator
    */
    const double z = R(2) + P(2) * recpgamma;

    if((z > position_m) && (z <= position_m + getElementLength())) {
        if(isAPepperPot_m) {

            /**
               ------------
               |(0)|  |(0)|
               ----   -----
               |    a)    |
               |          |
               ----   -----
               |(0)|  |(0)|
               yL------------
               xL
               |---| d
               |--| pitch
               Observation: the area in a) is much larger than the
               area(s) (0). In a) particles are lost in (0)
               particles they are not lost.

            */
            const double h  =   pitch_m;
            const double xL = - 0.5 * h * (nHolesX_m - 1);
            const double yL = - 0.5 * h * (nHolesY_m - 1);
            bool alive = false;

            for(unsigned int m = 0; (m < nHolesX_m && (!alive)); m++) {
                for(unsigned int n = 0; (n < nHolesY_m && (!alive)); n++) {
                    double x_m = xL  + (m * h);
                    double y_m = yL  + (n * h);
                    /** are we in a) ? */
                    double rr = std::pow((R(0) - x_m) / rHole_m, 2) + std::pow((R(1) - y_m) / rHole_m, 2);
                    alive = (rr < 1.0);
                }
            }
            return !alive;
        } else if(isASlit_m || isARColl_m) {
            return (std::abs(R(0)) >= getXsize() || std::abs(R(1)) >= getYsize());
        } else if(isAWire_m) {
            ERRORMSG("Not yet implemented");
        } else {
            // case of an elliptic collimator
            return (std::pow(R(0) / getXsize(), 2) + std::pow(R(1) / getYsize(), 2)) >= 1.0;
        }
    }
    return false;
}

bool Collimator::apply(const size_t &i, const double &t, double E[], double B[]) {
    Vector_t Ev(0, 0, 0), Bv(0, 0, 0);
    return apply(i, t, Ev, Bv);
}

bool Collimator::apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B) {
    const Vector_t &R = RefPartBunch_m->R[i] - Vector_t(dx_m, dy_m, ds_m); // including the missaligment
    const Vector_t &P = RefPartBunch_m->P[i];
    const double recpgamma = Physics::c * RefPartBunch_m->getdT() / sqrt(1.0  + dot(P, P));
    bool pdead = false;
    pdead = isInColl(R,P,recpgamma);
    if(pdead) {
      if (lossDs_m) {
	double frac = (R(2) - position_m) / P(2) * recpgamma;
	lossDs_m->addParticle(R,P, RefPartBunch_m->ID[i], t + frac * RefPartBunch_m->getdT(), 0);
      }
      ++losses_m;
    }
    return pdead;
}

bool Collimator::apply(const Vector_t &R, const Vector_t &centroid, const double &t, Vector_t &E, Vector_t &B) {
    return false;
}

bool Collimator::checkCollimator(Vector_t r, Vector_t rmin, Vector_t rmax) {

    double r_start = sqrt(xstart_m * xstart_m + ystart_m * ystart_m);
    double r_end = sqrt(xend_m * xend_m + yend_m * yend_m);
    double r1 = sqrt(rmax(0) * rmax(0) + rmax(1) * rmax(1));
    bool isDead = false;
    if(rmax(2) >= zstart_m && rmin(2) <= zend_m) {
        if( r1 > r_start - 10.0 && r1 < r_end + 10.0 ){
            if(r(2) < zend_m && r(2) > zstart_m ) {
                int pflag = checkPoint(r(0), r(1));
                isDead = (pflag != 0);
            }
        }
    }
    return isDead;
}


// rectangle collimators in cyclotron cyclindral coordiantes
// without surfacephysics, the particle hitting collimator is deleted directly
bool Collimator::checkCollimator(PartBunch &bunch, const int turnnumber, const double t, const double tstep) {

    bool flagNeedUpdate = false;
    Vector_t rmin, rmax;

    bunch.get_bounds(rmin, rmax);
    double r_start = sqrt(xstart_m * xstart_m + ystart_m * ystart_m);
    double r_end = sqrt(xend_m * xend_m + yend_m * yend_m);
    double r1 = sqrt(rmax(0) * rmax(0) + rmax(1) * rmax(1));

    if(rmax(2) >= zstart_m && rmin(2) <= zend_m) {
//        if( r1 > r_start - 10.0 && r1 < r_end + 10.0 ){
        if( r1 > r_start - 100.0 && r1 < r_end + 100.0 ){
            size_t tempnum = bunch.getLocalNum();
            int pflag = 0;
            for(unsigned int i = 0; i < tempnum; ++i) {
                if(bunch.PType[i] == ParticleType::REGULAR && bunch.R[i](2) < zend_m && bunch.R[i](2) > zstart_m ) {
                    pflag = checkPoint(bunch.R[i](0), bunch.R[i](1));
		    /// bunch.Bin[i] != -1 makes sure the partcile is not stored in more than one collimator
                    if ((pflag != 0) && (bunch.Bin[i] != -1))  {
                        if (!sphys_m)
                            lossDs_m->addParticle(bunch.R[i], bunch.P[i], bunch.ID[i]);
                        bunch.Bin[i] = -1;
                        flagNeedUpdate = true;
                    }
                }
            }
        }
    }
    reduce(&flagNeedUpdate, &flagNeedUpdate + 1, &flagNeedUpdate, OpBitwiseOrAssign());
    if (flagNeedUpdate && sphys_m) {
        sphys_m->apply(bunch);
    }
    return flagNeedUpdate;
}

void Collimator::initialise(PartBunch *bunch, double &startField, double &endField, const double &scaleFactor) {
    RefPartBunch_m = bunch;
    position_m = startField;
    endField = position_m + getElementLength();

    sphys_m = getSurfacePhysics();

    //if (!sphys_m) {
    if (filename_m == std::string(""))
      lossDs_m = std::unique_ptr<LossDataSink>(new LossDataSink(getName(), !Options::asciidump));
    else
      lossDs_m = std::unique_ptr<LossDataSink>(new LossDataSink(filename_m.substr(0, filename_m.rfind(".")), !Options::asciidump));
      //}

    goOnline(-1e6);
}

void Collimator::initialise(PartBunch *bunch, const double &scaleFactor) {
    RefPartBunch_m = bunch;

    sphys_m = getSurfacePhysics();

    //    if (!sphys_m) {
    if (filename_m == std::string(""))
      lossDs_m = std::unique_ptr<LossDataSink>(new LossDataSink(getName(), !Options::asciidump));
    else
      lossDs_m = std::unique_ptr<LossDataSink>(new LossDataSink(filename_m.substr(0, filename_m.rfind(".")), !Options::asciidump));
    //    }

    goOnline(-1e6);
}

void Collimator::finalise()
{
    if(online_m)
        goOffline();
    if(isAPepperPot_m)
      *gmsg << "* Finalize pepperpot " << getName() << endl;
    else if(isASlit_m)
      *gmsg << "* Finalize slit " << getName() << endl;
    else if(isARColl_m)
      *gmsg << "* Finalize collimator " << getName() << endl;
    else if(isACColl_m)
      *gmsg << "* Finalize cyclotron collimator " << getName() << endl;
    else if(isAWire_m)
      *gmsg << "* Finalize wire " << getName() << endl;
    else
      *gmsg << "* Finalize ecollimator " << getName() << endl;
}

void Collimator::goOnline(const double &) {
    print();

    PosX_m.reserve((int)(1.1 * RefPartBunch_m->getLocalNum()));
    PosY_m.reserve((int)(1.1 * RefPartBunch_m->getLocalNum()));
    PosZ_m.reserve((int)(1.1 * RefPartBunch_m->getLocalNum()));
    MomentumX_m.reserve((int)(1.1 * RefPartBunch_m->getLocalNum()));
    MomentumY_m.reserve((int)(1.1 * RefPartBunch_m->getLocalNum()));
    MomentumZ_m.reserve((int)(1.1 * RefPartBunch_m->getLocalNum()));
    time_m.reserve((int)(1.1 * RefPartBunch_m->getLocalNum()));
    id_m.reserve((int)(1.1 * RefPartBunch_m->getLocalNum()));
    online_m = true;
}

void Collimator::print() {
    if(RefPartBunch_m == NULL) {
        if(!informed_m) {
            std::string errormsg = Fieldmap::typeset_msg("BUNCH SIZE NOT SET", "warning");
            ERRORMSG(errormsg << endl);
            if(Ippl::myNode() == 0) {
                ofstream omsg("errormsg.txt", ios_base::app);
                omsg << errormsg << endl;
                omsg.close();
            }
            informed_m = true;
        }
        return;
    }

    *gmsg << level1;
    if(isAPepperPot_m)
      *gmsg << "* Pepperpot " + getName() + " x= " << a_m << " y= " << b_m << " r= " << rHole_m
              << " nx= " << nHolesX_m << " ny= " << nHolesY_m << " pitch= " << pitch_m << endl;
    else if(isASlit_m)
        *gmsg << "* Slit " + getName() + " x= " << getXsize() << " Slit y= " << getYsize()
              << " start= " << position_m << " fn= " << filename_m << endl;
    else if(isARColl_m)
      *gmsg << "* RCollimator " + getName() + " a= " << getXsize() << " b= " << getYsize()
              << " start= " << position_m << " fn= " << filename_m
              << " ny= " << nHolesY_m << " pitch= " << pitch_m << endl;
    else if(isACColl_m)
      *gmsg << "* CCollimator " + getName() + " angle start " << xstart_m << " angle end " << ystart_m << " all (Deg) "
              << "R start " << xend_m << " R rend " << yend_m << " all (mm)" << endl;
    else if(isAWire_m)
      *gmsg << "* Wire " + getName() + " x= " << x0_m << " y= " << y0_m << endl;
    else
      *gmsg << "* ECollimator "+ getName() + " a= " << getXsize() << " b= " << b_m << " start= " << position_m
              << " fn= " << filename_m << " ny= " << nHolesY_m << " pitch= " << pitch_m << endl;
}

void Collimator::goOffline() {
    if (online_m && lossDs_m)
        lossDs_m->save();
    online_m = false;
}

bool Collimator::bends() const {
    return false;
}

void Collimator::setOutputFN(std::string fn) {
    filename_m = fn;
}

string Collimator::getOutputFN() {
    if (filename_m == std::string(""))
        return getName();
    else
        return filename_m.substr(0, filename_m.rfind("."));
}

unsigned int Collimator::getLosses() const {
    return losses_m;
}

void Collimator::setXsize(double a) {
    a_m = a;
}

void Collimator::setYsize(double b) {
    b_m = b;
}

void Collimator::setXpos(double x0) {
    x0_m = x0;
}

void Collimator::setYpos(double y0) {
    y0_m = y0;
}


double Collimator::getXsize(double a) {
    return a_m;
}

double Collimator::getYsize(double b) {
    return b_m;
}

double Collimator::getXpos() {
    return x0_m;
}

double Collimator::getYpos() {
    return y0_m;

    // --------Cyclotron collimator
}
void Collimator::setXStart(double xstart) {
    xstart_m = xstart;
}

void Collimator::setXEnd(double xend) {
    xend_m = xend;
}

void Collimator::setYStart(double ystart) {
    ystart_m = ystart;
}

void Collimator::setYEnd(double yend) {
    yend_m = yend;
}

void Collimator::setZStart(double zstart) {
    zstart_m = zstart;
}

void Collimator::setZEnd(double zend) {
    zend_m = zend;
}

void Collimator::setWidth(double width) {
    width_m = width;
}

double Collimator::getXStart() {
    return xstart_m;
}

double Collimator::getXEnd() {
    return xend_m;
}

double Collimator::getYStart() {
    return ystart_m;
}

double Collimator::getYEnd() {
    return yend_m;
}

double Collimator::getZStart() {
    return zstart_m;
}

double Collimator::getZEnd() {
    return zend_m;
}

double Collimator::getWidth() {
    return width_m;
}

//-------------------------------

void Collimator::setRHole(double r) {
    rHole_m = r;
}
void Collimator::setNHoles(unsigned int nx, unsigned int ny) {
    nHolesX_m = nx;
    nHolesY_m = ny;
}
void Collimator::setPitch(double p) {
    pitch_m = p;
}


void Collimator::setPepperPot() {
    isAPepperPot_m = true;
}
void Collimator::setSlit() {
    isASlit_m = true;
}

void Collimator::setRColl() {
    isARColl_m = true;
}

void Collimator::setCColl() {
    isACColl_m = true;
}

void Collimator::setWire() {
    isAWire_m = true;
}
void Collimator::getDimensions(double &zBegin, double &zEnd) const {
    zBegin = position_m;
    zEnd = position_m + getElementLength();

    // zBegin = position_m - 0.005;
    //  zEnd = position_m + 0.005;

}

ElementBase::ElementType Collimator::getType() const {
    return COLLIMATOR;
}

string Collimator::getCollimatorShape() {
    if(isAPepperPot_m)
        return "PeperPot";
    else if(isASlit_m)
        return "Slit";
    else if(isARColl_m)
        return "RCollimator";
    else if(isACColl_m)
        return "CCollimator";
    else if(isAWire_m)
        return "Wire";
    else
        return "ECollimator";
}

void Collimator::setGeom() {

    double slope;
    if (xend_m == xstart_m)
        slope = 1.0e12;
    else
        slope = (yend_m - ystart_m) / (xend_m - xstart_m);

    double coeff2 = sqrt(1 + slope * slope);
    double coeff1 = slope / coeff2;
    double halfdist = width_m / 2.0;
    geom_m[0].x = xstart_m - halfdist * coeff1;
    geom_m[0].y = ystart_m + halfdist / coeff2;

    geom_m[1].x = xstart_m + halfdist * coeff1;
    geom_m[1].y = ystart_m - halfdist / coeff2;

    geom_m[2].x = xend_m + halfdist * coeff1;
    geom_m[2].y = yend_m - halfdist  / coeff2;

    geom_m[3].x = xend_m - halfdist * coeff1;
    geom_m[3].y = yend_m + halfdist / coeff2;

    geom_m[4].x = geom_m[0].x;
    geom_m[4].y = geom_m[0].y;

    if (zstart_m > zend_m){
        double tempz = 0.0;
        tempz = zstart_m;
        zstart_m = zend_m;
        zend_m = tempz;
    }
}


int Collimator::checkPoint(const double &x, const double &y) {
    int    cn = 0;

    for(int i = 0; i < 4; i++) {
        if(((geom_m[i].y <= y) && (geom_m[i+1].y > y))
           || ((geom_m[i].y > y) && (geom_m[i+1].y <= y))) {

            float vt = (float)(y - geom_m[i].y) / (geom_m[i+1].y - geom_m[i].y);
            if(x < geom_m[i].x + vt * (geom_m[i+1].x - geom_m[i].x))
                ++cn;
        }
    }
    return (cn & 1);  // 0 if even (out), and 1 if odd (in)
}