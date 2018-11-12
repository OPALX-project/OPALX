// ------------------------------------------------------------------------
// $RCSfile: CCollimator.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: CCollimator
//   Defines the abstract interface for a beam collimator for cyclotrons.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/CCollimator.h"
#include "Physics/Physics.h"
#include "Algorithms/PartBunchBase.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "Fields/Fieldmap.h"
#include "Structure/LossDataSink.h"
#include "Utilities/Options.h"
#include "Solvers/ParticleMatterInteractionHandler.hh"
#include "Utilities/Util.h"

#include <memory>

extern Inform *gmsg;

// Class Collimator
// ------------------------------------------------------------------------

CCollimator::CCollimator():
    Component(),
    filename_m(""),
    informed_m(false),
    losses_m(0),
    lossDs_m(nullptr),
    parmatint_m(NULL) {
    setDimensions(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
}


CCollimator::CCollimator(const CCollimator &right):
    Component(right),
    filename_m(right.filename_m),
    informed_m(right.informed_m),
    losses_m(0),
    lossDs_m(nullptr),
    parmatint_m(NULL)
{
    setDimensions(right.xstart_m, right.xend_m,
                  right.ystart_m, right.yend_m,
                  right.zstart_m, right.zend_m,
                  right.width_m);
    setGeom();
}


CCollimator::CCollimator(const std::string &name):
    Component(name),
    filename_m(""),
    informed_m(false),
    losses_m(0),
    lossDs_m(nullptr),
    parmatint_m(NULL) {
    setDimensions(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    setGeom();
}


CCollimator::~CCollimator() {
    if (online_m)
        goOffline();
}


void CCollimator::accept(BeamlineVisitor &visitor) const {
    visitor.visitCCollimator(*this);
}


bool CCollimator::apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B) {
    return false;
}

bool CCollimator::applyToReferenceParticle(const Vector_t &R, const Vector_t &P, const double &t, Vector_t &E, Vector_t &B) {
    return false;
}

bool CCollimator::checkCollimator(Vector_t r, Vector_t rmin, Vector_t rmax) {

    double r1 = std::hypot(rmax(0), rmax(1));
    bool isDead = false;
    if (rmax(2) >= zstart_m && rmin(2) <= zend_m) {
        if ( r1 > rstart_m - 10.0 && r1 < rend_m + 10.0 ){
            if (r(2) < zend_m && r(2) > zstart_m ) {
                int pflag = checkPoint(r(0), r(1));
                isDead = (pflag != 0);
            }
        }
    }
    return isDead;
}

// rectangle collimators in cyclotron cylindrical coordinates
// without particlematterinteraction, the particle hitting collimator is deleted directly
bool CCollimator::checkCollimator(PartBunchBase<double, 3> *bunch, const int /*turnnumber*/, const double /*t*/, const double /*tstep*/) {

    bool flagNeedUpdate = false;
    Vector_t rmin, rmax;
    bunch->get_bounds(rmin, rmax);

    if (rmax(2) >= zstart_m && rmin(2) <= zend_m) {
        // interested in absolute minimum and maximum
        double xmin = std::min(std::abs(rmin(0)), std::abs(rmax(0)));
        double xmax = std::max(std::abs(rmin(0)), std::abs(rmax(0)));
        double ymin = std::min(std::abs(rmin(1)), std::abs(rmax(1)));
        double ymax = std::max(std::abs(rmin(1)), std::abs(rmax(1)));
        double rbunch_min = std::hypot(xmin, ymin);
        double rbunch_max = std::hypot(xmax, ymax);

        if( rbunch_max > rmin_m && rbunch_min < rmax_m ){ // check similar to z
            size_t tempnum = bunch->getLocalNum();
            int pflag = 0;
            // now check each particle in bunch
            for (unsigned int i = 0; i < tempnum; ++i) {
                if (bunch->PType[i] == ParticleType::REGULAR && bunch->R[i](2) < zend_m && bunch->R[i](2) > zstart_m ) {
                    // only now careful check in r
                    pflag = checkPoint(bunch->R[i](0), bunch->R[i](1));
                    /// bunch->Bin[i] != -1 makes sure the particle is not stored in more than one collimator
                    if ((pflag != 0) && (bunch->Bin[i] != -1)) {
                        if (!parmatint_m)
                            lossDs_m->addParticle(bunch->R[i], bunch->P[i], bunch->ID[i]);
                        bunch->Bin[i] = -1;
                        flagNeedUpdate = true;
                    }
                }
            }
        }
    }
    reduce(&flagNeedUpdate, &flagNeedUpdate + 1, &flagNeedUpdate, OpBitwiseOrAssign());
    if (flagNeedUpdate && parmatint_m) {
        std::pair<Vector_t, double> boundingSphere;
        boundingSphere.first = 0.5 * (rmax + rmin);
        boundingSphere.second = euclidean_norm(rmax - boundingSphere.first);
        parmatint_m->apply(bunch, boundingSphere);
    }
    return flagNeedUpdate;
}

void CCollimator::initialise(PartBunchBase<double, 3> *bunch, double &startField, double &endField) {
    endField = startField + getElementLength();
    initialise(bunch);
}

void CCollimator::initialise(PartBunchBase<double, 3> *bunch) {
    RefPartBunch_m = bunch;

    parmatint_m = getParticleMatterInteraction();

    // if (!parmatint_m) {
    if (filename_m == std::string(""))
        lossDs_m = std::unique_ptr<LossDataSink>(new LossDataSink(getName(), !Options::asciidump));
    else
        lossDs_m = std::unique_ptr<LossDataSink>(new LossDataSink(filename_m.substr(0, filename_m.rfind(".")), !Options::asciidump));
    // }

    goOnline(-1e6);
}

void CCollimator::finalise()
{
    if (online_m)
        goOffline();
    *gmsg << "* Finalize cyclotron collimator " << getName() << endl;
}

void CCollimator::goOnline(const double &) {
    print();

    // PosX_m.reserve((int)(1.1 * RefPartBunch_m->getLocalNum()));
    // PosY_m.reserve((int)(1.1 * RefPartBunch_m->getLocalNum()));
    // PosZ_m.reserve((int)(1.1 * RefPartBunch_m->getLocalNum()));
    // MomentumX_m.reserve((int)(1.1 * RefPartBunch_m->getLocalNum()));
    // MomentumY_m.reserve((int)(1.1 * RefPartBunch_m->getLocalNum()));
    // MomentumZ_m.reserve((int)(1.1 * RefPartBunch_m->getLocalNum()));
    // time_m.reserve((int)(1.1 * RefPartBunch_m->getLocalNum()));
    // id_m.reserve((int)(1.1 * RefPartBunch_m->getLocalNum()));
    online_m = true;
}

void CCollimator::print() {
    if (RefPartBunch_m == NULL) {
        if (!informed_m) {
            std::string errormsg = Fieldmap::typeset_msg("BUNCH SIZE NOT SET", "warning");
            ERRORMSG(errormsg << endl);
            if (Ippl::myNode() == 0) {
                std::ofstream omsg("errormsg.txt", std::ios_base::app);
                omsg << errormsg << std::endl;
                omsg.close();
            }
            informed_m = true;
        }
        return;
    }

    *gmsg << "* CCollimator angle start " << xstart_m << " (Deg) angle end " << ystart_m << " (Deg) "
          << "R start " << xend_m << " (mm) R rend " << yend_m << " (mm)" << endl;
}

void CCollimator::goOffline() {
    if (online_m && lossDs_m)
        lossDs_m->save();
    lossDs_m.reset(0);
    online_m = false;
}

void CCollimator::setDimensions(double xstart, double xend, double ystart, double yend, double zstart, double zend, double width) {
    xstart_m = xstart;
    ystart_m = ystart;
    xend_m   = xend;
    yend_m   = yend;
    zstart_m = zstart;
    zend_m   = zend;
    width_m  = width;
    rstart_m = std::hypot(xstart, ystart);
    rend_m   = std::hypot(xend,   yend);
    // start position is the one with lowest radius
    if (rstart_m > rend_m) {
        std::swap(xstart_m, xend_m);
        std::swap(ystart_m, yend_m);
        std::swap(rstart_m, rend_m);
    }
    // zstart and zend are independent from x, y
    if (zstart_m > zend_m){
        std::swap(zstart_m, zend_m);
    }
}

bool CCollimator::bends() const {
    return false;
}

void CCollimator::setOutputFN(std::string fn) {
    filename_m = fn;
}

std::string CCollimator::getOutputFN() {
    if (filename_m == std::string(""))
        return getName();
    else
        return filename_m.substr(0, filename_m.rfind("."));
}

void CCollimator::getDimensions(double &zBegin, double &zEnd) const {
    zBegin = 0.0;
    zEnd = getElementLength();
}

ElementBase::ElementType CCollimator::getType() const {
    return CCOLLIMATOR;
}

std::string CCollimator::getCollimatorShape() {
    return "CCollimator";
}

void CCollimator::setGeom() {

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
    geom_m[2].y = yend_m - halfdist / coeff2;

    geom_m[3].x = xend_m - halfdist * coeff1;
    geom_m[3].y = yend_m + halfdist / coeff2;

    geom_m[4].x = geom_m[0].x;
    geom_m[4].y = geom_m[0].y;

    // calculate maximum and mininum r from these
    // current implementation not perfect minimum does not need to lie on corner, or in middle
    rmin_m = std::hypot(xstart_m, ystart_m);
    rmax_m = -std::numeric_limits<double>::max();
    for (int i=0; i<4; i++) {
        double r = std::hypot(geom_m[i].x, geom_m[i].y);
        rmin_m = std::min(rmin_m, r);
        rmax_m = std::max(rmax_m, r);
    }
    // some debug printout
    // for (int i=0; i<4; i++) {
    //     *gmsg << "point " << i << " ( " << geom_m[i].x << ", " << geom_m[i].y << ")" << endl;
    // }
    // *gmsg << "rmin " << rmin_m << " rmax " << rmax_m << endl;
}


int CCollimator::checkPoint(const double &x, const double &y) {
    int cn = 0;

    for (int i = 0; i < 4; i++) {
        if (((geom_m[i].y <= y) && (geom_m[i+1].y > y)) ||
            ((geom_m[i].y > y)  && (geom_m[i+1].y <= y))) {

            float vt = (float)(y - geom_m[i].y) / (geom_m[i+1].y - geom_m[i].y);
            if (x < geom_m[i].x + vt * (geom_m[i+1].x - geom_m[i].x))
                ++cn;
        }
    }
    return (cn & 1);  // 0 if even (out), and 1 if odd (in)
}