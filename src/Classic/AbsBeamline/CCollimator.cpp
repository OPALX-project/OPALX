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
    xstart_m(0.0),
    xend_m(0.0),
    ystart_m(0.0),
    yend_m(0.0),
    width_m(0.0),
    losses_m(0),
    lossDs_m(nullptr),
    parmatint_m(NULL)
{}


CCollimator::CCollimator(const CCollimator &right):
    Component(right),
    filename_m(right.filename_m),
    informed_m(right.informed_m),
    xstart_m(right.xstart_m),
    xend_m(right.xend_m),
    ystart_m(right.ystart_m),
    yend_m(right.yend_m),
    zstart_m(right.zstart_m),
    zend_m(right.zend_m),
    width_m(right.width_m),
    losses_m(0),
    lossDs_m(nullptr),
    parmatint_m(NULL)
{
    setGeom();
}


CCollimator::CCollimator(const std::string &name):
    Component(name),
    filename_m(""),
    informed_m(false),
    xstart_m(0.0),
    xend_m(0.0),
    ystart_m(0.0),
    yend_m(0.0),
    zstart_m(0.0),
    zend_m(0.0),
    width_m(0.0),
    losses_m(0),
    lossDs_m(nullptr),
    parmatint_m(NULL)
{}


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

    double r_start = sqrt(xstart_m * xstart_m + ystart_m * ystart_m);
    double r_end = sqrt(xend_m * xend_m + yend_m * yend_m);
    double r1 = sqrt(rmax(0) * rmax(0) + rmax(1) * rmax(1));
    bool isDead = false;
    if (rmax(2) >= zstart_m && rmin(2) <= zend_m) {
        if ( r1 > r_start - 10.0 && r1 < r_end + 10.0 ){
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
bool CCollimator::checkCollimator(PartBunchBase<double, 3> *bunch, const int turnnumber, const double t, const double tstep) {

    bool flagNeedUpdate = false;
    Vector_t rmin, rmax;

    bunch->get_bounds(rmin, rmax);
    std::pair<Vector_t, double> boundingSphere;
    boundingSphere.first = 0.5 * (rmax + rmin);
    boundingSphere.second = euclidean_norm(rmax - boundingSphere.first);

    if (rmax(2) >= zstart_m && rmin(2) <= zend_m) {
        double r_start = sqrt(xstart_m * xstart_m + ystart_m * ystart_m);
        double r_end   = sqrt(  xend_m * xend_m   +   yend_m * yend_m);
        double rbunch_min = sqrt(rmin(0) * rmin(0) + rmin(1) * rmin(1));
        double rbunch_max = sqrt(rmax(0) * rmax(0) + rmax(1) * rmax(1));
        // if ( rbunch_max > rstart_m - 100.0 && rbunch_max < rend_m + 100.0 ){ // old check, only checks outer bound
        if( rbunch_max > r_start && rbunch_min < r_end ){ // check similar to z
            size_t tempnum = bunch->getLocalNum();
            int pflag = 0;
            for (unsigned int i = 0; i < tempnum; ++i) {
                if (bunch->PType[i] == ParticleType::REGULAR && bunch->R[i](2) < zend_m && bunch->R[i](2) > zstart_m ) {
                    pflag = checkPoint(bunch->R[i](0), bunch->R[i](1));
                    // bunch->Bin[i] != -1 makes sure the partcile is not stored in more than one collimator
                    if ((pflag != 0) && (bunch->Bin[i] != -1))  {
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
    *gmsg << "* Finalize probe" << endl;
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

    if (zstart_m > zend_m){
        double tempz = zstart_m;
        zstart_m = zend_m;
        zend_m = tempz;
    }
}


int CCollimator::checkPoint(const double &x, const double &y) {
    int cn = 0;

    for (int i = 0; i < 4; i++) {
        if (((geom_m[i].y <= y) && (geom_m[i+1].y > y))
            || ((geom_m[i].y > y) && (geom_m[i+1].y <= y))) {

            float vt = (float)(y - geom_m[i].y) / (geom_m[i+1].y - geom_m[i].y);
            if (x < geom_m[i].x + vt * (geom_m[i+1].x - geom_m[i].x))
                ++cn;
        }
    }
    return (cn & 1);  // 0 if even (out), and 1 if odd (in)
}