#include "AbsBeamline/FlexibleCollimator.h"
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

// Class FlexibleCollimator
// ------------------------------------------------------------------------

FlexibleCollimator::FlexibleCollimator():
    Component(),
    filename_m(""),
    informed_m(false),
    losses_m(0),
    lossDs_m(nullptr),
    parmatint_m(NULL){}


FlexibleCollimator::FlexibleCollimator(const FlexibleCollimator &right):
    Component(right),
    //holes_m(right.holes_m.begin(), right.holes_m.end()),
    bb_m(right.bb_m),
    tree_m(right.tree_m),
    filename_m(right.filename_m),
    informed_m(right.informed_m),
    losses_m(0),
    losses1_m(0),
    losses2_m(0),
    losses3_m(0),
    lossDs_m(nullptr),
    parmatint_m(NULL)
{
    for (const mslang::Base *obj: right.holes_m) {
        holes_m.push_back(obj->clone());
    }
}


FlexibleCollimator::FlexibleCollimator(const std::string &name):
    Component(name),
    filename_m(""),
    informed_m(false),
    losses_m(0),
    lossDs_m(nullptr),
    parmatint_m(NULL)
{}


FlexibleCollimator::~FlexibleCollimator() {
    if (online_m)
        goOffline();

    for (mslang::Base *obj: holes_m) {
        delete obj;
    }
}


void FlexibleCollimator::accept(BeamlineVisitor &visitor) const {
    visitor.visitFlexibleCollimator(*this);
}


bool FlexibleCollimator::isStopped(const Vector_t &R, const Vector_t &P, double recpgamma) {
    const double z = R(2);

    if ((z < 0.0) ||
        (z > getElementLength()) ||
        (!isInsideTransverse(R))) return false;

    ++ losses3_m;
    if (!bb_m.isInside(R)) {
        ++losses1_m;
        return true;
    }

    if (!tree_m.isInside(R)) {
        ++ losses2_m;
        return true;
    }

    return  false;
}

bool FlexibleCollimator::apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B) {
    const Vector_t &R = RefPartBunch_m->R[i];
    const Vector_t &P = RefPartBunch_m->P[i];
    const double &dt = RefPartBunch_m->dt[i];
    const double recpgamma = Physics::c * dt / sqrt(1.0 + dot(P, P));
    bool pdead = isStopped(R, P, recpgamma);

    if (pdead) {
        if (lossDs_m) {
            double frac = -R(2) / P(2) * recpgamma;
            lossDs_m->addParticle(R, P,
                                  RefPartBunch_m->ID[i],
                                  t + frac * dt, 0);
        }
        ++ losses_m;
    }

    return pdead;
}

bool FlexibleCollimator::applyToReferenceParticle(const Vector_t &R, const Vector_t &P, const double &t, Vector_t &E, Vector_t &B) {
    return false;
}

// rectangle collimators in cyclotron cyclindral coordiantes
// without particlematterinteraction, the particle hitting collimator is deleted directly
bool FlexibleCollimator::checkCollimator(PartBunchBase<double, 3> *bunch, const int turnnumber, const double t, const double tstep) {

    return false;
}

void FlexibleCollimator::initialise(PartBunchBase<double, 3> *bunch, double &startField, double &endField) {
    RefPartBunch_m = bunch;
    endField = startField + getElementLength();

    parmatint_m = getParticleMatterInteraction();

    // if (!parmatint_m) {
    if (filename_m == std::string(""))
        lossDs_m = std::unique_ptr<LossDataSink>(new LossDataSink(getName(), !Options::asciidump));
    else
        lossDs_m = std::unique_ptr<LossDataSink>(new LossDataSink(filename_m.substr(0, filename_m.rfind(".")), !Options::asciidump));
    // }

    goOnline(-1e6);
}

void FlexibleCollimator::initialise(PartBunchBase<double, 3> *bunch) {
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

void FlexibleCollimator::finalise()
{
    if (online_m)
        goOffline();
    *gmsg << "* Finalize probe" << endl;
}

void FlexibleCollimator::goOnline(const double &) {
    print();

    online_m = true;
}

void FlexibleCollimator::print() {
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

    *gmsg << level3;
}

void FlexibleCollimator::goOffline() {
    if (online_m && lossDs_m)
        lossDs_m->save();
    lossDs_m.reset(0);
    online_m = false;
}

void FlexibleCollimator::getDimensions(double &zBegin, double &zEnd) const {
    zBegin = 0.0;
    zEnd = getElementLength();
}

ElementBase::ElementType FlexibleCollimator::getType() const {
    return COLLIMATOR;
}

void FlexibleCollimator::setDescription(const std::string &desc) {
    mslang::Function *fun;

    if (!mslang::parse(desc, fun))
        throw GeneralClassicException("FlexibleCollimator::setDescription",
                                      "Couldn't parse input file");

    fun->apply(holes_m);

    if (holes_m.size() == 0) return;

    for (auto it: holes_m) {
        it->computeBoundingBox();
    }

    mslang::Base *first = holes_m.front();
    const mslang::BoundingBox &bb = first->bb_m;

    Vector_t llc(bb.center_m[0] - 0.5 * bb.width_m,
                 bb.center_m[1] - 0.5 * bb.height_m,
                 0.0);
    Vector_t urc(bb.center_m[0] + 0.5 * bb.width_m,
                 bb.center_m[1] + 0.5 * bb.height_m,
                 0.0);

    for (auto it: holes_m) {
        const mslang::BoundingBox &bb = it->bb_m;
        llc[0] = std::min(llc[0], bb.center_m[0] - 0.5 * bb.width_m);
        llc[1] = std::min(llc[1], bb.center_m[1] - 0.5 * bb.height_m);
        urc[0] = std::max(urc[0], bb.center_m[0] + 0.5 * bb.width_m);
        urc[1] = std::max(urc[1], bb.center_m[1] + 0.5 * bb.height_m);
    }

    bb_m = mslang::BoundingBox(llc, urc);

    tree_m.bb_m = bb_m;
    tree_m.objects_m.insert(tree_m.objects_m.end(), holes_m.begin(), holes_m.end());
    tree_m.buildUp();

    std::ofstream out("quadtree.gpl");
    tree_m.writeGnuplot(out);
}