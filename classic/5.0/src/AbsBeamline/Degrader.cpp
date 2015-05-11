// ------------------------------------------------------------------------
// $RCSfile: Degrader.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Degrader
//   Defines the abstract interface for a beam Degrader.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Degrader.h"
#include "Algorithms/PartBunch.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "Fields/Fieldmap.hh"
#include "Structure/LossDataSink.h"
#include "Utilities/Options.h"
#include "Physics/Physics.h"
#include <memory>

extern Inform *gmsg;

using namespace std;

// Class Degrader
// ------------------------------------------------------------------------

Degrader::Degrader():
    Component(),
    filename_m(""),
    position_m(0.0),
    deg_width_m(0.0),
    PosX_m(0),
    PosY_m(0),
    PosZ_m(0),
    MomentumX_m(0),
    MomentumY_m(0),
    MomentumZ_m(0),
    time_m(0),
    id_m(0),
    informed_m(false)
{}

Degrader::Degrader(const Degrader &right):
    Component(right),
    filename_m(right.filename_m),
    position_m(right.position_m),
    deg_width_m(right.deg_width_m),
    PosX_m(right.PosX_m),
    PosY_m(right.PosY_m),
    PosZ_m(right.PosZ_m),
    MomentumX_m(right.MomentumX_m),
    MomentumY_m(right.MomentumY_m),
    MomentumZ_m(right.MomentumZ_m),
    time_m(right.time_m),
    id_m(right.id_m),
    informed_m(right.informed_m),
    zstart_m(right.zstart_m),
    zend_m(right.zend_m)
{}

Degrader::Degrader(const std::string &name):
    Component(name),
    filename_m(""),
    position_m(0.0),
    deg_width_m(0.0),
    PosX_m(0),
    PosY_m(0),
    PosZ_m(0),
    MomentumX_m(0),
    MomentumY_m(0),
    MomentumZ_m(0),
    time_m(0),
    id_m(0),
    informed_m(false),
    zstart_m(0.0),
    zend_m(0.0)
{}


Degrader::~Degrader() {

  if(online_m)
    goOffline();
}


void Degrader::accept(BeamlineVisitor &visitor) const {
    visitor.visitDegrader(*this);
}


inline bool Degrader::isInMaterial(double z ) {
 /**
     check if the particle is in the degarder material

  */
    return ((z > position_m) && (z <= position_m + getZSize())); //getElementLength()));
}

bool Degrader::apply(const size_t &i, const double &t, double E[], double B[]) {
    Vector_t Ev(0, 0, 0), Bv(0, 0, 0);
    return apply(i, t, Ev, Bv);
}

bool Degrader::apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B) {

    const Vector_t &R = RefPartBunch_m->R[i] - Vector_t(dx_m, dy_m, ds_m); // including the missaligment
    const Vector_t &P = RefPartBunch_m->P[i];
    const double recpgamma = Physics::c * RefPartBunch_m->getdT() / sqrt(1.0  + dot(P, P));

    bool pdead = false;
    pdead = isInMaterial(R(2));

    if(pdead) {
      RefPartBunch_m->Bin[i] = -1;
      double frac = (R(2) - position_m) / P(2) * recpgamma;
      PosX_m.push_back(R(0));
      PosY_m.push_back(R(1));
      PosZ_m.push_back(R(2));
      MomentumX_m.push_back(P(0));
      MomentumY_m.push_back(P(1));
      MomentumZ_m.push_back(P(2));
      time_m.push_back(t + frac * RefPartBunch_m->getdT());
      id_m.push_back(i);
    }
    return false;
}

bool Degrader::apply(const Vector_t &R, const Vector_t &centroid, const double &t, Vector_t &E, Vector_t &B) {
    return false;
}


void Degrader::initialise(PartBunch *bunch, double &startField, double &endField, const double &scaleFactor) {
    RefPartBunch_m = bunch;
    position_m = startField;
    endField = position_m + getElementLength();

    if (filename_m == std::string(""))
        lossDs_m = std::unique_ptr<LossDataSink>(new LossDataSink(getName(), !Options::asciidump));
    else
        lossDs_m = std::unique_ptr<LossDataSink>(new LossDataSink(filename_m.substr(0, filename_m.rfind(".")), !Options::asciidump));
}

void Degrader::initialise(PartBunch *bunch, const double &scaleFactor) {
    RefPartBunch_m = bunch;
    if (filename_m == std::string(""))
        lossDs_m = std::unique_ptr<LossDataSink>(new LossDataSink(getName(), !Options::asciidump));
    else
        lossDs_m = std::unique_ptr<LossDataSink>(new LossDataSink(filename_m.substr(0, filename_m.rfind(".")), !Options::asciidump));
}


void Degrader::finalise()
{
  *gmsg << "* Finalize Degrader" << endl;
}

void Degrader::goOnline(const double &) {
 Inform msg("Degrader::goOnline ");
   if(RefPartBunch_m == NULL) {
        if(!informed_m) {
            std::string errormsg = Fieldmap::typeset_msg("BUNCH SIZE NOT SET", "warning");
            msg << errormsg << "\n"
                << endl;
            if(Ippl::myNode() == 0) {
                ofstream omsg("errormsg.txt", ios_base::app);
                omsg << errormsg << endl;
                omsg.close();
            }
            informed_m = true;
        }
        return;
    }
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

void Degrader::goOffline() {
    Inform msg("Degrader::goOffline ");
    online_m = false;
    lossDs_m->save();
    msg << " done..." << endl;
}

bool Degrader::bends() const {
    return false;
}

void Degrader::setOutputFN(std::string fn) {
    filename_m = fn;
}

string Degrader::getOutputFN() {
    if (filename_m == std::string(""))
        return getName();
    else
        return filename_m.substr(0, filename_m.rfind("."));
}

double Degrader::getZStart() {
    return zstart_m;
}

double Degrader::getZEnd() {
    return zend_m;
}

void Degrader::getDimensions(double &zBegin, double &zEnd) const {
    zBegin = position_m;
    zEnd = position_m + getElementLength();

}

const std::string &Degrader::getType() const {
    static const std::string type("DEGRADER");
    return type;
}

string Degrader::getDegraderShape() {
    return "DEGRADER";

}

void Degrader::setZSize(double z) {
    deg_width_m = z;
}

double Degrader::getZSize() {
    return deg_width_m;
}