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

#include "Physics/Physics.h"
#include "AbsBeamline/Degrader.h"
#include "Algorithms/PartBunch.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "Fields/Fieldmap.hh"
#include "Structure/LossDataSink.h"
#include "H5hut.h"
#include <memory>

extern Inform *gmsg;

using namespace std;

// Class Degrader
// ------------------------------------------------------------------------

Degrader::Degrader():
    Component(),
    H5file_m(NULL),
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
    pitch_m(0.0)
{}


Degrader::Degrader(const Degrader &right):
    Component(right),
    H5file_m(NULL),
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
    pitch_m(right.pitch_m)
{

}


Degrader::Degrader(const string &name):
    Component(name),
    H5file_m(NULL),
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
    pitch_m(0.0)
{}


Degrader::~Degrader() {
}


void Degrader::accept(BeamlineVisitor &visitor) const {
    visitor.visitDegrader(*this);
}


inline bool Degrader::isInColl(Vector_t R, Vector_t P, double recpgamma) {
    Inform m ("Degrader::isInColl ");
 /**
     check if we are in the longitudinal
     range of the collimator
  */
  const double z = R(2) + P(2) * recpgamma;
  
  if((z > position_m) && (z <= position_m + getElementLength())) {
      const double trm1 = ((R(0)*R(0))/(getXsize()*getXsize()));
      const double trm2 = ((R(1)*R(1))/(getYsize()*getYsize()));                                 
      return (trm1 + trm2) > 1.0;
      m << R << endl;
  }
  return false;
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
    pdead = isInColl(R,P,recpgamma);

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
    return pdead;
}

bool Degrader::apply(const Vector_t &R, const Vector_t &centroid, const double &t, Vector_t &E, Vector_t &B) {
    return false;
}


// rectangle collimators in cyclotron cyclindral coordiantes
// without surfacephysics, the particle hitting collimator is deleted directly
bool Degrader::checkCollimator(PartBunch &bunch, const int turnnumber, const double t, const double tstep) {

    bool flagNeedUpdate = false;
    Vector_t rmin, rmax;

    bunch.get_bounds(rmin, rmax);
    
    if(rmax(2) >= zstart_m && rmin(2) <= zend_m) {
        flagNeedUpdate = true;  
        
        for(unsigned int i = 0; i < bunch.getLocalNum(); ++i) {
            if(bunch.PType[i] == 0 && bunch.R[i](2) < zend_m && bunch.R[i](2) > zstart_m ) {
                lossDs_m->addParticle(bunch.R[i], bunch.P[i], bunch.ID[i]);
                bunch.Bin[i] = -1;                    
            }
        }
    }
  reduce(&flagNeedUpdate, &flagNeedUpdate + 1, &flagNeedUpdate, OpBitwiseOrAssign());
  if(flagNeedUpdate) lossDs_m->save(getName());
  return flagNeedUpdate;
}


void Degrader::initialise(PartBunch *bunch, double &startField, double &endField, const double &scaleFactor) {
    RefPartBunch_m = bunch;
    position_m = startField;
    endField = position_m + getElementLength();
    // initialize DataSink with H5Part output enabled
    bool doH5 = false;
    lossDs_m = new LossDataSink(bunch->getTotalNum(), doH5);
    lossDs_m->openH5(getName());
}

void Degrader::initialise(PartBunch *bunch, const double &scaleFactor) {
    RefPartBunch_m = bunch;
    // initialize DataSink with H5Part output enabled
    bool doH5 = false;
    lossDs_m = new LossDataSink(bunch->getTotalNum(), doH5);
    lossDs_m->openH5(getName());
}


void Degrader::finalise()
{
  *gmsg << "Finalize probe" << endl;
  if(lossDs_m)
    delete lossDs_m;
}

void Degrader::goOnline() {
    Inform msg("DEGRADER ");
    if(RefPartBunch_m == NULL) {
        if(!informed_m) {
            string errormsg = Fieldmap::typeset_msg("BUNCH SIZE NOT SET", "warning");
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

    msg << "DEGRADER a= " << getXsize() << " b= " << b_m << " start= " << position_m << " fn= " << filename_m << " ny= " << nHolesY_m << " pitch= " << pitch_m << endl;

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
    h5_int64_t rc;
    reduce(online_m, online_m, OpOr());

    if(online_m) {
        online_m = false;
        if(filename_m == "") return;

        /**
           Check how much particles we have lost
        */
        int localLost = PosX_m.size();
        int globLost = 0;
        reduce(localLost, globLost, OpAddAssign());

        /**
           Check if we have to append of
           if we write a new file
        */
        if(globLost != 0) {
            ifstream inp;
            inp.open(filename_m.c_str(), ifstream::in);
            inp.close();
            if(inp.fail()) {
#ifdef PARALLEL_IO
                H5file_m = H5OpenFile(filename_m.c_str(), H5_O_WRONLY, Ippl::getComm());
#else
                H5file_m = H5OpenFile(filename_m.c_str(), H5_O_WRONLY, 0);
#endif
                rc = H5SetStep(H5file_m, 0);
                if(rc != H5_SUCCESS)
                    ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
                msg << "Lost " << globLost << " partices at collimator/pepperpot " << getName() << " write step 0" << endl;
                rc = H5WriteFileAttribString(H5file_m, "timeUnit", "s");
                if(rc != H5_SUCCESS)
                    ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
                rc = H5WriteFileAttribString(H5file_m, "xUnit", "m");
                if(rc != H5_SUCCESS)
                    ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
                rc = H5WriteFileAttribString(H5file_m, "yUnit", "m");
                if(rc != H5_SUCCESS)
                    ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
                rc = H5WriteFileAttribString(H5file_m, "zUnit", "m");
                if(rc != H5_SUCCESS)
                    ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
                rc = H5WriteFileAttribString(H5file_m, "pxUnit", "#beta#gamma");
                if(rc != H5_SUCCESS)
                    ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
                rc = H5WriteFileAttribString(H5file_m, "pyUnit", "#beta#gamma");
                if(rc != H5_SUCCESS)
                    ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
                rc = H5WriteFileAttribString(H5file_m, "pzUnit", "#beta#gamma");
                if(rc != H5_SUCCESS)
                    ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
            } else {
#ifdef PARALLEL_IO
                H5file_m = H5OpenFile(filename_m.c_str(), H5_O_APPEND, Ippl::getComm());
#else
                H5file_m = H5OpenFile(filename_m.c_str(), H5P_O_APPEND, 0);
#endif
                int numStepsInFile = H5GetNumSteps(H5file_m);
                rc = H5SetStep(H5file_m, numStepsInFile);
                if(rc != H5_SUCCESS)
                    ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
                msg << "Lost " << globLost << " partices at collimator/slit/pepperpot " << getName() << " append  step " << numStepsInFile << endl;
            }
        } else {
            msg << "collimator/pepperpot not used " << getName() << endl;
            return;
        }

        if(PosX_m.size() == 0) {
            rc = H5PartSetNumParticles(H5file_m, 0);
            if(rc != H5_SUCCESS)
                ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        } else {
            rc = H5PartSetNumParticles(H5file_m, PosX_m.size());
            if(rc != H5_SUCCESS)
                ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

            std::unique_ptr<char[]> varray(new char[PosX_m.size() * sizeof(double)]);
            double *fvalues = reinterpret_cast<double*>(varray.get());
            h5_int64_t *ids = reinterpret_cast<h5_int64_t*>(varray.get());

            int i = 0;
            vector<double>::iterator it;

            for(it = PosX_m.begin(); it != PosX_m.end(); ++it)
                fvalues[i++] = *it;
            rc = H5PartWriteDataFloat64(H5file_m, "x", fvalues);
            if(rc != H5_SUCCESS)
                ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

            i = 0;
            for(it = PosY_m.begin(); it != PosY_m.end(); ++it)
                fvalues[i++] = *it;
            rc = H5PartWriteDataFloat64(H5file_m, "y", fvalues);
            if(rc != H5_SUCCESS)
                ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

            i = 0;
            for(it = PosZ_m.begin(); it != PosZ_m.end(); ++it)
                fvalues[i++] = *it;
            rc = H5PartWriteDataFloat64(H5file_m, "z", fvalues);
            if(rc != H5_SUCCESS)
                ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

            i = 0;
            for(it = MomentumX_m.begin(); it != MomentumX_m.end(); ++it)
                fvalues[i++] = *it;
            rc = H5PartWriteDataFloat64(H5file_m, "px", fvalues);
            if(rc != H5_SUCCESS)
                ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

            i = 0;
            for(it = MomentumY_m.begin(); it != MomentumY_m.end(); ++it)
                fvalues[i++] = *it;
            rc = H5PartWriteDataFloat64(H5file_m, "py", fvalues);
            if(rc != H5_SUCCESS)
                ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

            i = 0;
            for(it = MomentumZ_m.begin(); it != MomentumZ_m.end(); ++it)
                fvalues[i++] = *it;
            rc = H5PartWriteDataFloat64(H5file_m, "pz", fvalues);
            if(rc != H5_SUCCESS)
                ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

            i = 0;
            for(it = time_m.begin(); it != time_m.end(); ++it)
                fvalues[i++] = *it;
            rc = H5PartWriteDataFloat64(H5file_m, "time", fvalues);
            if(rc != H5_SUCCESS)
                ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

            i = 0;
            for(vector<int>::iterator int_it = id_m.begin(); int_it != id_m.end(); ++int_it)
                ids[i++] = *int_it;
            rc = H5PartWriteDataInt64(H5file_m, "id", ids);
            if(rc != H5_SUCCESS)
                ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

        }
        rc = H5CloseFile(H5file_m);
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        PosX_m.clear();
        PosY_m.clear();
        PosZ_m.clear();
        MomentumX_m.clear();
        MomentumY_m.clear();
        MomentumZ_m.clear();
        time_m.clear();
        id_m.clear();
    }
}

bool Degrader::bends() const {
    return false;
}

void Degrader::setOutputFN(string fn) {
    filename_m = fn;
}

string Degrader::getOutputFN() {
    return  filename_m;
}

void Degrader::setXsize(double a) {
    a_m = a;
}

void Degrader::setYsize(double b) {
    b_m = b;
}

void Degrader::setXpos(double x0) {
    x0_m = x0;
}

void Degrader::setYpos(double y0) {
    y0_m = y0;
}


double Degrader::getXsize(double a) {
    return a_m;
}

double Degrader::getYsize(double b) {
    return b_m;
}

double Degrader::getXpos() {
    return x0_m;
}

double Degrader::getYpos() {
    return y0_m;

    // --------Cyclotron collimator
}
void Degrader::setXStart(double xstart) {
    xstart_m = xstart;
}

void Degrader::setXEnd(double xend) {
    xend_m = xend;
}

void Degrader::setYStart(double ystart) {
    ystart_m = ystart;
}

void Degrader::setYEnd(double yend) {
    yend_m = yend;
}

void Degrader::setZStart(double zstart) {
    zstart_m = zstart;
}

void Degrader::setZEnd(double zend) {
    zend_m = zend;
}

void Degrader::setWidth(double width) {
    width_m = width;
}

double Degrader::getXStart() {
    return xstart_m;
}

double Degrader::getXEnd() {
    return xend_m;
}

double Degrader::getYStart() {
    return ystart_m;
}

double Degrader::getYEnd() {
    return yend_m;
}

double Degrader::getZStart() {
    return zstart_m;
}

double Degrader::getZEnd() {
    return zend_m;
}

double Degrader::getWidth() {
    return width_m;
}

//-------------------------------

void Degrader::setRHole(double r) {
    rHole_m = r;
}
void Degrader::setNHoles(unsigned int nx, unsigned int ny) {
    nHolesX_m = nx;
    nHolesY_m = ny;
}
void Degrader::setPitch(double p) {
    pitch_m = p;
}


void Degrader::setPepperPot() {
    isAPepperPot_m = true;
}
void Degrader::setSlit() {
    isASlit_m = true;
}

void Degrader::setRColl() {
    isARColl_m = true;
}

void Degrader::setCColl() {
    isACColl_m = true;
}

void Degrader::setWire() {
    isAWire_m = true;
}
void Degrader::getDimensions(double &zBegin, double &zEnd) const {
    zBegin = position_m;
    zEnd = position_m + getElementLength();

}

const string &Degrader::getType() const {
    static const string type("DEGRADER");
    return type;
}

string Degrader::getDegraderShape() {
    return "DEGRADER";

}

