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
    position_m(0.0),
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
    H5file_m(NULL),
    filename_m(right.filename_m),
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
    zstart_m(right.zstart_m),
    zend_m(right.zend_m)
{}

Degrader::Degrader(const string &name):
    Component(name),
    H5file_m(NULL),
    filename_m(""),
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
    zstart_m(0.0),
    zend_m(0.0)
{}


Degrader::~Degrader() {
}


void Degrader::accept(BeamlineVisitor &visitor) const {
    visitor.visitDegrader(*this);
}


inline bool Degrader::isInMaterial(double z ) {
 /**
     check if the particle is in the degarder material
   
  */
  return ((z > position_m) && (z <= position_m + getElementLength()));
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
 Inform msg("Degrader::goOnline ");
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

const string &Degrader::getType() const {
    static const string type("DEGRADER");
    return type;
}

string Degrader::getDegraderShape() {
    return "DEGRADER";

}
