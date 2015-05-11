// ------------------------------------------------------------------------
// $RCSfile: Monitor.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Monitor
//   Defines the abstract interface for a beam position monitor.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------
#include "AbsBeamline/Monitor.h"
#include "Physics/Physics.h"
#include "Algorithms/PartBunch.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "Fields/Fieldmap.hh"
#include "Structure/LossDataSink.h"
#include "Utilities/Options.h"

#include <memory>

#include <fstream>
#include <memory>

using namespace std;

// Class Monitor
// ------------------------------------------------------------------------

Monitor::Monitor():
    Component(),
    filename_m(""),
    plane_m(OFF),
    position_m(0.0),
    informed_m(false)
{}


Monitor::Monitor(const Monitor &right):
    Component(right),
    filename_m(right.filename_m),
    plane_m(right.plane_m),
    position_m(right.position_m),
    informed_m(right.informed_m)
{}


Monitor::Monitor(const std::string &name):
    Component(name),
    filename_m(""),
    plane_m(OFF),
    position_m(0.0),
    informed_m(false)
{}


Monitor::~Monitor()
{}


void Monitor::accept(BeamlineVisitor &visitor) const {
    visitor.visitMonitor(*this);
}

bool Monitor::apply(const size_t &i, const double &t, double E[], double B[]) {
    Vector_t Ev(0, 0, 0), Bv(0, 0, 0);
    return apply(i, t, Ev, Bv);
}

bool Monitor::apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B) {
    const Vector_t &R = RefPartBunch_m->R[i];
    const Vector_t &P = RefPartBunch_m->P[i];
    const double recpgamma = Physics::c * RefPartBunch_m->getdT() / sqrt(1.0  + dot(P, P));
    if(online_m && R(2) < position_m && R(2) + P(2) * recpgamma > position_m) {
        double frac = (position_m - R(2)) / (P(2) * recpgamma);

        lossDs_m->addParticle(Vector_t(R(0) + frac * P(0) * recpgamma, R(1) + frac * P(1) * recpgamma, position_m),
                                   P, RefPartBunch_m->ID[i], t + frac * RefPartBunch_m->getdT(), 0);
    }

    return false;
}

bool Monitor::apply(const Vector_t &R, const Vector_t &centroid, const double &t, Vector_t &E, Vector_t &B) {
    return false;
}

void Monitor::initialise(PartBunch *bunch, double &startField, double &endField, const double &scaleFactor) {
    RefPartBunch_m = bunch;
    position_m = startField;
    startField -= 0.005;
    endField = position_m + 0.005;
    if (filename_m == std::string(""))
        lossDs_m = std::unique_ptr<LossDataSink>(new LossDataSink(getName(), !Options::asciidump));
    else
        lossDs_m = std::unique_ptr<LossDataSink>(new LossDataSink(filename_m.substr(0, filename_m.rfind(".")), !Options::asciidump));

}

void Monitor::finalise() {

}

void Monitor::goOnline(const double &) {
    if(RefPartBunch_m == NULL) {
        if(!informed_m) {
            Inform msg("Monitor ");
            std::string errormsg;
            errormsg = Fieldmap::typeset_msg("BUNCH SIZE NOT SET", "warning");
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

    if(Monitor::h5pfiles_s.find(filename_m) == Monitor::h5pfiles_s.end()) {
        Monitor::h5pfiles_s.insert(pair<string, unsigned int>(filename_m, 1));
        step_m = 0;
    } else {
        step_m = (*Monitor::h5pfiles_s.find(filename_m)).second ++;
    }
    online_m = true;
}

void Monitor::goOffline() {
    lossDs_m->save();
}

bool Monitor::bends() const {
    return false;
}

void Monitor::setOutputFN(std::string fn) {
    filename_m = fn;
}

void Monitor::getDimensions(double &zBegin, double &zEnd) const {
    zBegin = position_m - 0.005;
    zEnd = position_m + 0.005;
}


const std::string &Monitor::getType() const {
    static const std::string type("Monitor");
    return type;
}

void Monitor::moveBy(const double &dz) {
    position_m += dz;
}

map<string, unsigned int> Monitor::h5pfiles_s = map<string, unsigned int>();




    /*
    if (!Options::enableHDF5) return;

	reduce(online_m, online_m, OpOr());

    if(online_m) {
        online_m = false;
        if(filename_m == "") return;

        unsigned long nLoc = PosX_m.size();
        unsigned long i = 0;
        h5_file_t *H5file;
        h5_int64_t rc;
        if(step_m == 0) {
#ifdef PARALLEL_IO
            H5file = H5OpenFile(filename_m.c_str(), H5_O_WRONLY, Ippl::getComm());
#else
            H5file = H5OpenFile(filename_m.c_str(), H5_O_WRONLY, 0);
#endif
            rc = H5WriteFileAttribString(H5file, "timeUnit", "s");
            if(rc != H5_SUCCESS)
                ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
            rc = H5WriteFileAttribString(H5file, "xUnit", "m");
            if(rc != H5_SUCCESS)
                ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
            rc = H5WriteFileAttribString(H5file, "yUnit", "m");
            if(rc != H5_SUCCESS)
                ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
            rc = H5WriteFileAttribString(H5file, "pxUnit", "#beta#gamma");
            if(rc != H5_SUCCESS)
                ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
            rc = H5WriteFileAttribString(H5file, "pyUnit", "#beta#gamma");
            if(rc != H5_SUCCESS)
                ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
            rc = H5WriteFileAttribString(H5file, "pzUnit", "#beta#gamma");
            if(rc != H5_SUCCESS)
                ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
            rc = H5WriteFileAttribString(H5file, "SPOSUnit", "m");
            if(rc != H5_SUCCESS)
                ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        } else {
#ifdef PARALLEL_IO
            H5file = H5OpenFile(filename_m.c_str(), H5_O_APPEND, Ippl::getComm());
#else
            H5file = H5OpenFile(filename_m.c_str(), H5_O_APPEND, 0);
#endif
        }

        rc = H5SetStep(H5file, step_m);
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        rc = H5WriteStepAttribFloat64(H5file, "SPOS", &position_m, 1);
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        rc = H5PartSetNumParticles(H5file, PosX_m.size());
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

        std::unique_ptr<char> varray(new char[nLoc * sizeof(double)]);
        double *fvalues = reinterpret_cast<double*>(varray.get());
        h5_int64_t *ids = reinterpret_cast<h5_int64_t*>(varray.get());

	  FixMe: if I write with nLoc==0 -> rc == -2



	if (nLoc > 0) {
	  for(i = 0; i < nLoc; ++i) {
            fvalues[i] = PosX_m.front();
            PosX_m.pop_front();
	  }
	  rc = H5PartWriteDataFloat64(H5file, "x", fvalues);
	  if(rc != H5_SUCCESS)
	    ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << " nloc= " << nLoc << " fn= " << filename_m << endl);
	  for(i = 0; i < nLoc; ++i) {
            fvalues[i] = PosY_m.front();
            PosY_m.pop_front();
	  }
	  rc = H5PartWriteDataFloat64(H5file, "y", fvalues);
	  if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
	  for(i = 0; i < nLoc; ++i) {
            fvalues[i] = MomentumX_m.front();
            MomentumX_m.pop_front();
	  }
	  rc = H5PartWriteDataFloat64(H5file, "px", fvalues);
	  if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
	  for(i = 0; i < nLoc; ++i) {
            fvalues[i] = MomentumY_m.front();
            MomentumY_m.pop_front();
	  }
	  rc = H5PartWriteDataFloat64(H5file, "py", fvalues);
	  if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
	  for(i = 0; i < nLoc; ++i) {
            fvalues[i] = MomentumZ_m.front();
            MomentumZ_m.pop_front();
	  }
	  rc = H5PartWriteDataFloat64(H5file, "pz", fvalues);
	  if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
	  for(i = 0; i < nLoc; ++i) {
            fvalues[i] = time_m.front();
            time_m.pop_front();
	  }
	  rc = H5PartWriteDataFloat64(H5file, "time", fvalues);
	  if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
	  for(i = 0; i < nLoc; ++i) {
            ids[i] = id_m.front();
            id_m.pop_front();
	  }
	  rc = H5PartWriteDataInt64(H5file, "id", ids);
	  if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
	}
        rc = H5CloseFile(H5file);
        if(rc != H5_SUCCESS)
	  ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    }
    */
