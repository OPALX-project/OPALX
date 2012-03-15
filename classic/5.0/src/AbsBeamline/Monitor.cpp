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

#include "Physics/Physics.h"
#include "AbsBeamline/Monitor.h"
#include "Fields/Fieldmap.hh"
#include "AbsBeamline/BeamlineVisitor.h"


// Class Monitor
// ------------------------------------------------------------------------

Monitor::Monitor():
    Component(),
    filename_m(""),
    plane_m(OFF),
    position_m(0.0),
    PosX_m(0),
    PosY_m(0),
    MomentumX_m(0),
    MomentumY_m(0),
    MomentumZ_m(0),
    time_m(0),
    id_m(0),
    informed_m(false)
{}


Monitor::Monitor(const Monitor &right):
    Component(right),
    filename_m(right.filename_m),
    plane_m(right.plane_m),
    position_m(right.position_m),
    PosX_m(right.PosX_m),
    PosY_m(right.PosY_m),
    MomentumX_m(right.MomentumX_m),
    MomentumY_m(right.MomentumY_m),
    MomentumZ_m(right.MomentumZ_m),
    time_m(right.time_m),
    id_m(right.id_m),
    informed_m(right.informed_m)
{}


Monitor::Monitor(const string &name):
    Component(name),
    filename_m(""),
    plane_m(OFF),
    position_m(0.0),
    PosX_m(0),
    PosY_m(0),
    MomentumX_m(0),
    MomentumY_m(0),
    MomentumZ_m(0),
    time_m(0),
    id_m(0),
    informed_m(false)
{}


Monitor::~Monitor()
{}


void Monitor::accept(BeamlineVisitor &visitor) const {
    visitor.visitMonitor(*this);
}

bool Monitor::apply(const int &i, const double &t, double E[], double B[]) {
    Vector_t Ev(0, 0, 0), Bv(0, 0, 0);
    return apply(i, t, Ev, Bv);
}

bool Monitor::apply(const int &i, const double &t, Vector_t &E, Vector_t &B) {
    const Vector_t &R = RefPartBunch_m->R[i];
    const Vector_t &P = RefPartBunch_m->P[i];
    const double recpgamma = Physics::c * RefPartBunch_m->getdT() / sqrt(1.0  + dot(P, P));
    if(online_m && R(2) < position_m && R(2) + P(2) * recpgamma > position_m) {
        double frac = (position_m - R(2)) / (P(2) * recpgamma);
        PosX_m.push_back(R(0) + frac * P(0) * recpgamma);
        PosY_m.push_back(R(1) + frac * P(1) * recpgamma);
        MomentumX_m.push_back(P(0));
        MomentumY_m.push_back(P(1));
        MomentumZ_m.push_back(P(2));
        time_m.push_back(t + frac * RefPartBunch_m->getdT());
        id_m.push_back(i);
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
}

void Monitor::finalise() {

}

void Monitor::goOnline() {
    if(RefPartBunch_m == NULL) {
        if(!informed_m) {
            Inform msg("Monitor ");
            string errormsg;
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

    if (Monitor::h5pfiles_s.find(filename_m) == Monitor::h5pfiles_s.end()) {
        Monitor::h5pfiles_s.insert(pair<string, unsigned int>(filename_m, 1));
        step_m = 0;
    } else {
        step_m = (*Monitor::h5pfiles_s.find(filename_m)).second ++;
    }
    online_m = true;
}

void Monitor::goOffline() {
    reduce(online_m, online_m, OpOr());

    if(online_m) {
        online_m = false;
        if(filename_m == "") return;
        
        unsigned long nLoc = PosX_m.size();
        unsigned long i = 0;
        H5PartFile *H5file;
        if (step_m == 0) {
#ifdef PARALLEL_IO
            H5file = H5PartOpenFileParallel(filename_m.c_str(), H5PART_WRITE, MPI_COMM_WORLD);
#else
            H5file = H5PartOpenFile(filename_m.c_str(), H5PART_WRITE);
#endif
            H5PartWriteFileAttribString(H5file, "timeUnit", "s");
            H5PartWriteFileAttribString(H5file, "xUnit", "m");
            H5PartWriteFileAttribString(H5file, "yUnit", "m");
            H5PartWriteFileAttribString(H5file, "pxUnit", "#beta#gamma");
            H5PartWriteFileAttribString(H5file, "pyUnit", "#beta#gamma");
            H5PartWriteFileAttribString(H5file, "pzUnit", "#beta#gamma");
            H5PartWriteFileAttribString(H5file, "SPOSUnit", "m");
        } else {
#ifdef PARALLEL_IO
            H5file = H5PartOpenFileParallel(filename_m.c_str(), H5PART_APPEND, MPI_COMM_WORLD);
#else
            H5file = H5PartOpenFile(filename_m.c_str(), H5PART_APPEND);
#endif
        }

        H5PartSetStep(H5file, step_m);
        H5PartWriteStepAttrib(H5file, "SPOS", H5T_NATIVE_DOUBLE, &position_m, 1);
        H5PartSetNumParticles(H5file, PosX_m.size());

        void *varray = malloc(nLoc * sizeof(double));
        double *fvalues = (double *)varray;
        h5part_int64_t *ids = (h5part_int64_t *)varray;

        for(i = 0; i < nLoc; ++i) {
            fvalues[i] = PosX_m.front();
            PosX_m.pop_front();
        }
        H5PartWriteDataFloat64(H5file, "x", fvalues);

        for(i = 0; i < nLoc; ++i) {
            fvalues[i] = PosY_m.front();
            PosY_m.pop_front();
        }
        H5PartWriteDataFloat64(H5file, "y", fvalues);

        for(i = 0; i < nLoc; ++i) {
            fvalues[i] = MomentumX_m.front();
            MomentumX_m.pop_front();
        }
        H5PartWriteDataFloat64(H5file, "px", fvalues);

        for(i = 0; i < nLoc; ++i) {
            fvalues[i] = MomentumY_m.front();
            MomentumY_m.pop_front();
        }
        H5PartWriteDataFloat64(H5file, "py", fvalues);

        for(i = 0; i < nLoc; ++i) {
            fvalues[i] = MomentumZ_m.front();
            MomentumZ_m.pop_front();
        }
        H5PartWriteDataFloat64(H5file, "pz", fvalues);

        for(i = 0; i < nLoc; ++i) {
            fvalues[i] = time_m.front();
            time_m.pop_front();
        }
        H5PartWriteDataFloat64(H5file, "time", fvalues);

        for(i = 0; i < nLoc; ++i) {
            ids[i] = id_m.front();
            id_m.pop_front();
        }
        H5PartWriteDataInt64(H5file, "id", ids);

        H5Fflush(H5file->file, H5F_SCOPE_GLOBAL);

        free(varray);

        H5PartCloseFile(H5file);
    }
}

bool Monitor::bends() const {
    return false;
}

void Monitor::setOutputFN(string fn) {
    filename_m = fn;
}

void Monitor::getDimensions(double &zBegin, double &zEnd) const {
    zBegin = position_m - 0.005;
    zEnd = position_m + 0.005;
}


const string &Monitor::getType() const {
    static const string type("Monitor");
    return type;
}

void Monitor::moveBy(const double & dz) {
    position_m += dz;
}

map<string, unsigned int> Monitor::h5pfiles_s = map<string, unsigned int>();
