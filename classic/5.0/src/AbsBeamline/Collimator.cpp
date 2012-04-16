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

#include "Physics/Physics.h"
#include "AbsBeamline/Collimator.h"
#include "Algorithms/PartBunch.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "Fields/Fieldmap.hh"
#include "Structure/LossDataSink.h"
#include <memory>

extern Inform *gmsg;

using namespace std;

// Class Collimator
// ------------------------------------------------------------------------

Collimator::Collimator():
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
    as_m(0.0),
    rs_m(0.0),
    ae_m(0.0),
    re_m(0.0),
    w_m(0.0),
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


Collimator::Collimator(const Collimator &right):
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
    as_m(right.as_m),
    rs_m(right.rs_m),
    ae_m(right.ae_m),
    re_m(right.re_m),
    w_m(right.w_m),
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
  setGeom();
}


Collimator::Collimator(const string &name):
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
    as_m(0.0),
    rs_m(0.0),
    ae_m(0.0),
    re_m(0.0),
    w_m(0.0),
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


Collimator::~Collimator() {
}


void Collimator::accept(BeamlineVisitor &visitor) const {
    visitor.visitCollimator(*this);
}

bool Collimator::apply(const int &i, const double &t, double E[], double B[]) {
    Vector_t Ev(0, 0, 0), Bv(0, 0, 0);
    return apply(i, t, Ev, Bv);
}

bool Collimator::apply(const int &i, const double &t, Vector_t &E, Vector_t &B) {
    const Vector_t &R = RefPartBunch_m->R[i] - Vector_t(dx_m, dy_m, ds_m); // including the missaligment
    const Vector_t &P = RefPartBunch_m->P[i];
    const double recpgamma = Physics::c * RefPartBunch_m->getdT() / sqrt(1.0  + dot(P, P));

    /**
    check if we are in the longitudinal
    range of the collimator
    */

    const double z = R(2) + P(2) * recpgamma;

    // particle is not dead
    bool pdead = false;

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
            pdead = !alive;
        } else if(isASlit_m) {
            //      if ( (abs(R(0) >= getXsize()) || (abs(R(1) >= getYsize()))))
            if(R(0) <= -getXsize() || R(1) <= -getYsize() || R(0) >= getXsize() || R(1) >= getYsize())
                pdead = true;
        } else {
            // case of an elliptic collimator
        }




        if(pdead) {
            double frac = (R(2) - position_m) / P(2) * recpgamma;
            PosX_m.push_back(R(0));
            PosY_m.push_back(R(1));
            PosZ_m.push_back(z);
            MomentumX_m.push_back(P(0));
            MomentumY_m.push_back(P(1));
            MomentumZ_m.push_back(P(2));
            time_m.push_back(t + frac * RefPartBunch_m->getdT());
            id_m.push_back(i);
        }
    }
    return pdead;
}

bool Collimator::apply(const Vector_t &R, const Vector_t &centroid, const double &t, Vector_t &E, Vector_t &B) {
    return false;
}


// rectangle collimators in cyclotron cyclindral coordiantes
// without surfacephysics, the particle hitting collimator is deleted directly
bool Collimator::checkCollimator(PartBunch &bunch, const int turnnumber, const double t, const double tstep) {

    bool flagNeedUpdate = false;
    Vector_t rmin, rmax;
    bunch.get_bounds(rmin, rmax);
    double r1 = sqrt(rmax(0) * rmax(0) + rmax(1) * rmax(1));

    if(r1 > rs_m - 50.0 && r1 < re_m + 50.0 ){
      size_t tempnum = bunch.getLocalNum();
      int pflag = 0;

      for(unsigned int i = 0; i < tempnum; ++i) {
	if(bunch.PType[i] == 0) {
	  pflag = checkPoint(bunch.R[i](0), bunch.R[i](1));
	  if(pflag != 0) {
	    lossDs_m->addParticle(bunch.R[i], bunch.P[i], bunch.ID[i]);
	    bunch.Bin[i] = -1;
	    flagNeedUpdate = true;
	  }
	}
      }
    }
    reduce(&flagNeedUpdate, &flagNeedUpdate + 1, &flagNeedUpdate, OpBitwiseOrAssign());
    if(flagNeedUpdate) lossDs_m->save(getName());
    return flagNeedUpdate;
}


void Collimator::initialise(PartBunch *bunch, double &startField, double &endField, const double &scaleFactor) {
    RefPartBunch_m = bunch;
    position_m = startField;
    endField = position_m + getElementLength();
    // initialize DataSink with H5Part output enabled
    bool doH5 = false;
    lossDs_m = new LossDataSink(bunch->getTotalNum(), doH5);
    lossDs_m->openH5(getName());
}

void Collimator::initialise(PartBunch *bunch, const double &scaleFactor) {
    RefPartBunch_m = bunch;
    // initialize DataSink with H5Part output enabled
    bool doH5 = false;
    lossDs_m = new LossDataSink(bunch->getTotalNum(), doH5);
    lossDs_m->openH5(getName());
}


void Collimator::finalise()
{
  *gmsg << "Finalize probe" << endl;
  if(lossDs_m)
    delete lossDs_m;
}

void Collimator::goOnline() {
    Inform msg("Collimator ");
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

    if(isAPepperPot_m)
        msg << "Pepperpot x= " << a_m << " y= " << b_m << " r= " << rHole_m << " nx= " << nHolesX_m << " ny= " << nHolesY_m << " pitch= " << pitch_m << endl;
    else if(isASlit_m)
        msg << "Slit x= " << getXsize() << " Slit y= " << getYsize() << " start= " << position_m << " fn= " << filename_m << endl;
    else if(isARColl_m)
        msg << "RCollimator a= " << getXsize() << " b= " << b_m << " start= " << position_m << " fn= " << filename_m << " ny= " << nHolesY_m << " pitch= " << pitch_m << endl;
    else if(isACColl_m)
        msg << "CCollimator angstart= " << as_m << " angend " << ae_m << " rstart " << rs_m << " rend " << re_m << endl;
    else if(isAWire_m)
        msg << "Wire x= " << x0_m << " y= " << y0_m << endl;
    else
        msg << "ECollimator a= " << getXsize() << " b= " << b_m << " start= " << position_m << " fn= " << filename_m << " ny= " << nHolesY_m << " pitch= " << pitch_m << endl;

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

void Collimator::goOffline() {
    Inform msg("Collimator::goOffline ");
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

bool Collimator::bends() const {
    return false;
}

void Collimator::setOutputFN(string fn) {
    filename_m = fn;
}

string Collimator::getOutputFN() {
    return  filename_m;
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
void Collimator::setAngStart(double as) {
    as_m = as;
}

void Collimator::setRStart(double rs) {
    rs_m = rs;
}

void Collimator::setAngEnd(double ae) {
    ae_m = ae;
}

void Collimator::setREnd(double re) {
    re_m = re;
}

void Collimator::setWidth(double w) {
    w_m = w;
}

double Collimator::getAngStart() {
    return as_m;
}

double Collimator::getRStart() {
    return rs_m;
}

double Collimator::getAngEnd() {
    return ae_m;
}

double Collimator::getREnd() {
    return re_m;
}

double Collimator::getWidth() {
    return w_m;
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

const string &Collimator::getType() const {
    static const string type("Collimator");
    return type;
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

    geom_m[0].x = rs_m*cos(as_m);
    geom_m[0].y = rs_m*sin(as_m);

    geom_m[1].x = re_m*cos(as_m);
    geom_m[1].y = re_m*sin(as_m);

    geom_m[2].x = re_m*cos(ae_m);
    geom_m[2].y = re_m*sin(ae_m);

    geom_m[3].x = rs_m*cos(ae_m);
    geom_m[3].y = rs_m*sin(ae_m);

    geom_m[4].x = geom_m[0].x;
    geom_m[4].y = geom_m[0].y;


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
