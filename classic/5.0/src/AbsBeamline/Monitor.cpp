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
#include "AbsBeamline/BeamlineVisitor.h"


// Class Monitor
// ------------------------------------------------------------------------

Monitor::Monitor():
  Component(),
  H5file_m(NULL),
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
  H5file_m(NULL),
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
  H5file_m(NULL),
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


void Monitor::accept(BeamlineVisitor &visitor) const
{
  visitor.visitMonitor(*this);
}

bool Monitor::apply(const int &i, const double &t, double E[], double B[])
{
  Vector_t Ev(0,0,0), Bv(0,0,0);
  return apply(i,t,Ev,Bv);
}

bool Monitor::apply(const int &i, const double &t, Vector_t &E, Vector_t &B)
{
  const Vector_t &R = RefPartBunch_m->R[i];
  const Vector_t &P = RefPartBunch_m->P[i];
  const double recpgamma = Physics::c * RefPartBunch_m->getdT() / sqrt(1.0  + dot(P,P));
  if (R(2) < position_m && R(2) + P(2) * recpgamma > position_m)
    { 
      double frac = (R(2) - position_m) / P(2) * recpgamma;
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

bool Monitor::apply(const Vector_t &R, const Vector_t &centroid, const double &t, Vector_t &E, Vector_t &B)
{ 
  return false;
}

void Monitor::initialise(const PartBunch *bunch, double &startField, double &endField, const double &scaleFactor)
{
  RefPartBunch_m = bunch;
  position_m = startField;
  startField -= 0.005;
  endField = position_m + 0.005;
}

void Monitor::finalise()
{

}

void Monitor::goOnline()
{
  if (RefPartBunch_m == NULL)
    {
      if (!informed_m)
        {
          Inform msg("Monitor ");
          msg << "* ************** W A R N I N G *****************************************************" << endl;
          msg << "* BUNCH SIZE NOT SET " << endl;
          msg << "* **********************************************************************************" << endl;
          informed_m = true;
        }
      return;
    }
  
  PosX_m.reserve((int)(1.1 * RefPartBunch_m->getLocalNum()));
  PosY_m.reserve((int)(1.1 * RefPartBunch_m->getLocalNum()));
  MomentumX_m.reserve((int)(1.1 * RefPartBunch_m->getLocalNum()));
  MomentumY_m.reserve((int)(1.1 * RefPartBunch_m->getLocalNum()));
  MomentumZ_m.reserve((int)(1.1 * RefPartBunch_m->getLocalNum()));
  time_m.reserve((int)(1.1 * RefPartBunch_m->getLocalNum()));
  id_m.reserve((int)(1.1 * RefPartBunch_m->getLocalNum()));
//   PosX_m.resize(RefPartBunch_m->getLocalNum());
//   PosY_m.resize(RefPartBunch_m->getLocalNum());
//   MomentumX_m.resize(RefPartBunch_m->getLocalNum());
//   MomentumY_m.resize(RefPartBunch_m->getLocalNum());
//   MomentumZ_m.resize(RefPartBunch_m->getLocalNum());
//   time_m.resize(RefPartBunch_m->getLocalNum());
//   id_m.resize(RefPartBunch_m->getLocalNum());
  online_m = true;
}

void Monitor::goOffline()
{

  reduce(online_m,online_m,OpOr());

  if (online_m) {
    online_m = false;
    if (filename_m == "") return;

#ifdef PARALLEL_IO
    H5file_m=H5PartOpenFileParallel(filename_m.c_str(),H5PART_WRITE,MPI_COMM_WORLD);
#else
    H5file_m=H5PartOpenFile(filename_m.c_str(),H5PART_WRITE);
#endif
    H5PartWriteFileAttribString(H5file_m,"timeUnit","s");
    H5PartWriteFileAttribString(H5file_m,"xUnit","m");
    H5PartWriteFileAttribString(H5file_m,"yUnit","m");
    H5PartWriteFileAttribString(H5file_m,"pxUnit","#beta#gamma");
    H5PartWriteFileAttribString(H5file_m,"pyUnit","#beta#gamma");
    H5PartWriteFileAttribString(H5file_m,"pzUnit","#beta#gamma");

    H5PartSetStep(H5file_m,0);
    H5PartSetNumParticles(H5file_m,PosX_m.size()); 

    void *varray = malloc(PosX_m.size()*sizeof(double));
    double *fvalues = (double*)varray;
    h5part_int64_t *ids = (h5part_int64_t *)varray;
    
    int i = 0;
    vector<double>::iterator it;

    for (it = PosX_m.begin(); it != PosX_m.end(); ++it)
      fvalues[i++] = *it;
    H5PartWriteDataFloat64(H5file_m,"x",fvalues); 
  
    i = 0;
    for (it = PosY_m.begin(); it != PosY_m.end(); ++it)
      fvalues[i++] = *it;
    H5PartWriteDataFloat64(H5file_m,"y",fvalues); 
    
    i = 0;
    for (it = MomentumX_m.begin(); it != MomentumX_m.end(); ++it)
      fvalues[i++] = *it;
    H5PartWriteDataFloat64(H5file_m,"px",fvalues); 

    i = 0;
    for (it = MomentumY_m.begin(); it != MomentumY_m.end(); ++it)
      fvalues[i++] = *it;
    H5PartWriteDataFloat64(H5file_m,"py",fvalues); 

    i = 0;
    for (it = MomentumZ_m.begin(); it != MomentumZ_m.end(); ++it)
      fvalues[i++] = *it;
    H5PartWriteDataFloat64(H5file_m,"pz",fvalues); 

    i = 0;
    for (it = time_m.begin(); it != time_m.end(); ++it)
      fvalues[i++] = *it;
    H5PartWriteDataFloat64(H5file_m,"time",fvalues); 

    i = 0;
    for (vector<int>::iterator int_it = id_m.begin(); int_it != id_m.end(); ++int_it)
      ids[i++] = *int_it;
    H5PartWriteDataInt64(H5file_m,"id",ids);  
    
    H5Fflush(H5file_m->file,H5F_SCOPE_GLOBAL);
    
    free(varray);

    H5PartCloseFile(H5file_m);  
    PosX_m.clear();
    PosY_m.clear();
    MomentumX_m.clear();
    MomentumY_m.clear();
    MomentumZ_m.clear();
    time_m.clear();
    id_m.clear();
  }
}

void Monitor::rescaleFieldMap(const double &scaleFactor)
{}

bool Monitor::bends() const
{
  return false;
}

void Monitor::setOutputFN(string fn)
{
  filename_m = fn;
}

void Monitor::getDimensions(double &zBegin, double &zEnd) const
{
  zBegin = position_m - 0.005;
  zEnd = position_m + 0.005;
}


const string& Monitor::getType() const
{
    static const string type("Monitor");
    return type;
}

