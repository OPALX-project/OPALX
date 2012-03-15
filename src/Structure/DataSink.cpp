// ------------------------------------------------------------------------
// $RCSfile: DataSink.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.3 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: DataSink
//
// ------------------------------------------------------------------------
//
// Revision History:
// $Date: 2004/06/02 19:38:54 $
// $Author: adelmann $
// $Log: DataSink.cpp,v $
// Revision 1.3  2004/06/02 19:38:54  adelmann
// add tabs in the stat file
//
// Revision 1.2  2003/12/12 10:46:39  adelmann
// Remove old stuff
//
// Revision 1.1.1.2  2003/07/04 12:52:00  adelmann
// July 4 2003
//
// Revision 3.2  2001/09/08 16:34:08  adelmann
// Add writeout of: moments, dispersion
//
// Revision 3.1  2001/08/30 11:28:26  adelmann
// Remove .madprogress and write to gnuplot
// add courant snyder parametres to statfile
//
// Revision 3.0  2001/08/22 14:41:33  adelmann
// The stable Version
//
// Revision 2.17  2001/08/11 05:30:22  adelmann
// Production Version
//
// Revision 1.6  2001/08/11 05:21:47  adelmann
// August 11 2001 V2.17
//
// Add #-tag to string in statfile
//
// Revision 1.5  2001/07/05 20:52:53  adelmann
// Mad9p V2.12
//
// Used reduction on one processor only -> program hangs
//
// Revision 1.4  2001/02/23 12:28:34  adelmann
// Add space between the column of the STAMARKER
//
// Revision 1.3  2001/02/11 17:14:27  adelmann
// Remove .mad9pprogress
//
// Revision 1.2  2001/01/26 12:42:49  adelmann
// This is MAD9p Version: Tue Jan 25 2001 V1.28
//
// Revision 1.1.1.1  2000/11/30 20:29:48  adelmann
// g++ and KCPP
//
// Revision 1.5  2000/11/09 15:24:24  adelmann
// - no output in .progress file anymore
//   the function is NOT removed so one can use
//   this to dump some debug information
//
// - collect all reference data drom PartData
//
// Revision 1.4  2000/11/05 04:24:47  adelmann
// Add mean valuse of momentum to the statistics file LANL Nov 2000
//
// Revision 1.3  2000/10/26 05:17:58  adelmann
// Remove DX stuff and add timestam and title
//
// Revision 1.2  2000/08/10 10:56:34  adelmann
// Some cleanup and add a option dx (data explorer) !
// The new printall skript extracts data from this stat file format
//
// Revision 1.1.1.1  2000/07/14 07:20:54  adelmann
// linux version Fri Jul 14 09:15:27 CEST 2000
//
// Revision 1.1.1.1  2000/05/20 11:13:58  adelmann
// Initial working version, thick elements, without: field dump and collimators
//
// Revision 1.3  2000/01/28 07:22:40  adelmann
// Fixrd some bugs with the Mad9pOutput
//
// Revision 1.2  2000/01/27 14:13:27  adelmann
// - Add  bunch->dataSink_m.saveStatDataGnuplotFormat( . )
//        DiscParticle write out
//        updateDotProgress
//
// Revision 1.1.1.1  2000/01/06 07:33:27  adelmann
// linux version works with gcc 991007 V2.96
//
// Revision 2.2  1999/10/29 05:02:05  adelmann
// *** empty log message ***
//
// Revision 2.1  1999/10/27 06:36:59  adelmann
// SGI-LINUX g++, with RF-Gap, REVSCATTER and read distribution from file
//
// Revision 1.1.1.1  1999/10/26 04:22:18  adelmann
// Classic 2.1 (with p)
//
// Revision 1.1.1.1  1999/10/26 04:14:36  adelmann
// Classic 2.1 (with p)
//
//
// ------------------------------------------------------------------------

#include "DataSink.h"
#include "AbstractObjects/OpalData.h"
#include "Physics/Physics.h"
#include "ValueDefinitions/RealVariable.h"

#include <hdf5.h>
#include "H5Part.h"
extern "C" {
#include "H5Block.h"
}

using namespace std;
//using Physics::m_e;

/////////////////////////////////////////////////////////////////////////////

DataSink::DataSink()
{

  H5PartTimer_m = IpplTimings::getTimer("H5PartTimer");
  StatMarkerTimer_m = IpplTimings::getTimer("StatMarkerTimer");

  // Because we write only once a header into the statfile. This flag gets
  // reset on the first write to the statfile
  firstWriteToStat_m = true;
  firstWriteH5part_m = true;
 
  string fn = OPAL.getInputFn(); 

  INFOMSG(fn<<endl);

  int pos=fn.find(string("."),0);  
  fn.erase(pos,fn.size()-pos);
  
  fn += string(".h5");

#ifdef PARALLEL_IO
  H5file_m=H5PartOpenFileParallel(fn.c_str(),H5PART_WRITE,MPI_COMM_WORLD);
#else
  H5file_m=H5PartOpenFile(fn.c_str(),H5PART_WRITE);
#endif

  if(!H5file_m) {
    ERRORMSG("File open failed:  exiting!" << endl);
    exit(0);
  }

  H5PartWriteFileAttribString(H5file_m,"tUnit","s");
  H5PartWriteFileAttribString(H5file_m,"xUnit","m");
  H5PartWriteFileAttribString(H5file_m,"yUnit","m");
  H5PartWriteFileAttribString(H5file_m,"zUnit","m");
  H5PartWriteFileAttribString(H5file_m,"pxUnit","#beta#gamma");
  H5PartWriteFileAttribString(H5file_m,"pyUnit","#beta#gamma");
  H5PartWriteFileAttribString(H5file_m,"pzUnit","#beta#gamma");
  H5PartWriteFileAttribString(H5file_m,"idUnit","1");
  H5PartWriteFileAttribString(H5file_m,"SPOSUnit","m");
  H5PartWriteFileAttribString(H5file_m,"TIMEUnit","s");
  H5PartWriteFileAttribString(H5file_m,"#gammaUnit","1");
  H5PartWriteFileAttribString(H5file_m,"ENERGYUnit","MeV");
  H5PartWriteFileAttribString(H5file_m,"#varepsilonUnit","m rad");
  H5PartWriteFileAttribString(H5file_m,"#varepsilonrUnit","m rad");

  H5PartWriteFileAttribString(H5file_m,"#varepsilon-geomUnit","m rad");

  H5PartWriteFileAttribString(H5file_m,"#sigmaUnit"," ");
  H5PartWriteFileAttribString(H5file_m,"RMSXUnit","m");
  H5PartWriteFileAttribString(H5file_m,"RMSRUnit","m");
  H5PartWriteFileAttribString(H5file_m,"RMSPUnit","#beta#gamma");

  H5PartWriteFileAttribString(H5file_m,"maxdEUnit","MeV");
  H5PartWriteFileAttribString(H5file_m,"max#phiUnit","deg");

  H5PartWriteFileAttribString(H5file_m,"phizUnit","deg");
  H5PartWriteFileAttribString(H5file_m,"enezUnit","keV");

  // for fields of head/ref particle/tail
  H5PartWriteFileAttribString(H5file_m,"spos-headUnit","m");
  H5PartWriteFileAttribString(H5file_m,"E-headUnit","MV/m");
  H5PartWriteFileAttribString(H5file_m,"B-headUnit","T");
  H5PartWriteFileAttribString(H5file_m,"spos-refUnit","m");
  H5PartWriteFileAttribString(H5file_m,"E-refUnit","MV/m");
  H5PartWriteFileAttribString(H5file_m,"B-refUnit","T");
  H5PartWriteFileAttribString(H5file_m,"spos-tailUnit","m");
  H5PartWriteFileAttribString(H5file_m,"E-tailUnit","MV/m");
  H5PartWriteFileAttribString(H5file_m,"B-tailUnit","T");
  H5PartWriteFileAttribString(H5file_m,"LPATHUnit","m");
  H5PartWriteFileAttribString(H5file_m,"StepUnit"," ");
  H5PartWriteFileAttribString(H5file_m,"TrackStepUnit"," ");
  H5PartWriteFileAttribString(H5file_m,"NumBunchStepUnit"," ");
  
  H5call_m = 0;

  sshift_m  = 0.0;
  RealVariable *ar = dynamic_cast<RealVariable *>(OPAL.find("SSHIFT"));
  if (ar) {
    sshift_m = ar->getReal();
    *gmsg << "SSHIFT " << sshift_m << endl;
  } 
  
}


DataSink::DataSink(int restartStep)
{
  H5PartTimer_m = IpplTimings::getTimer("H5PartTimer");
  StatMarkerTimer_m = IpplTimings::getTimer("StatMarkerTimer");

  firstWriteToStat_m = true;
  firstWriteH5part_m = false;

  string fn = OPAL.getInputFn();

  INFOMSG(fn<<endl);

  int pos=fn.find(string("."),0);
  fn.erase(pos,fn.size()-pos);

  fn += string(".h5");

#ifdef PARALLEL_IO
    H5file_m=H5PartOpenFileParallel(fn.c_str(),H5PART_APPEND,MPI_COMM_WORLD);
#else
    H5file_m=H5PartOpenFile(fn.c_str(),H5PART_APPEND);
#endif

  if(!H5file_m) {
    ERRORMSG("File open failed:  exiting!" << endl);
    exit(0);
  }

  int numStepsInFile = H5PartGetNumSteps(H5file_m);
  if(numStepsInFile < restartStep) {
    ERRORMSG("Cannot restart from a step > steps available in file: exiting!" << endl);
    exit(0);
  }
  if(numStepsInFile > restartStep) {
    ERRORMSG("Cannot restart from a step smaller than the last step in the restart file: exiting!" << endl);
    exit(0);
  }

  H5call_m = restartStep;

  sshift_m  = 0.0;
  RealVariable *ar = dynamic_cast<RealVariable *>(OPAL.find("SSHIFT"));
  if (ar) {
    sshift_m = ar->getReal();
    *gmsg << "SSHIFT " << sshift_m << endl;
  } 




}

/////////////////////////////////////////////////////////////////////////////

DataSink::~DataSink()
{
    H5PartCloseFile(H5file_m);
    Ippl::Comm->barrier();
}


void DataSink::writePhaseSpace(PartBunch &beam, Vector_t FDext[], double sposHead, double sposRef, double sposTail) {

  IpplTimings::startTimer(H5PartTimer_m); 

  double               actPos   = beam.get_sPos() +   sshift_m ;
  double               t        = beam.getT();
  Vektor< double, 3 >  rmin     = beam.get_origin();
  Vektor< double, 3 >  rmax     = beam.get_maxExtend();
  Vektor< double, 3 >  centroid = beam.get_centroid();

  size_t nTot                   = beam.getTotalNum();
  size_t nLoc                   = beam.getLocalNum();

  Vektor< double, 3 >  maxP(0.0);
  Vektor< double, 3 >  minP(0.0);

  beam.calcBeamParameters();

  Vektor<double, 3 > xsigma = beam.get_rrms();
  Vektor<double, 3 > psigma = beam.get_prms();
  Vektor<double, 3 > vareps = beam.get_emit();
  Vektor<double, 3 > geomvareps = beam.get_norm_emit();
  double meanEnergy = beam.get_meanEnergy();
  //  double *energy = beam.get_energy();


  double sigma = ((xsigma[0]*xsigma[0])+(xsigma[1]*xsigma[1])) / 
    (2.0*beam.get_gamma()*17.0e3*((geomvareps[0]*geomvareps[0]) + (geomvareps[1]*geomvareps[1]))); 


  beam.get_PBounds(minP,maxP);

  void *varray = malloc(nLoc*sizeof(double));
  double *farray = (double*)varray;
  h5part_int64_t *larray = (h5part_int64_t *)varray;
 
  /* ------------------------------------------------------------------------ 
     Get the particle decomposition from all the nodes
  */
  size_t *locN = (size_t *) malloc(Ippl::getNodes()*sizeof(size_t));
  size_t  *globN = (size_t*) malloc(Ippl::getNodes()*sizeof(size_t));
  
  for(int i=0; i<Ippl::getNodes(); i++) {
    globN[i] = locN[i]=0;
  }
  locN[Ippl::myNode()] = nLoc;
  reduce(locN, locN + Ippl::getNodes(), globN, OpAddAssign());
  
  /* ------------------------------------------------------------------------ */
  H5PartSetStep(H5file_m,H5call_m); 
  H5PartSetNumParticles(H5file_m,nLoc); 

  /* write scalar data i.e the header */
  // long long step = H5call_m;
  // H5PartWriteStepAttrib(H5file_m,"Step", H5T_NATIVE_INT64,&H5call_m,1);

  H5PartWriteStepAttrib(H5file_m,"SPOS",     H5T_NATIVE_DOUBLE,&actPos,1);
  H5PartWriteStepAttrib(H5file_m,"#sigma",     H5T_NATIVE_DOUBLE,&sigma,1);
  H5PartWriteStepAttrib(H5file_m,"RMSX",     H5T_NATIVE_DOUBLE,&xsigma,3); //sigma
  H5PartWriteStepAttrib(H5file_m,"RMSP",     H5T_NATIVE_DOUBLE,&psigma,3); //sigma
  H5PartWriteStepAttrib(H5file_m,"maxX",     H5T_NATIVE_DOUBLE,&rmax,3);
  H5PartWriteStepAttrib(H5file_m,"minX",     H5T_NATIVE_DOUBLE,&rmin,3);
  H5PartWriteStepAttrib(H5file_m,"maxP",     H5T_NATIVE_DOUBLE,&maxP,3);
  H5PartWriteStepAttrib(H5file_m,"minP",     H5T_NATIVE_DOUBLE,&minP,3);
  H5PartWriteStepAttrib(H5file_m,"centroid", H5T_NATIVE_DOUBLE,&centroid,3);
  H5PartWriteStepAttrib(H5file_m,"TIME",     H5T_NATIVE_DOUBLE,&t,1);
  H5PartWriteStepAttrib(H5file_m,"ENERGY",   H5T_NATIVE_DOUBLE,&meanEnergy,1);

  FDext[1] *= 1e-6;
  FDext[3] *= 1e-6;
  FDext[5] *= 1e-6;
  H5PartWriteStepAttrib(H5file_m,"B-head", H5T_NATIVE_DOUBLE,&FDext[0],3);
  H5PartWriteStepAttrib(H5file_m,"E-head", H5T_NATIVE_DOUBLE,&FDext[1],3);
  H5PartWriteStepAttrib(H5file_m,"B-ref",  H5T_NATIVE_DOUBLE,&FDext[2],3);
  H5PartWriteStepAttrib(H5file_m,"E-ref",  H5T_NATIVE_DOUBLE,&FDext[3],3);
  H5PartWriteStepAttrib(H5file_m,"B-tail", H5T_NATIVE_DOUBLE,&FDext[4],3);
  H5PartWriteStepAttrib(H5file_m,"E-tail", H5T_NATIVE_DOUBLE,&FDext[5],3);

  sposHead +=   sshift_m ;
  sposRef  +=   sshift_m ;
  sposTail +=   sshift_m ;

  H5PartWriteStepAttrib(H5file_m,"spos-head", H5T_NATIVE_DOUBLE, &sposHead, 1);
  H5PartWriteStepAttrib(H5file_m,"spos-ref",  H5T_NATIVE_DOUBLE, &sposRef, 1);
  H5PartWriteStepAttrib(H5file_m,"spos-tail", H5T_NATIVE_DOUBLE, &sposTail, 1);
   
  H5PartWriteStepAttrib(H5file_m,"nloc",H5T_NATIVE_INT64, globN, Ippl::getNodes());

  //not provided by PartBunch at the moment 
  //H5PartWriteStepAttrib(H5file_m,"#varepsilonr", H5T_NATIVE_DOUBLE, ,1);

  //unnormalized
  H5PartWriteStepAttrib(H5file_m,"#varepsilon", H5T_NATIVE_DOUBLE, &vareps, 3);
  //normalized
  H5PartWriteStepAttrib(H5file_m,"#varepsilon-geom", H5T_NATIVE_DOUBLE, &geomvareps, 3);

  for (size_t i=0; i<nLoc;i++)
    farray[i] =  beam.R[i](0);
  H5PartWriteDataFloat64(H5file_m,"x",farray); 
  
  for (size_t i=0; i<nLoc;i++)
    farray[i] =  beam.R[i](1);
  H5PartWriteDataFloat64(H5file_m,"y",farray); 
  
  for (size_t i=0; i<nLoc;i++)
    farray[i] =  beam.R[i](2);
  H5PartWriteDataFloat64(H5file_m,"z",farray); 
    
  for (size_t i=0; i<nLoc;i++)
    farray[i] =  beam.P[i](0);
  H5PartWriteDataFloat64(H5file_m,"px",farray); 
  
  for (size_t i=0; i<nLoc;i++)
    farray[i] =  beam.P[i](1);
  H5PartWriteDataFloat64(H5file_m,"py",farray); 
  
  for (size_t i=0; i<nLoc;i++)
    farray[i] =  beam.P[i](2);
  H5PartWriteDataFloat64(H5file_m,"pz",farray); 
  
  //for (size_t i=0; i<nLoc;i++)
  //  farray[i] =  phiz;
  //H5PartWriteDataFloat64(H5file_m,"phiz",farray); 
  
  //  H5PartWriteDataFloat64(H5file_m,"enez",energy); 
  
  for (size_t i=0; i<nLoc;i++)
    larray[i] =  beam.ID[i];
  H5PartWriteDataInt64(H5file_m,"id",larray);  

#ifdef H5BLOCKSAVE
  h5part_int64_t l[6];

  NDIndex<3> idx = getFieldLayout().getLocalNDIndex();
  NDIndex<3> elem;
  h5part_int64_t herr = H5BlockDefine3DFieldLayout (
						    H5file_m,
						    idx[0].min(), idx[0].max(),
						    idx[1].min(), idx[1].max(),
						    idx[2].min(), idx[2].max());
  if (herr < 0)
    gmsg << "H5BlockDefine3DFieldLayout err " << herr << endl;

  h5part_float64_t *data = (h5part_float64_t *) malloc ( (idx[0].max()+1)  * (idx[1].max()+1) * (idx[2].max()+1) * sizeof ( *data ) );
  
  int ii = 0;
  for (int i=idx[0].min(); i<=idx[0].max(); ++i) {
    elem[0] = Index(i,i);
    for (int j=idx[1].min(); j<=idx[1].max(); ++j) {
      elem[1] = Index(j,j);
      for (int k=idx[2].min(); k<=idx[2].max(); ++k) {
	elem[2] = Index(k,k);
	data[ii] = 0.0; // EFDMag_m.localElement(elem);
	ii++;
      }
    }
  }

  herr = H5Block3dWriteScalarField (H5file_m, "EFmag", data );
  if (herr < 0)
    gmsg << "H5Block3dWriteScalarField err " << herr << endl;

  herr = H5Block3dSetFieldSpacing (H5file_m,"EFmag",
				   (h5part_float64_t)hr_m[0],
				   (h5part_float64_t)hr_m[1],
				   (h5part_float64_t)hr_m[2]);

  if (Ippl::myNode() == 0) {
    for (h5part_int64_t p=0; p<Ippl::getNodes();p++) {
      herr = H5Block3dGetPartitionOfProc(H5file_m, p, &l[0], &l[1], &l[2], &l[3], &l[4], &l[5]);
      stringstream lstr;
      lstr << "layout" << p;
      H5BlockWriteFieldAttrib (H5file_m,"EFmag", lstr.str().c_str(), H5PART_INT64,l,6);
    }
  }
  if(data)
    free(data);
#endif

  H5Fflush(H5file_m->file,H5F_SCOPE_GLOBAL);

  H5call_m++;

  if(varray)  
    free(varray);

  IpplTimings::stopTimer(H5PartTimer_m);  
}



void DataSink::writePhaseSpace_cycl(PartBunch &beam, Vector_t FDext[]) {

  IpplTimings::startTimer(H5PartTimer_m); 

  // double               actPos   = beam.get_sPos();
  double               t        = beam.getT();
 
  Vektor< double, 3 >  rmin     = beam.get_origin();
  Vektor< double, 3 >  rmax     = beam.get_maxExtend();
  Vektor< double, 3 >  centroid = beam.get_centroid();

  size_t nTot                   = beam.getTotalNum();
  size_t nLoc                   = beam.getLocalNum();

  Vektor< double, 3 >  maxP(0.0);
  Vektor< double, 3 >  minP(0.0);

  beam.calcBeamParameters_cycl();

  Vektor<double, 3 > xsigma = beam.get_rrms();
  Vektor<double, 3 > psigma = beam.get_prms();
  Vektor<double, 3 > geomvareps = beam.get_emit();
  Vektor<double, 3 > vareps = beam.get_norm_emit();
  double meanEnergy = beam.get_meanEnergy();

  beam.get_PBounds(minP,maxP);
  
  void *varray = malloc(nLoc*sizeof(double));
  double *farray = (double*)varray;
  h5part_int64_t *larray = (h5part_int64_t *)varray;

  double  pathLength = beam.getLPath();
  h5part_int64_t trackStep =(h5part_int64_t)beam.getTrackStep();
  h5part_int64_t numBunch =(h5part_int64_t)beam.getNumBunch();

  /* ------------------------------------------------------------------------ 
     Get the particle decomposition from all the nodes
  */
  size_t *locN = (size_t *) malloc(Ippl::getNodes()*sizeof(size_t));
  size_t  *globN = (size_t*) malloc(Ippl::getNodes()*sizeof(size_t));
  
  for(int i=0; i<Ippl::getNodes(); i++) {
    globN[i] = locN[i]=0;
  }
  locN[Ippl::myNode()] = nLoc;
  reduce(locN, locN + Ippl::getNodes(), globN, OpAddAssign());
  
  /* ------------------------------------------------------------------------ */
  H5PartSetStep(H5file_m,H5call_m); 
  H5PartSetNumParticles(H5file_m,nLoc); 

  /* write scalar data i.e the header */

  // H5PartWriteStepAttrib(H5file_m,"SPOS",     H5T_NATIVE_DOUBLE,&actPos,1);
  H5PartWriteStepAttrib(H5file_m,"RMSX",     H5T_NATIVE_DOUBLE,&xsigma,3); //sigma
  H5PartWriteStepAttrib(H5file_m,"RMSP",     H5T_NATIVE_DOUBLE,&psigma,3); //sigma
  H5PartWriteStepAttrib(H5file_m,"maxX",     H5T_NATIVE_DOUBLE,&rmax,3);
  H5PartWriteStepAttrib(H5file_m,"minX",     H5T_NATIVE_DOUBLE,&rmin,3);
  H5PartWriteStepAttrib(H5file_m,"maxP",     H5T_NATIVE_DOUBLE,&maxP,3);
  H5PartWriteStepAttrib(H5file_m,"minP",     H5T_NATIVE_DOUBLE,&minP,3);
  H5PartWriteStepAttrib(H5file_m,"centroid", H5T_NATIVE_DOUBLE,&centroid,3);
  H5PartWriteStepAttrib(H5file_m,"TIME",     H5T_NATIVE_DOUBLE,&t,1);
  H5PartWriteStepAttrib(H5file_m,"ENERGY",   H5T_NATIVE_DOUBLE,&meanEnergy,1);

  // H5PartWriteStepAttrib(H5file_m,"B-head", H5T_NATIVE_DOUBLE,&FDext[0],3);
  // H5PartWriteStepAttrib(H5file_m,"E-head", H5T_NATIVE_DOUBLE,&FDext[1],3);
  H5PartWriteStepAttrib(H5file_m,"B-ref",  H5T_NATIVE_DOUBLE,&FDext[0],3);
  H5PartWriteStepAttrib(H5file_m,"E-ref",  H5T_NATIVE_DOUBLE,&FDext[1],3);
  // H5PartWriteStepAttrib(H5file_m,"B-tail", H5T_NATIVE_DOUBLE,&FDext[4],3);
  // H5PartWriteStepAttrib(H5file_m,"E-tail", H5T_NATIVE_DOUBLE,&FDext[5],3);

  //  H5PartWriteStepAttrib(H5file_m,"spos-head", H5T_NATIVE_DOUBLE, &sposHead, 1);
  //  H5PartWriteStepAttrib(H5file_m,"spos-ref",  H5T_NATIVE_DOUBLE, &sposRef, 1);
  //  H5PartWriteStepAttrib(H5file_m,"spos-tail", H5T_NATIVE_DOUBLE, &sposTail, 1);
   
  H5PartWriteStepAttrib(H5file_m,"nloc",H5T_NATIVE_INT64, globN, Ippl::getNodes());

  H5PartWriteStepAttrib(H5file_m,"LPATH",     H5T_NATIVE_DOUBLE,&pathLength,1);
  
  //not provided by PartBunch at the moment 
  //H5PartWriteStepAttrib(H5file_m,"#varepsilonr", H5T_NATIVE_DOUBLE, ,1);

  //unnormalized
  H5PartWriteStepAttrib(H5file_m,"#varepsilon", H5T_NATIVE_DOUBLE, &vareps, 3);
  //normalized
  H5PartWriteStepAttrib(H5file_m,"#varepsilon-geom", H5T_NATIVE_DOUBLE, &geomvareps, 3);

  H5PartWriteStepAttrib(H5file_m,"TrackStep",     H5T_NATIVE_INT64, &trackStep,1);

  H5PartWriteStepAttrib(H5file_m,"NumBunch",     H5T_NATIVE_INT64, &numBunch,1);
  
  for (size_t i=0; i<nLoc;i++)
    farray[i] =  beam.R[i](0);
  H5PartWriteDataFloat64(H5file_m,"x",farray); 
  
  for (size_t i=0; i<nLoc;i++)
    farray[i] =  beam.R[i](1);
  H5PartWriteDataFloat64(H5file_m,"y",farray); 
  
  for (size_t i=0; i<nLoc;i++)
    farray[i] =  beam.R[i](2);
  H5PartWriteDataFloat64(H5file_m,"z",farray); 
    
  for (size_t i=0; i<nLoc;i++)
    farray[i] =  beam.P[i](0);
  H5PartWriteDataFloat64(H5file_m,"px",farray); 
  
  for (size_t i=0; i<nLoc;i++)
    farray[i] =  beam.P[i](1);
  H5PartWriteDataFloat64(H5file_m,"py",farray); 
  
  for (size_t i=0; i<nLoc;i++)
    farray[i] =  beam.P[i](2);
  H5PartWriteDataFloat64(H5file_m,"pz",farray); 
  
  //for (size_t i=0; i<nLoc;i++)
  //  farray[i] =  phiz;
  //H5PartWriteDataFloat64(H5file_m,"phiz",farray); 
  
  //  H5PartWriteDataFloat64(H5file_m,"enez",energy); 
  
  for (size_t i=0; i<nLoc;i++)
    larray[i] =  beam.ID[i];
  H5PartWriteDataInt64(H5file_m,"id",larray);  

#ifdef H5BLOCKSAVE
  h5part_int64_t l[6];

  NDIndex<3> idx = getFieldLayout().getLocalNDIndex();
  NDIndex<3> elem;
  h5part_int64_t herr = H5BlockDefine3DFieldLayout (
						    H5file_m,
						    idx[0].min(), idx[0].max(),
						    idx[1].min(), idx[1].max(),
						    idx[2].min(), idx[2].max());
  if (herr < 0)
    gmsg << "H5BlockDefine3DFieldLayout err " << herr << endl;

  h5part_float64_t *data = (h5part_float64_t *) malloc ( (idx[0].max()+1)  * (idx[1].max()+1) * (idx[2].max()+1) * sizeof ( *data ) );
  
  int ii = 0;
  for (int i=idx[0].min(); i<=idx[0].max(); ++i) {
    elem[0] = Index(i,i);
    for (int j=idx[1].min(); j<=idx[1].max(); ++j) {
      elem[1] = Index(j,j);
      for (int k=idx[2].min(); k<=idx[2].max(); ++k) {
	elem[2] = Index(k,k);
	data[ii] = 0.0; // EFDMag_m.localElement(elem);
	ii++;
      }
    }
  }

  herr = H5Block3dWriteScalarField (H5file_m, "EFmag", data );
  if (herr < 0)
    gmsg << "H5Block3dWriteScalarField err " << herr << endl;

  herr = H5Block3dSetFieldSpacing (H5file_m,"EFmag",
				   (h5part_float64_t)hr_m[0],
				   (h5part_float64_t)hr_m[1],
				   (h5part_float64_t)hr_m[2]);

  if (Ippl::myNode() == 0) {
    for (h5part_int64_t p=0; p<Ippl::getNodes();p++) {
      herr = H5Block3dGetPartitionOfProc(H5file_m, p, &l[0], &l[1], &l[2], &l[3], &l[4], &l[5]);
      stringstream lstr;
      lstr << "layout" << p;
      H5BlockWriteFieldAttrib (H5file_m,"EFmag", lstr.str().c_str(), H5PART_INT64,l,6);
    }
  }
  if(data)
    free(data);
#endif

  H5Fflush(H5file_m->file,H5F_SCOPE_GLOBAL);

  H5call_m++;

  if(varray)  
    free(varray);

  IpplTimings::stopTimer(H5PartTimer_m);  
}

void DataSink::writeStatData(PartBunch &beam,
			    const string fname, 
			    const string el) {
  ofstream os_statData;
  ofstream momentsOf_m;

  IpplTimings::startTimer(StatMarkerTimer_m); 

  unsigned int pwi = 10;
  
  beam.gatherLoadBalanceStatistics();

  beam.calcBeamParameters();

  if(Ippl::myNode() == 0) {
    if (firstWriteToStat_m) {
      if (OPAL.inRestartRun()) {
	os_statData.open(fname.c_str(),ios::app);
	os_statData.precision(15);
	os_statData.setf(ios::scientific, ios::floatfield);
      }
      else {
	os_statData.open(fname.c_str(),ios::out);
	os_statData.precision(15);
	os_statData.setf(ios::scientific, ios::floatfield);
	// writeSDDSHeader(os_statData,OPAL.getTitle(),qq,i0,N);
      }
      firstWriteToStat_m = false;
    }
    else {
      os_statData.open(fname.c_str(),ios::app);
      os_statData.precision(15);
      os_statData.setf(ios::scientific, ios::floatfield);
    }
  }

  double spos = beam.get_sPos();
  
  if(Ippl::myNode() == 0) {
    os_statData << spos  << setw(pwi) << "\t";   // 1
    os_statData << beam.get_phase() << setw(pwi)<< "\t";    // 2
    os_statData << beam.get_gamma() << setw(pwi)<< "\t";    // 3
    
    os_statData << beam.get_rrms()(0) << setw(pwi) << "\t";             // 4
    os_statData << beam.get_rrms()(1) << setw(pwi) << "\t";             // 5
    os_statData << beam.get_rrms()(2) << setw(pwi) << "\t";             // 6
    
    os_statData << beam.get_prms()(0) << setw(pwi) << "\t";             // 7
    os_statData << beam.get_prms()(1) << setw(pwi) << "\t";             // 8
    os_statData << beam.get_prms()(2) << setw(pwi) << "\t";             // 9
    
    for (int i = 0 ; i < 3 ; ++i) 
      os_statData << beam.get_emit()(i)   << setw(pwi)<< "\t";          // 10,11,12

    for (int i = 0; i < 3; ++i)
      os_statData << beam.get_rmean()(i)  << setw(pwi)<< "\t"; // 13,14,15
    
    for (unsigned int i = 0 ; i < 3 ; ++i)                  // 16,17,18
      os_statData << beam.get_maxExtend()(i) << setw(pwi) << "\t";           

    double dummy = 0.0;
    // write out Courant Snyder parameter 
    os_statData << dummy  << setw(pwi) << "\t";           // 19
    os_statData << dummy  << setw(pwi) << "\t";           // 20
    
    os_statData << dummy << setw(pwi) << "\t";           // 21
    os_statData << dummy << setw(pwi) << "\t";           // 22
    
    // write out dispersion 
    os_statData << dummy << setw(pwi) << "\t";    // 23
    os_statData << dummy << setw(pwi) << "\t";   // 24
    os_statData << dummy << setw(pwi) << "\t";    // 25
    os_statData << dummy << setw(pwi) << "\t";   // 26

    for (int p=0; p<Ippl::getNodes(); p++)
      os_statData << beam.getLoadBalance(p)  << setw(pwi) << "\t";   // 27 ....
    
    os_statData << endl;
    os_statData.close();
  }
  
  IpplTimings::stopTimer(StatMarkerTimer_m); 
}
/***************************************************************************
 * $RCSfile: DataSink.cpp,v $   $Author: adelmann $
 * $Revision: 1.3 $   $Date: 2004/06/02 19:38:54 $
 ***************************************************************************/

