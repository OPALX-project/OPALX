// ------------------------------------------------------------------------
// $RCSfile: twostream-2.cpp,v $
// $Header: $
// ------------------------------------------------------------------------
// $Revision:  $
// $State : $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
// Description: 
//  usage: 
// ------------------------------------------------------------------------
// Class category: 
// ------------------------------------------------------------------------
//
/*
  Units: 
  
  Physics::PMASS    [kg]  == 1.6726231e-27
  Physics::EMASS    [kg]  == 9.1093897e-31
  Physics::q_e  [C]   == 1.6021773349e-19

  n:  particls to macroparticle ratio

  
  bunch->R [m]
  bunch->P [m/s] 
  bunch->E [V/m]
  dt       [s]
  bunch->Q [C ]  == (+/-) n * Physics::q_e
  bunch->M [kg ] ==       n * Physics::PMASS or Physics::EMASS


Problems:

  - Have to initialize
    eBunch->elemNo=0; 
    pBunch->elemNo=0;
  
    before calling:
     assignAccelerarorElementToPosition(pBunch,eBunch,lattice); 


ToDo:
  
   - Window for statistics


*/

#define PWI 24
#define DEGF 4

/*
  FIXME ada uggly
  Must be befor Ippl.h otherwise the 
  compiler freakes out.
  Has to do with ParticleLayoutFromGrid.h.
*/

extern "C" {
#include "cmee.h"
}

#include "Ippl.h"

#include "GTConst.hh"
#include "GTConfigure.hh"
#include "GTElemData.hh"

#include "Physics.hh"
#include "Algorithms/PartData.h"
#include "GTTimer.hh"

#include "Parser/FileStream.h"

#include <fstream>
#include <new>
#include <exception>

#include <vector>
#include <map>
#include <set> 
#include <iostream>
#include <sstream>

#include "GTChargedParticles.hh"
#include "GTDistribution.hh"

//#include "GTTimeIntegration.hh"


#ifdef GTHDF5
#include <hdf5.h>
#include "H5Part.hh"
#endif

typedef ParticleSpatialLayout<double,3>::SingleParticlePos_t Vector_t;
typedef ParticleSpatialLayout<double,3>::SingleParticlePos_t Vector3;
typedef ParticleSpatialLayout<double,3> playout_t;
typedef UniformCartesian<3,double> Mesh_t;
typedef Cell                                       Center_t;
typedef CenteredFieldLayout<3, Mesh_t, Center_t> FieldLayout_t;
typedef Field<double, 3, Mesh_t, Center_t>       Field_t;
typedef Field<Vector_t, 3, Mesh_t, Center_t>     VField_t;

using namespace std;
//using namespace SPECIS;

class H5Dump {

 public:
  
  H5Dump (string fn, ChargedParticles<playout_t> *ebunch, ChargedParticles<playout_t> *pbunch):
    H5call_m(0),
    ebunch_m(ebunch),
    pbunch_m(pbunch)
  { 
#ifdef PARALLEL_IO
    H5file_m=H5PartOpenFileParallel(fn.c_str(),H5PART_WRITE,MPI_COMM_WORLD);
#else
    H5file_m=H5PartOpenFile(fn.c_str(),H5PART_WRITE);
#endif
    if(!H5file_m) {
      ERRORMSG("File open failed:  exiting!" << endl);
      exit(0);
    }
  }
  ~H5Dump() { 

  }
  
  void closeFiles() {H5PartCloseFile(H5file_m);}

  void saveHDF5(string markerName) {
    
    Inform msg("savePhaseSpaceDataH5 ");
    Inform msg2all("savePhaseSpaceDataH5 ",INFORM_ALL_NODES);
  
    long long  Ne = ebunch_m->getLocalNum();
    long long  Np = pbunch_m->getLocalNum();
    long long   N = Ne + Np;
    
    void *varray = malloc(N*sizeof(double));
    double *farray = (double*)varray;
    long long int *larray = (long long int *)varray;
  
    /* 
       Get the particle decomposition from all the nodes
    */
    long long  *locN  = (long long*) malloc(Ippl::getNodes()*sizeof(long long));
    long long  *globN = (long long*) malloc(Ippl::getNodes()*sizeof(long long));
  
    for(int i=0; i<Ippl::getNodes(); i++) {
      globN[i] = locN[i]=0;
    }
    locN[Ippl::myNode()] = N;
    reduce(locN, locN + Ippl::getNodes(), globN, OpAddAssign());
  
    /* ------------------------------------------------------------------------ */
  
    double   spos    = pbunch_m->get_spos();
    double   periodL = pbunch_m->get_period_length();
    Vector_t origin  = ebunch_m->get_origin();

    Vector_t rmax    = ebunch_m->get_maxExtend();

    Vector_t centroid(0.0);

    Vector_t pmax,pmin;
    //pbunch_m->bounds(pbunch_m->P,pmin,pmax);

    H5PartSetStep(H5file_m,H5call_m);  
    H5PartSetNumParticles(H5file_m,N); 
  
    /* write scalar data i.e the header */
    long long step = H5call_m;

    H5PartWriteStepAttrib(H5file_m,"Step", H5T_NATIVE_INT64,&H5call_m,1);

    H5PartWriteAttrib(H5file_m,"Spos",     H5T_NATIVE_DOUBLE,&spos,1);
    H5PartWriteAttrib(H5file_m,"structLen",H5T_NATIVE_DOUBLE,&periodL,1);  
    H5PartWriteAttrib(H5file_m,"org",      H5T_NATIVE_DOUBLE,&origin,3);
    H5PartWriteAttrib(H5file_m,"maxX",     H5T_NATIVE_DOUBLE,&rmax,3);
    H5PartWriteAttrib(H5file_m,"minX",     H5T_NATIVE_DOUBLE,&origin,3);
    H5PartWriteAttrib(H5file_m,"maxP",     H5T_NATIVE_DOUBLE,&pmax,3);
    H5PartWriteAttrib(H5file_m,"minP",     H5T_NATIVE_DOUBLE,&pmin,3);
    H5PartWriteAttrib(H5file_m,"centroid", H5T_NATIVE_DOUBLE,&centroid,3);

    H5PartWriteAttrib(H5file_m,"nloc",H5T_NATIVE_INT64, globN, Ippl::getNodes());

    for (long long i=0; i<Ne;i++)
      farray[i]    =  ebunch_m->R[i](0);
    for (long long i=0; i<Np;i++)
      farray[i+Ne] =  pbunch_m->R[i](0);
    H5PartWriteDataFloat64(H5file_m,"x",farray); 

    for (long long i=0; i<Ne;i++)
      farray[i]    =  ebunch_m->R[i](1);
    for (long long i=0; i<Np;i++)
      farray[i+Ne] =  pbunch_m->R[i](1);
    H5PartWriteDataFloat64(H5file_m,"y",farray); 

    for (long long i=0; i<Ne;i++)
      farray[i]    =  ebunch_m->R[i](2);
    for (long long i=0; i<Np;i++)
      farray[i+Ne] =  pbunch_m->R[i](2);
    H5PartWriteDataFloat64(H5file_m,"z",farray); 

    for (long long i=0; i<Ne;i++)
      farray[i]    =  ebunch_m->P[i](0);
    for (long long i=0; i<Np;i++)
      farray[i+Ne] =  pbunch_m->P[i](0);
    H5PartWriteDataFloat64(H5file_m,"px",farray); 

    for (long long i=0; i<Ne;i++)
      farray[i]    =  ebunch_m->P[i](1);
    for (long long i=0; i<Np;i++)
      farray[i+Ne] =  pbunch_m->P[i](1);
    H5PartWriteDataFloat64(H5file_m,"py",farray); 

    for (long long i=0; i<Ne;i++)
      farray[i]    =  ebunch_m->P[i](2);
    for (long long i=0; i<Np;i++)
      farray[i+Ne] =  pbunch_m->P[i](2);
    H5PartWriteDataFloat64(H5file_m,"pz",farray); 

    for (long long i=0; i<Ne;i++)
      larray[i]    =   -1*(long long int)ebunch_m->ID[i];
    for (long long i=0; i<Np;i++)
      larray[i+Ne] =  pbunch_m->ID[i];
    H5PartWriteDataInt64(H5file_m,"id",larray);  
    
    H5call_m++;

    if(varray)  
      free(varray);
  }
 private:  
  H5PartFile *H5file_m;
  unsigned    H5call_m;

  ChargedParticles<playout_t> *ebunch_m;
  ChargedParticles<playout_t> *pbunch_m;
};

void setUpUnitsScalingAndOtherFunStuff(SimCfgData &simCfg,
                                       double *gamma,
                                       double *beta,
                                       double *P0,
                                       double *q,
                                       double *m,
                                       double *couplingConstant)
{
  if(simCfg.interaction==COULOMB) {
    
    // in classic 3.3 will this will be more difficult
#ifdef MAPS
    FTpsT::setGlobalTruncOrder(simCfg.order);
#endif
    /*

    Bunched Beams:
    
       Given some of the parameters: RF Frequency/rep rate f,
       charge per bunch Qtot, and the average current Iave:
       make a consistend set of data.

       Each particle has a attribute q for its charge in
       units of Coulomb. We consider only ONE bunch (so far):
                 q = Qtot/N
       where N are the number of particles.

       Given: Qtot -> Iave = Qtot * f
              Iave -> Qtot = Iave / f
    DC Beams:

    */
    if (simCfg.interaction!=COULOMB) {
      if (simCfg.I0 < 0.0) // not defined
        simCfg.I0=simCfg.Qtot*simCfg.frequ;
      else
        simCfg.Qtot=simCfg.I0/simCfg.frequ;
    }


    /*
      Physics::m_e, Physics::m_p rest mass GeV
      Ek: kinetic energy in Gev
    */

    if (simCfg.specis==PROTONS) {
      *gamma = 1.0 + (simCfg.Ek / Physics::m_p);
      *beta = sqrt(1-(1/(*gamma * *gamma)));
      *P0 = *beta * *gamma *Physics::m_p;
      *q = 1.0;
      *m = Physics::m_p;
    }
    else if (simCfg.specis==ELECTRONS) {
      *gamma = 1.0 + (simCfg.Ek / Physics::m_e);
      *beta = sqrt(1-(1/(*gamma * *gamma)));
      *P0 = *beta * *gamma *Physics::m_e;
      *q = -1.0;
      *m = Physics::m_e;
    }
    *couplingConstant = 1.0/(4.0*Physics::pi*Physics::epsilon_0);
  }
  else if (simCfg.interaction==COULOMB) {
    /*
      Astrophysical setup requested
    */
  }
  else if (simCfg.interaction==ELECTRONCLOUD) {
    /*
      MKS exept P == \beta
      like POSINST
    */
    
    // drive beam
    *gamma = 1.0 + (simCfg.Ek / Physics::m_p);
    *beta = sqrt(1-(1/(*gamma * *gamma)));
    *P0 = *beta * *gamma *Physics::m_p;
    *q = 1.0;
    *m = Physics::m_p;
    // electons
  }
}

double readInLattice(SimCfgData &simCfg, DistrData &distrData,
		     vector<PartData> &pdataVec,
		     vector< ElMap > &lattice,
		     vector<vector< ElMap >::iterator> &latItList,
		     string inpFn)
{
  Inform msg("readInLattice ");
  double gamma,beta,q,p,mass,couplingConstant;
  
  simCfg.simInputFileName = inpFn;
  FileStream::setEcho(false);
  FileStream *is;
  
  try {
    is = new FileStream(inpFn);
  } catch (...) {
    is = 0;
    cerr << "No startup file \"" << inpFn << "\" found." << endl;
  }
  
  if (!FSConfigure(is, &simCfg, &distrData))
    ERRORMSG("Error in configuration section of file \"" << inpFn << "\" found." << endl);
  
  setUpUnitsScalingAndOtherFunStuff(simCfg,&gamma,&beta,&p,&q,&mass,&couplingConstant);

  double perLength = buildLattice (is,lattice,simCfg,pdataVec[0])/simCfg.periods;
  
  PartData pd(q,mass*1.0e9,p*1.0e9);

  pdataVec[0] = pd;
  
  if (simCfg.interaction==ELECTRONCLOUD) {
    double ga = 1.0 + (simCfg.Ek / Physics::m_e);
    double be = sqrt(1-(1/(ga * ga)));
    double P0 = be * ga *Physics::m_e;
    double q = -1.0;
    double m = Physics::m_e;
    PartData pde(q,m*1.0e9,P0*1.0e9);
    pdataVec[1]= pde;
  }

  vector< ElMap >::iterator latticeIt;
  
  for(latticeIt=lattice.begin(); latticeIt!=lattice.end(); latticeIt++) {
    if (((*latticeIt).elemData.elType != STATMARKER) &&
	((*latticeIt).elemData.elType != PARTMARKER)) {
      latItList.push_back(latticeIt);
    }
  }  
  return perLength;
}    



void reInitializeTimesteps(ChargedParticles<playout_t> *eBunch, unsigned int subTSteps, 
			 ChargedParticles<playout_t> *pBunch)
{ 
  pBunch->setdtleft(1);
  pBunch->resetStatCounters();
  
  eBunch->setdtleft(subTSteps);
  eBunch->resetStatCounters(); 
}

void assignAccelerarorElementToPosition(ChargedParticles<playout_t> *eBunch, 
					ChargedParticles<playout_t> *pBunch,
					vector< ElMap > &lattice)
{ 
  vector< ElMap >::iterator latticeIt;
  int pos = 0;
  
  for(latticeIt=lattice.begin(); latticeIt!=lattice.end(); latticeIt++) {
    if (((*latticeIt).elemData.elType != STATMARKER) &&
	((*latticeIt).elemData.elType != PARTMARKER)) {
      double beginpos = (*latticeIt).elemData.entrance;
      double endpos   =  beginpos +(*latticeIt).elemData.length; 
      INFOMSG("pos= " << pos << " beginpos= " << beginpos << " endpos= " << endpos << endl);
      for (unsigned i=0; i<pBunch->getLocalNum(); i++) {
	if ((pBunch->R[i](2)>=beginpos) && (pBunch->R[i](2)<endpos)) {
	  if (pos > pBunch->getnumberOfLatticeElements())
	    ERRORMSG("pos " << pos << " > numberOfLatticeElements_m " << pBunch->getnumberOfLatticeElements() << endl);
	  pBunch->elemNo=pos;
	  pBunch->stat = ALIVE;
	}
      }
      pos++;
    }
  }
}

void printStatus(bool full, 
		 ChargedParticles<playout_t> *eBunch, 
		 ChargedParticles<playout_t> *pBunch,
		 SimCfgData &simCfg,
		 DistrData &distrData,
		 vector<PartData> &pdataVec,
		 vector< ElMap > &lattice,
		 string fn1, string fn2)
{
  Inform m("PARSEC ");  
  if (full) {
    m << endl << endl;
    m << "---------------------- PARSEC -----------------------------------------------------" << endl;
    m << "Lattice input file " << fn1 << " secondary emission configuration " << fn2 << endl;
    m << "Number of processors " << Ippl::getNodes() << " serial dimension is z " << endl;
    m << "Base filename for data save " << simCfg.diagfile << endl << endl; 
    m << "Drive beam gamma= " << pdataVec[0].getGamma() << endl;
  
    m << "Distribution informations " << endl;
    m << "------------------------- " << endl;      
    m << " Protons Gauss:  sigx= " << distrData.sigx << " [m]\t sigpx= " <<  distrData.sigpx << " [m/s]\t xpx-corel= " << distrData.corxpx << " [1]" << endl; 
    m << "                 sigy= " << distrData.sigy << " [m]\t sigpy= " <<  distrData.sigpy << " [m/s]\t ypy-corel= " << distrData.corypy << " [1]" << endl; 
    m << "                 sigz= " << distrData.sigt << " [m]\t sigpz= " <<  distrData.sigpt << " [m/s]\t zpz-corel= " << distrData.cortpt << " [1]" << endl; 
    m << " Electrons uniform: max(|E|) = 100 [eV] " << endl << endl; 

    m << "Lattice informations " << endl;
    m << "-------------------- " << endl;
        
    vector< ElMap >::iterator latticeIt;
    unsigned int latItListLen = 0;
    
    for(latticeIt=lattice.begin(); latticeIt!=lattice.end(); latticeIt++) {
      if (((*latticeIt).elemData.elType != STATMARKER) &&
	  ((*latticeIt).elemData.elType != PARTMARKER)) {
	latItListLen++;
	m << (*latticeIt).elemData.convStr();
	m << " entrance at (" << (*latticeIt).elemData.entrance;
	m << " exit at " << (*latticeIt).elemData.entrance+(*latticeIt).elemData.length << "(" 
	  << " a= "  << (*latticeIt).elemData.a 
	  << " b= "  << (*latticeIt).elemData.b 
	  << " L= "  << (*latticeIt).elemData.length 
	  << " kx= " << (*latticeIt).elemData.K0 
	  << " ky= " << (*latticeIt).elemData.K1
	  << " kz= " << (*latticeIt).elemData.K2 << endl << endl;	
      }
    }
    
    m << "Lattice indirection list length = " << latItListLen << endl << endl;
     
    m << "Parameters for electron proton ratio " << endl
      << "------------------------------------ " << endl
      << " Proton line density= " << simCfg.lambdaP << "\t electron line density= " << simCfg.neutralfac*pow(distrData.sigy/distrData.sigx,2.0)*simCfg.lambdaP << endl
      << " Proton neutralisation factor= " << simCfg.neutralfac << endl
      << " Np= " << pBunch->getTotalNum() << " \t Qp= " << pBunch->getQTot() << endl
      << " Ne= " << eBunch->getTotalNum() << " \t Qe= " << eBunch->getQTot() << endl << endl;
  
    m << "Parameters for time integration " << endl
      << "------------------------------- " << endl
      << " Tperiod= " << pBunch->get_period_length()/(pBunch->getBeta()*Physics::c) << " [s]\t NdT= " << simCfg.ndt << " [1]\t dTp= " << simCfg.dt 
      << " [s]\t Nperiods= " << simCfg.periods << endl
      << " Sub-timesteps for electrons M= " << simCfg.subTSteps << " [1]\t dTe= " << simCfg.dt/simCfg.subTSteps << " [s]" << endl << endl;

    m << "Parameters from the electron and proton bunches " << endl
      << "----------------------------------------------- " << endl
      << " Electrons maxX= " << eBunch->get_maxExtend() << " [m]\t minX= " <<  eBunch->get_origin() << " [m]\t"
      << "average Ekin= " << eBunch->get_meanEkineV() - (Physics::m_e*1E9)   << " [eV]" << endl
      << " Protons maxX= " << pBunch->get_maxExtend() << " [m]\t minX= " <<  pBunch->get_origin() << " [m]\t"
      << "average Etot= " << pBunch->get_meanEkineV()/1E9 << " [GeV]" << endl << endl;
  }
  else {

 
  }
}

void renorm_ellip(Vector_t &x, double a, double b, Vector_t &xn, Vector_t &tg)
/* 
   Given a point Q(x0,y0) near the ellipse (x/a)**2+(y/b)**2=1, finds the 
   coordinates (x,y) of the point P on the ellipse that is closest to Q. 
   The calculation is first-order in the distance from P to Q, so 
   it is valid only is Q is close enough to the ellipse (but Q is allowed 
   to be either inside or outside the ellipse).
   Then computes the inward normal and CCW tangent unit vectors
   Then substitute output coordinates in place of the old ones
   Then, renormalize output (x,y) so that it is just inside the ellipse
*/
{	

  const double pi=acos(-1.0);
  const double twopi=2.0*pi;
  
  // for a chamber of radius=2 cm, the value acc=5e-5 will reposition the 
  // electron 1 micron inside the surface
  double acc=5e-5;
  double renorm=1.0-acc;
  
  double phi;
  
  // first, determine angle phi such that x=a*r*cos(phi), y=b*r*sin(phi),
  // and r=sqrt((x/a)**2+(y/b)**2));
  // special cases: four corners (these are extremely infrequent cases)

  if(x[0]>0 && x[1]==0) 
    phi=0;
  else if(x[0]==0 && x[1]>0)
    phi=pi/2;
  else if(x[0]< 0 && x[1]==0)
    phi=pi;
  else if(x[0]==0 && x[0]<0)
    phi=3*pi/2;
  
  // generic cases:
  double rho=abs((a*x[1])/(b*x[0]));
  
  if(x[0]>0 && x[1]>0)	        // 1st quadrant
    phi=atan(rho);
  else if(x[0]<0 && x[1]>0)	// 2nd quadrant
    phi=pi-atan(rho);
  else if(x[0]<0 && x[1]<0)	// 3rd quadrant
    phi=pi+atan(rho);
  else if(x[0]>0 && x[1]<0)	// 4th quadrant
    phi=twopi-atan(rho);
  
  
  if(a==b) {
    x[0]=a*cos(phi);
    x[1]=b*sin(phi);
  }
  else if(a<=b) {
    double rsq=(x[0]*x[0])/(a*a)+(x[1]*x[1])/(b*b);
    double eps=(rsq-1)/(sqrt(rsq)+1);
    double gsq=(a+b)*(a-b);
    double alpha=-eps*gsq*x[0]*x[1]*a*b/(a*a*a*a*x[1]*x[1]+b*b*b*b*x[0]*x[0]);
    x[0]=a*cos(phi+alpha);
    x[1]=b*sin(phi+alpha);
  }
  xn[0]=-x[0]/(a*a);
  xn[1]=-x[1]/(b*b);
  double xnorm=sqrt(xn[0]*xn[0] + xn[1]*xn[1]);
  xn[0]=xn[0]/xnorm;
  xn[1]=xn[1]/xnorm;	// inward normal unit vector
  
  tg[0]=xn[0];		// CCW tangent unit vector
  tg[1]=-xn[1];
  // if you do not renormalize, round-off will put some eletrons slightly
  // outside the chamber, others slightly inside
  x[0]=renorm*x[0];
  x[1]=renorm*x[1];
}


inline void secElec(ChargedParticles<playout_t> *bunch, vector<vector< ElMap >::iterator> latItList)
{
  Inform m ("secElec ");
  // m << "inc x,y= " << x << " " << y << " Einc= " << inc_energy << " Rrenor " <<  bunch->R[i] << " ns= " << ns[0] << endl;

  unsigned nei        = bunch->getTotalNum();
  unsigned neAtWall   = 0;
  double dr           = 1.0;                     // radial distance to chamber
  double dt           = bunch->getdt() / 2.0 ;   // initial dt

  double xf,yf;
  double rGeom;
  double px;
  double py;
  double  x;
  double  y;
  /*
    cemee stuff
  */
  int ns[1];
  double bn[MAXSEC];
  double bt[MAXSEC];
  double bz[MAXSEC];

  double maxEsecElec =  0.0;
  double eSecElec    =  0.0;
  unsigned nsTot     = 0;
  
  /* The material the electron hits is described 
   * by an integer (1=copper, 2=stainless steel)
   */
  int mat_number = 1;

  for(unsigned i=0; i < bunch->getLocalNum(); i++) {          
    if(bunch->stat[i] == ATTHEWALL) {
      neAtWall++;
      rGeom = (*latItList[bunch->elemNo[i]]).elemData.a;
      px = bunch->P[i](0);
      py = bunch->P[i](1);
      x = bunch->R[i](0);
      y = bunch->R[i](1);
      
      /*
	Next calculate incident electron energy in eV  and 
	angle of incidence expressed as cosine of the angle measured 
	with angle zero as normal (costheta = 1. at normal incidence) 
      */
      
      double betainc2   = dot(bunch->P[i],bunch->P[i])/(Physics::c*Physics::c);
      double gammainc   = 1.0/sqrt(1.0 - betainc2);
      double inc_energy = 1E9*Physics::m_e*betainc2*gammainc*gammainc/(1.0+gammainc);
      
      double bnin = x*px/Physics::c + y*py/Physics::c;
      
      Vector_t xn; Vector_t tg;
      renorm_ellip(bunch->R[i], rGeom, rGeom, xn, tg);
      
      double costheta = (bunch->R[i](0)*xn[0]+bunch->R[i](1)*xn[1])/sqrt(betainc2);
      
      // nsec fills the return arrays ns, bn, bt and bz
      nsec(inc_energy, costheta, mat_number, ns, bn, bt, bz);
      nsTot += ns[0];

      unsigned bottom = bunch->getLocalNum();
      bunch->create(ns[0]);
      for(unsigned int k=0;k<ns[0];k++) {
	bunch->R[bottom+k]      = bunch->R[i];
	bunch->P[bottom+k]      = Physics::c * Vector_t(bn[k],bt[k],bz[k]);
	bunch->dtleft[bottom+k]  = bunch->dtleft[i];
	bunch->ndt[bottom+k]     = bunch->ndt[i];
	bunch->elemNo[bottom+k] = bunch->elemNo[i];
	bunch->stat[bottom+k]   = ALIVE;
	
	maxEsecElec = max(maxEsecElec,inc_energy);
	eSecElec   += inc_energy;
      }
      bunch->destroy(1,i);
    }
  }
  bunch->boundp();
  reduce(nsTot,nsTot,OpAddAssign());
  reduce(neAtWall,neAtWall,OpAddAssign());
  reduce(eSecElec,eSecElec,OpAddAssign());
  reduce(maxEsecElec,maxEsecElec,OpMaxAssign());

  bunch->storeEStat(neAtWall,eSecElec,maxEsecElec,nsTot); // will be saved with next writeStatistics
  if (neAtWall > 0) 
    m << "tacr= " << bunch->getTime() 
      << " [s] nE@Wall= " << neAtWall << "\t   maxEsecElec= " << maxEsecElec 
      << " [eV] \tNetot= " << bunch->getTotalNum() 
      << " Necre= " << nsTot << endl;
}

inline void constFoc(ChargedParticles<playout_t> *bunch, vector<vector< ElMap >::iterator> latItList)
{
  
  /*
    Units: 

     Physics::PMASS    [kg]  == 1.6726231e-27
     Physics::EMASS    [kg]  == 9.1093897e-31
     Physics::q_e  [C]   == 1.6021773349e-19

     n:  particls to macroparticle ratio


     bunch->R [m]
     bunch->P [m/s] 
     bunch->E [V/m]
     dt       [s]
     bunch->Q [C ]  == (+/-) n * Physics::q_e
     bunch->M [kg ] ==       n * Physics::PMASS or Physics::EMASS

  */
  const double Tper = bunch->get_period_length()/(bunch->getBeta()*Physics::c);
  const double omega0 = 2.0*Physics::pi/Tper;

  int latElem = -1;
  double Qx, Qy, Qz;
  double omegax,omegay,omegaz;

  for(unsigned i=0; i < bunch->getLocalNum(); i++) {          
    
    if (bunch->elemNo[i] != latElem) {
      latElem = bunch->elemNo[i];
      Qx = (*latItList[latElem]).elemData.K0; 
      Qy = (*latItList[latElem]).elemData.K1; 
      Qz = (*latItList[latElem]).elemData.K2; 
      omegax = Qx*omega0;
      omegay = Qy*omega0;
      omegaz = Qz*omega0;
    }
    const double px = bunch->P[i](0);
    const double py = bunch->P[i](1);
    const double pz = bunch->P[i](2);
    const double  x = bunch->R[i](0);
    const double  y = bunch->R[i](1);
    const double  z = bunch->R[i](2);

    const double dt = bunch->getdt() * bunch->ndt[i];
    
    bunch->R[i](0) =         x*cos(omegax*dt) + 1.0/(omegax)*px*sin(omegax*dt); 
    bunch->P[i](0) = -omegax*x*sin(omegax*dt) +              px*cos(omegax*dt);
    
    bunch->R[i](1) =         y*cos(omegay*dt) + 1.0/(omegay)*py*sin(omegay*dt);  
    bunch->P[i](1) = -omegay*y*sin(omegay*dt) +              py*cos(omegay*dt);      
    
    bunch->R[i](2) =         z*cos(omegaz*dt) + 1.0/(omegaz)*pz*sin(omegaz*dt);  
    bunch->P[i](2) = -omegaz*z*sin(omegaz*dt) +              pz*cos(omegaz*dt);      
 
    /*
      checkBC == true
      particle hits the wall and we have
      
      status == ATTHEWALL
      R and P old values,
      we will determine LATER the exact location
      and do the secondary emission because of that
      we copy back the phase space coordinates.
      
    */
    
    if (bunch->checkBC(latItList,i)) {
      bunch->R[i](0) = x;
      bunch->P[i](0) = px;
      bunch->R[i](1) = y;
      bunch->P[i](1) = py;
      bunch->R[i](2) = z;
      bunch->P[i](2) = pz;
    }
    else {
      bunch->dtleft[i]--;
      if (bunch->dtleft[i] == 0)
	bunch->stat = DONE;
    }
  }
  bunch->boundp();  
}



int main(int argc, char ** argv)
{
  Ippl ippl(argc, argv);
  Inform msg("twostream-2 ");  
  Inform msg2all("twostream-2 ",INFORM_ALL_NODES);  
  IpplRandom.SetSeed(234244113131719);

  SimCfgData simCfg;
  ElemData elemData; 
  DistrData distrData;
  vector< ElMap > lattice;
  vector<PartData> pdataVec(2);

  vector<vector< ElMap >::iterator> latItList;

  double perLength = readInLattice(simCfg, distrData, pdataVec, lattice, latItList, string(argv[1]));
  
  double gamma = pdataVec[0].getGamma();

  double a =  (*latItList[0]).elemData.a; 
  double b =  (*latItList[0]).elemData.b; 
  double Lg = gamma*perLength;
  
  double nRealProtons   = simCfg.lambdaP*perLength;
  double protonMacroPartRatio = nRealProtons/simCfg.nPInit;
  
  double lambdaE = simCfg.neutralfac*pow(distrData.sigy/distrData.sigx,2.0)*simCfg.lambdaP;
  double nRealElectrons = lambdaE*perLength;
  double electronMacroPartRatio = nRealElectrons/simCfg.nEInit;
              
  /* ------------------------------------------------------------------------------------------------------------------ */
  /* Create a bunch                                                                                                     */
  /* ------------------------------------------------------------------------------------------------------------------ */

  e_dim_tag decomp[3];
  int serialDim = 2;
  Vektor<int,3> nr(64);
  
  Mesh_t *mesh;
  FieldLayout_t *FL;
  ChargedParticles<playout_t> *P;
  
  NDIndex<3> domain;
  for(int i=0; i<3; i++)
    domain[i] = domain[i] = Index(nr[i] + 1);
  
  for (int d=0; d < 3; ++d)
    decomp[d] = (d == serialDim) ? SERIAL : PARALLEL;
  
  // create mesh and layout objects for this problem domain
  mesh          = new Mesh_t(domain);
  FL            = new FieldLayout_t(*mesh, decomp);
  playout_t* PL = new playout_t(*FL, *mesh);

  ChargedParticles<playout_t> *eBunch = new ChargedParticles<playout_t>(PL,decomp,
									electronMacroPartRatio,
									-Physics::q_e,Physics::EMASS,ELECTRONS,
									gamma,
									perLength,1);

  ChargedParticles<playout_t> *pBunch = new ChargedParticles<playout_t>(PL,decomp,
									protonMacroPartRatio,
									Physics::q_e,Physics::PMASS,PROTONS,
									gamma,
									perLength,1);
  
  /* ------------------------------------------------------------------------------------------------------------------ */
  /* Create Distributions                                                                                               */
  /* ------------------------------------------------------------------------------------------------------------------ */

  Distribution *distr = new Distribution(&distrData,perLength);

  /*
    initial electron background
  */
  
  double maxE = 270.0; // maximal energy of generated electrons [eV]
  distr->createUniformElectronBackground(eBunch,simCfg.nEInit, a, perLength, maxE);
  
  /*
    drive beam (protons) 
  */
  Vector_t xmean   = Vector_t(0.0,0.0,0.5);
  Vector_t pmean (0.0);
  Vector_t xstddev = Vector_t(distrData.sigx, distrData.sigy, distrData.sigt);
  Vector_t pstddev = Vector_t(distrData.sigpx,distrData.sigpy,distrData.sigpt);
  Vector_t angle (0.0);
  distr->generateGaussian(pBunch,simCfg.nPInit, xmean, pmean, xstddev, pstddev, angle);

  /* ------------------------------------------------------------------------------------------------------------------ */
  /* Save initial conditions                                                                                            */
  /* ------------------------------------------------------------------------------------------------------------------ */

  pBunch->openFiles((simCfg.diagfile+string(".p")),simCfg.title);
  pBunch->writeStatistics();

  eBunch->openFiles((simCfg.diagfile+string(".e")),simCfg.title);
  eBunch->writeStatistics();
  
  pBunch->writePhaseSpaceSDDS((simCfg.diagfile+string(".p")));
  eBunch->writePhaseSpaceSDDS((simCfg.diagfile+string(".e")));

  H5Dump *h5dump = new H5Dump((simCfg.diagfile+string(".h5")),eBunch,pBunch);
  h5dump->saveHDF5(string("H5TEST"));

  /* ------------------------------------------------------------------------------------------------------------------ */
  /* Prepare particle integration                                                                                       */
  /* ------------------------------------------------------------------------------------------------------------------ */

  /*
    setdT sets the macro intervall for the timeintegration

    The dT then is subdivided by M which defines a micro intervall
    to respect the different timescales of the particles 

    After each dT a space charge calculation takes place.
  */


  double Tper = perLength/(pBunch->getBeta()*Physics::c);
  
  simCfg.dt = Tper/simCfg.ndt;
  
  pBunch->setdt(simCfg.dt); 
  pBunch->setdtleft(1);
  pBunch->setndt(1);
  pBunch->resetStatCounters();
  
  eBunch->setdt(simCfg.dt/simCfg.subTSteps);
  eBunch->setdtleft(simCfg.subTSteps);
  eBunch->setndt(1);
  eBunch->resetStatCounters(); 
  
  assignAccelerarorElementToPosition(pBunch,eBunch,lattice); 

  
  bool printStatusFull = true;
  printStatus(printStatusFull,eBunch,pBunch,simCfg,distrData,pdataVec,lattice,string(argv[1]),string(argv[2]));
      
  /* ------------------------------------------------------------------------------------------------------------------ */
  /* Do particle integration                                                                                            */
  /* ------------------------------------------------------------------------------------------------------------------ */

  while (pBunch->getTime() < Tper*simCfg.periods) {

    reInitializeTimesteps(eBunch,simCfg.subTSteps,pBunch); 
    
    // kick only protons
    constFoc(pBunch,latItList); 
    
    // kick all electrons until all have the same time 
    unsigned nLeft = 0;
    do {
      constFoc(eBunch,latItList); 
      secElec(eBunch,latItList); 
      nLeft = eBunch->allDone();
    } while (nLeft != 0);

    eBunch->incTime(eBunch->getdt());
    pBunch->incTime(pBunch->getdt());
    pBunch->writeStatistics();
    eBunch->writeStatistics();

    pBunch->writePhaseSpaceSDDS((simCfg.diagfile+string(".p")));
    eBunch->writePhaseSpaceSDDS((simCfg.diagfile+string(".e")));
    h5dump->saveHDF5(string("H5TEST"));
    
    /*
      First calculate the electric field of the possible relativistic proton
      beam:
	  
      x' = L(x)
      \phi' = \nabla^2 \rho'
      E' = \nable\phi'
      E = L^{-1}(E')
	  
      first solve for the relativistic protons
      then solve for the electrons superposition of th
      fields done in gradPhi(bunch)
    */
      
    // kick electrons with leap frog second part and 
    // apply space charge kick to protons 
    //kickEP(bunch,latItList,perLength,simCfg.dt); 
	
    // kick only protons

    reInitializeTimesteps(eBunch,simCfg.subTSteps,pBunch); 

    // kick only protons
    constFoc(pBunch,latItList); 
    
    // kick all electrons until all have the same time 

    do {
      constFoc(eBunch,latItList); 
      secElec(eBunch,latItList); 
    } while (eBunch->allDone() != 0);
  
    eBunch->incTime(eBunch->getdt());
    pBunch->incTime(pBunch->getdt());
    pBunch->writeStatistics();
    eBunch->writeStatistics();

    msg << "@ sact= " << pBunch->get_spos() << endl;
  }
  
  h5dump->closeFiles();
  pBunch->closeFiles();
  msg << "all done ....  " << endl;
}
