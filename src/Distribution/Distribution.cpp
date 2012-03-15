// ------------------------------------------------------------------------
// $RCSfile: Distribution.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.3.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Distribution
//   The class for the OPAL Distribution command.
//
// ------------------------------------------------------------------------

#include "Distribution/Distribution.h"
#include "BasicActions/Option.h"
#include "AbstractObjects/OpalData.h"
#include <fstream>
#include <string>
#include <cmath>
#include <cfloat>
#include <vector>
#include <iostream>  // Neeeded for stream I/O
#include <fstream>   // Needed for file I/O
#include <iomanip>   // Needed for I/O manipulators


#include "Distribution/Bins.h"

using namespace Expressions;
using namespace Physics;
using namespace Attributes;

// Class Distribution
// ------------------------------------------------------------------------

// The attributes of class Distribution.

namespace {
  enum {
    // DESCRIPTION OF THE DISTRIBUTION:
    DISTRIBUTION,  
    FNAME,
    XMULT,
    YMULT,
    TMULT,
    PXMULT,
    PYMULT,
    PTMULT,
    BETAX,
    BETAY,
    ALPHAX,
    ALPHAY,
    MX,
    MY,
    MT,
    DX,
    DDX,
    DY,
    DDY,
    SIGMAX,
    SIGMAY,
    SIGMAT,
    SIGMAPX,
    SIGMAPY,
    SIGMAPT,
    CORRX,
    CORRY,
    CORRT,
    TEMISSION,
    NBIN,
    DEBIN,
    SIZE
  };
}

Distribution::Distribution():
  Definition(SIZE,"DISTRIBUTION","The DISTRIBUTION statement defines data for the 6D particle distr.")
{
  itsAttr[DISTRIBUTION] = makeString("DISTRIBUTION", "Distribution type: GAUSS, BINOMIAL, FROMFILE, GUNGAUSS, UNIFORMXYZ ", "GAUSS");

  itsAttr[FNAME] = makeString("FNAME", "File for read in 6D particle coordinates");

  itsAttr[XMULT] = makeReal("XMULT", "Multiplier for X",1.0);
  itsAttr[YMULT] = makeReal("YMULT", "Multiplier for Y",1.0);
  itsAttr[TMULT] = makeReal("TMULT", "Multiplier for T",1.0);

  itsAttr[PXMULT]= makeReal("PXMULT", "Multiplier for PX",1.0);
  itsAttr[PYMULT]= makeReal("PYMULT", "Multiplier for PY",1.0);
  itsAttr[PTMULT]= makeReal("PTMULT", "Multiplier for PT",1.0);
    
  itsAttr[ALPHAX]= makeReal("ALPHAX", "Courant Sinder parameter",1.0);
  itsAttr[ALPHAY]= makeReal("ALPHAY", "Courant Sinder parameter",1.0);

  itsAttr[BETAX] = makeReal("BETAX", "Courant Sinder parameter",-1.0);   
  itsAttr[BETAY] = makeReal("BETAY", "Courant Sinder parameter",1.0);

  itsAttr[MX]    = makeReal("MX", "Defines the distribution in x, 0+eps .. inf",1.0);
  itsAttr[MY]    = makeReal("MY", "Defines the distribution in y, 0+eps .. inf",1.0);
  itsAttr[MT]    = makeReal("MT", "Defines the distribution in t, 0+eps .. inf",1.0);

  itsAttr[DX]    = makeReal("DX", "Dispersion in x",0.0);      
  itsAttr[DDX]   = makeReal("DDX", "First derivative of Dx",0.0);
 
  itsAttr[DY]    = makeReal("DY", "Dispersion in x",0.0);      
  itsAttr[DDY]   = makeReal("DDY", "First derivative of Dy",0.0);
  
  itsAttr[SIGMAX] = makeReal("SIGMAX", "SIGMAx",1.0e-2);   
  itsAttr[SIGMAY] = makeReal("SIGMAY", "SIGNAy",1.0e-2);
  itsAttr[SIGMAT] = makeReal("SIGMAT", "SIGMAt",1.0e-2);

  itsAttr[SIGMAPX] = makeReal("SIGMAPX", "SIGMApx",0.0);   
  itsAttr[SIGMAPY] = makeReal("SIGMAPY", "SIGNApy",0.0);
  itsAttr[SIGMAPT] = makeReal("SIGMAPT", "SIGMApt",0.0);

  itsAttr[CORRX] = makeReal("CORRX", "CORRx",-0.5);
  itsAttr[CORRY] = makeReal("CORRY", "CORRy", 0.5);
  itsAttr[CORRT] = makeReal("CORRT", "CORRt", 0.0);

  itsAttr[TEMISSION] = makeReal("TEMISSION", "Time in seconds in which we have emission",  0.0);
  itsAttr[NBIN]      = makeReal("NBIN", "In case of emission how many bins should we use", 0.0);
  itsAttr[DEBIN]     = makeReal("DEBIN", "Energy band for a bin in keV, defines the rebinning", 1000000.0);

  // Set up default beam.
  Distribution *defDistribution = clone("UNNAMED_Distribution");
  defDistribution->builtin = true;

  try {
    defDistribution->update();
    OPAL.define(defDistribution);
  } catch (...) {
    delete defDistribution;
  }
  pbin_m = 0;
}

Distribution::Distribution(const string &name, Distribution *parent):
  Definition(name, parent),
  reference(parent->reference)
{}

Distribution::~Distribution()
{
  if (pbin_m)
    delete pbin_m;
}

bool Distribution::canReplaceBy(Object *object)
{
  // Can replace only by another Distribution.
  return dynamic_cast<Distribution *>(object) != 0;
}


Distribution *Distribution::clone(const string &name)
{
  return new Distribution(name, this);
}


void Distribution::execute()
{
  update();
}

void Distribution::createBinom(Vector_t emit, Vector_t alpha, Vector_t beta, Vector_t gamma, 
			       Vector_t bincoef, bool gauss, PartBunch &beam, size_t particles,
			       bool isBinned) 
{
    const double two_pi = 2.0 * 4.0*atan(1.0);

    Vector_t M;
    Vector_t PM;
    Vector_t COSCHI;
    Vector_t SINCHI;
    Vector_t CHI;
    Vector_t AMI;
    Vector_t L;
    Vector_t PL;

    size_t pc = 0;
    size_t count = 0;      

    for (int i=0;i<3;i++) {
      gamma[i] *= 4.0;
      beta[i] *= 4.0;
      M[i]       =  sqrt(emit[i]*beta[i]);
      PM[i]      =  sqrt(emit[i]*gamma[i]);
      COSCHI[i]  =  sqrt(1.0/(1.0+alpha[i]*alpha[i]));
      SINCHI[i]  = -alpha[i]*COSCHI[i];
      CHI[i]     =  atan2(SINCHI[i],COSCHI[i]);
      AMI[i]     =  1.0/bincoef[i];
      L[i]       =  sqrt((bincoef[i]+1.0)/2.0)*M[i];
      PL[i]      =  sqrt((bincoef[i]+1.0)/2.0)*PM[i];
    }
    
    Vector_t x;
    Vector_t p;
    
    if (!gauss) {
      for (size_t n=0; n < particles; ++n) {
	double S1,S2;
	double A,AL,U,V;
	for (int i=0;i<3;i++) {            
	  S1=IpplRandom(); S2=IpplRandom();
	  A=sqrt(1.0 - pow(S1,AMI[i]));
	  AL=two_pi*S2;
	  U=A*cos(AL);
	  V=A*sin(AL);
	  x[i] = L[i]*U;                    
	  p[i] = PL[i]*(U*SINCHI[i]+V*COSCHI[i]);
	}
	if (pc == Ippl::myNode()) {
	  if (isBinned) {
	    vector<double> tmp;
	    tmp.push_back(x[0]);
	    tmp.push_back(x[1]);
	    tmp.push_back(x[2]);
	    tmp.push_back(0.0);
	    tmp.push_back(0.0);
	    tmp.push_back(p[2]);
	    tmp.push_back(0);
	    pbin_m->fill(tmp);
	  }
	  else {
	    beam.create(1);
	    beam.R[count] = x;
	    beam.P[count] = p;
	  }
	  count++; 
	}
	pc++;
	if (pc == Ippl::getNodes())
	  pc = 0;
      }
    }
    else {
      for (size_t n=0; n < particles; ++n) {
	double S1,S2;
	double A,AL,U,V;
	for (int i=0;i<3;i++) {       
	  S1=IpplRandom(); S2=IpplRandom();
	  A=sqrt(2.0)/2.0 * sqrt(-log(S1));
	  AL=two_pi*S2;
	  U=A*cos(AL);
	  V=A*sin(AL);
	  x[i]=M[i]*U;
	  p[i]=PM[i]*(U*SINCHI[i]+V*COSCHI[i]);
	}
	if (pc == Ippl::myNode()) {
	  if (isBinned) {
	    vector<double> tmp;
	    tmp.push_back(x[0]);
	    tmp.push_back(x[1]);
	    tmp.push_back(x[2]);
	    tmp.push_back(0.0);
	    tmp.push_back(0.0);
	    tmp.push_back(p[2]);
	    tmp.push_back(0);
	    pbin_m->fill(tmp);
	  }
	  else {
	    beam.create(1);
	    beam.R[count] = x;
	    beam.P[count] = p;
	  }
	  count++; 
	}
	pc++;
	if (pc == Ippl::getNodes())
	  pc = 0;
      }
    }
    if (isBinned) {
      pbin_m->sortArray();  
      // now copy this over to the bunch
      // so that we can emmit the particles  
      beam.setPBins(pbin_m);
    }
}

void Distribution::create(PartBunch &beam, size_t Np) {

  Inform msg("Distribution::create ");

  int ebins = (int) Attributes::getReal(itsAttr[NBIN]);
  bool isBinned = (ebins > 0);

  if (isBinned)
    pbin_m = new PartBins((int) Attributes::getReal(itsAttr[NBIN]));
  else
    pbin_m = NULL;
  
  beam.setTEmission(Attributes::getReal(itsAttr[TEMISSION]));
  beam.setNumBunch(1);

  if (Attributes::getString(itsAttr[DISTRIBUTION]) == "GUNGAUSS") {
    binnDistribution(beam, Np);
  }
  else if (Attributes::getString(itsAttr[DISTRIBUTION]) == "BINOMIAL") {

    Vector_t corr(Attributes::getReal(itsAttr[CORRX]),
		  Attributes::getReal(itsAttr[CORRY]),
		  Attributes::getReal(itsAttr[CORRT]));  
    
    Vector_t sigX(Attributes::getReal(itsAttr[SIGMAX]),
		  Attributes::getReal(itsAttr[SIGMAY]),
		  Attributes::getReal(itsAttr[SIGMAT]));  
    
    Vector_t sigPX(Attributes::getReal(itsAttr[SIGMAPX]),
		   Attributes::getReal(itsAttr[SIGMAPY]),
		   Attributes::getReal(itsAttr[SIGMAPT]));  

    Vector_t bincoef(Attributes::getReal(itsAttr[MX]),
		     Attributes::getReal(itsAttr[MY]),
		     Attributes::getReal(itsAttr[MT]));  
    
    Vector_t emit;
    Vector_t alpha;
    Vector_t beta;
    Vector_t gamma; 
    
    for(int j=0; j<3; j++) {
      double chi = asin(corr[j]);    
      emit[j] = sigX[j]*sigPX[j]*cos(chi);
    }
    for(int j=0; j<3; j++) {
      beta[j]  = sigX[j]*sigX[j]/emit[j];
      gamma[j] = sigPX[j]*sigPX[j]/emit[j];
      alpha[j] = -corr[j]*sqrt(beta[j]*abs(gamma[j]));
    }
    
    /// if M is larger than 10000 we switch to a real Gauss
    bool gauss = (bincoef[0]>10000 && bincoef[1]>10000 && bincoef[2]>10000);
    
    createBinom(emit, alpha, beta, gamma, bincoef, gauss, beam, Np, isBinned);
  }
  else {
    if (Attributes::getString(itsAttr[DISTRIBUTION]) == "GAUSS") {
      
      double corr[3];

      corr[0] = Attributes::getReal(itsAttr[CORRX]);
      corr[1] = Attributes::getReal(itsAttr[CORRY]);
      corr[2] = Attributes::getReal(itsAttr[CORRT]);
      
      double Hs2a = Attributes::getReal(itsAttr[SIGMAX]);
      double Hs2b = Attributes::getReal(itsAttr[SIGMAPX]);
      
      double Vs2a = Attributes::getReal(itsAttr[SIGMAY]);
      double Vs2b = Attributes::getReal(itsAttr[SIGMAPY]);
      
      double Ls2a = Attributes::getReal(itsAttr[SIGMAT]);
      double Ls2b = Attributes::getReal(itsAttr[SIGMAPT]);
      
      RANLIB_class *rGen = new RANLIB_class(265314159,4);
 
      unsigned int pc = 0;
      unsigned int count = 0;      
  
      if (beam.isZPeriodic()) 
        INFOMSG("Distribution: uniform in z (periodic BC) z = " << -0.5*beam.getGaBeLa() << " ... " << 0.5*beam.getGaBeLa() << endl);
      
      for(int i=0; i<Np; i++) {        
	
        double x,y;       // generate independent Gaussians, then correlate and finaly scale 
	
        x  = rGen->gauss(0.0,1.0);
        y  = rGen->gauss(0.0,1.0);
	
        double px0  = x*corr[0] + y*sqrt(1.0 - corr[0]*corr[0]);
        double x0   =  x*Hs2a;
        px0 *= Hs2b;
	
        x  = rGen->gauss(0.0,1.0);
        y  = rGen->gauss(0.0,1.0);
        double py0  = x*corr[1] + y*sqrt(1.0 - corr[1]*corr[1]);
        double y0   =  x*Vs2a;
        py0 *= Vs2b;
	
        double del0;
        double psi0;
	
        if (beam.isZPeriodic()) {
          /*
            create uniform distribution
          */
          //del0 = (IpplRandom()-0.5)*p[2];
          //	psi0 = (IpplRandom()*beam.getGaBeLa()) - (0.5*beam.getGaBeLa());
	  
        } else {
          x  = rGen->gauss(0.0,1.0);
          y  = rGen->gauss(0.0,1.0);
          del0  = x*corr[2] + y*sqrt(1.0 - corr[2]*corr[2]);
          psi0   =  x*Ls2a;
          del0 *= Ls2b;
        }      
        if (pc == Ippl::myNode()) {
          beam.create(1);
          beam.R[count] = Vector_t(x0,y0,psi0);
          beam.P[count] = Vector_t(px0,py0,del0);
          beam.Bin[count] = 0; // not initialized
          count++;
        }
        pc++;
        if (pc == Ippl::getNodes())
          pc = 0;      
      }
    } else if (Attributes::getString(itsAttr[DISTRIBUTION]) == "UNIFORMXYZ"){

      double corr[3];

      corr[0] = Attributes::getReal(itsAttr[CORRX]);
      corr[1] = Attributes::getReal(itsAttr[CORRY]);
      corr[2] = Attributes::getReal(itsAttr[CORRT]);

      double Hs2a = Attributes::getReal(itsAttr[SIGMAX]);
      double Hs2b = Attributes::getReal(itsAttr[SIGMAPX]);
  
      double Vs2a = Attributes::getReal(itsAttr[SIGMAY]);
      double Vs2b = Attributes::getReal(itsAttr[SIGMAPY]);
  
      double Ls2a = Attributes::getReal(itsAttr[SIGMAT]);
      double Ls2b = Attributes::getReal(itsAttr[SIGMAPT]);

      double nBins = Attributes::getReal(itsAttr[NBIN]);
  
      RANLIB_class *rGen = new RANLIB_class(265314159,4);
 
      unsigned int pc = 0;
      unsigned int count = 0;      

      double gamma = 1 + (Ls2b*1.0E-6/0.511);
      double beta  = sqrt(1- (1/(gamma*gamma)));
      double bega  = beta*gamma;
  
      *gmsg << "* Gamma= " << gamma << " Beta= " << beta << endl;
  
      for(int i=0; i<Np; i++) {        

        double x,y,z,px,py,pz;       // generate independent Gaussians, then correlate and finaly scale 
        double R =3;

        while(R >1){
      
          x = rGen->uniform(-1.0,1.0);
          y = rGen->uniform(-1.0,1.0);
          z = rGen->uniform(-1.0,1.0);

          px = rGen->uniform(-1.0,1.0);
          py = rGen->uniform(-1.0,1.0);
          pz = rGen->uniform(-1.0,1.0);

          R = sqrt(x*x + y*y + z*z);

          // or can generate uniform distribution in 6D phase space by
          // using following line instead of above one. 
          // R = sqrt(x*x + y*y + z*z + px*px + py*py +pz*pz);
        }

        px  = x*corr[0] + px*sqrt(1.0 - corr[0]*corr[0]);
        x   = x*Hs2a;
        px *= Hs2b;
	    
        py  = y*corr[1] + py*sqrt(1.0 - corr[1]*corr[1]);
        y   = y*Vs2a;
        py *= Vs2b;

        pz  = z*corr[2] + pz*sqrt(1.0 - corr[2]*corr[2]);
        z   = z*Ls2a;
        pz *= Ls2b;

        if (pc == Ippl::myNode()) {
          beam.create(1);
          beam.R[count] = Vector_t(x,y,z);
          beam.P[count] = Vector_t(px,py,pz );
          beam.Bin[count] = 0; // not initialized
          count++;
        }
        pc++;
        if (pc == Ippl::getNodes())
          pc = 0;      
      }
      
    }else if (Attributes::getString(itsAttr[DISTRIBUTION]) == "FROMFILE") {

      *gmsg <<"---------------------------------------------"<<endl;
      *gmsg <<"     READ ININITAL DISTRIBUTION FROM FILE    "<<endl;
      const string filename=Attributes::getString(itsAttr[FNAME]);

      double x0,px0,y0,py0,psi0,del0;

      std::ifstream fs;
      fs.open(filename.c_str());

      if(fs.fail()){
        throw OpalException("Distribution::Create()",
                            "Open file operation failed,please check if \""
                            + filename +  "\" really exist.");
      }
    
      fs >> Np;
      if ( Np <= 0 ){
        throw OpalException("Distribution::Create()",
                            " The particle number should big than zero! please check,please check the first line of \""
                            + filename +  "\".");
      }
      unsigned int pc = 0;
      unsigned int count = 0;
      // read in particles number from the first line of the file.
    
      for(unsigned int i=0; i<Np; i++) {        

        if ( !fs.eof() ){
          if (pc == Ippl::myNode()) {
            // create 1 particle
            beam.create(1);

            fs>> x0 >> px0 >> y0 >> py0 >> psi0 >> del0; 
            // reasd in r and p
            beam.R[count] = Vector_t(x0,y0,psi0);
            beam.P[count] = Vector_t(px0,py0,del0);
            beam.Bin[count] = 0; // not initialized
            count++;
          }
        }
        else {
          throw OpalException("Distribution::Create()",
                              "File end reached before all particles imported, please check \""
                              + filename +  "\".");
          return;
        }

      }
      if (count != Np) {
      
        throw OpalException("Distribution::Create()",
                            "The actual readed particle number does NOT match Np (number of first line),please check \""
                            + filename +  "\".");
      }
    
      *gmsg <<"                     DONE !                  "<<endl;
      *gmsg <<"---------------------------------------------"<<endl;
      fs.close();

    }
  }	

  /*
    In the case of a binned distribution (gun)
    we have to do the boundp after emission.
  */
  if (!(isBinned)) {
    beam.boundp();
    beam.LastSection=0;
  }
}


void Distribution::binnDistribution(PartBunch &beam, size_t Np)
{
  double corr[3];

  corr[0] = Attributes::getReal(itsAttr[CORRX]);
  corr[1] = Attributes::getReal(itsAttr[CORRY]);
  corr[2] = Attributes::getReal(itsAttr[CORRT]);

  double Hs2a = Attributes::getReal(itsAttr[SIGMAX]);
  double Hs2b = Attributes::getReal(itsAttr[SIGMAPX]);
  
  double Vs2a = Attributes::getReal(itsAttr[SIGMAY]);
  double Vs2b = Attributes::getReal(itsAttr[SIGMAPY]);
  
  double Ls2a = Attributes::getReal(itsAttr[SIGMAT]);
  double Ls2b = Attributes::getReal(itsAttr[SIGMAPT]);

  double nBins = Attributes::getReal(itsAttr[NBIN]);

  double dEBins = Attributes::getReal(itsAttr[DEBIN]);
  
  RANLIB_class *rGen = new RANLIB_class(265314159,4);
 
  unsigned int pc = 0;
  unsigned int count = 0;      

  double gamma = 1 + (Ls2b*1.0E-6/0.511);
  double beta  = sqrt(1- (1/(gamma*gamma)));
  double bega  = beta*gamma;

  *gmsg << "* Gamma= " << gamma << " Beta= " << beta << endl;
  
  for(int i=0; i<Np; i++) {        

    double x,y;       // generate independent Gaussians, then correlate and finaly scale 
    double xy = 2;
    //x  = rGen->gauss(0.0,1.0);
    //y  = rGen->gauss(0.0,1.0);
    while (xy > 1){
      x  = rGen->uniform(-1.0,1.0);
      y  = rGen->uniform(-1.0,1.0);
      xy = sqrt(x*x + y*y);
    }
                
    //x  = rGen->gauss(0.0,1.0);
    //y  = rGen->gauss(0.0,1.0);
    double px0  = 0.0; //x*corr[0] + y*sqrt(1.0 - corr[0]*corr[0]);
    double x0   =  x*Hs2a;
    px0 *= Hs2b;
    
    //x  = rGen->gauss(0.0,1.0);
    //y  = rGen->gauss(0.0,1.0);
    double py0  = 0.0; //x*corr[1] + y*sqrt(1.0 - corr[1]*corr[1]);
    double y0   =  y*Vs2a;
    py0 *= Vs2b;

    double del0;
    double psi0;
    
    x  = rGen->gauss(0.0,1.0);
    y  = rGen->gauss(0.0,1.0);
   
    del0  = x*corr[2] + y*sqrt(1.0 - corr[2]*corr[2]);
    psi0  = x*Ls2a;
    del0 *= Ls2b;

    if (pc == Ippl::myNode()) {
      vector<double> tmp;
      tmp.push_back(x0);
      tmp.push_back(y0);
      tmp.push_back(psi0);
      //      tmp.push_back(px0);
      // tmp.push_back(py0);
      tmp.push_back(0.0);
      tmp.push_back(0.0);
      tmp.push_back(bega);
      tmp.push_back(0);
      pbin_m->fill(tmp);
      count++;
    }
    /*
      In case of binned distribution 
      we create all particles on node 0
      pc++;
      if (pc == Ippl::getNodes())
      pc = 0;      
      
    */
  }

  if (pc == Ippl::myNode())
    pbin_m->sortArray();  

  pbin_m->setRebinEnergy(dEBins);
  
  // now copy this over to the bunch
  // so that we can emmit the particles
  beam.setPBins(pbin_m);
}

void Distribution::doRestart(PartBunch &beam, size_t Np, size_t restartStep)
{

  H5PartFile *H5file;
  string fn = OPAL.getInputFn(); 
  int pos=fn.find(string("."),0);  
  fn.erase(pos,fn.size()-pos);

  beam.setTEmission(Attributes::getReal(itsAttr[TEMISSION]));
  
  fn += string(".h5");

#ifdef PARALLEL_IO
  H5file=H5PartOpenFileParallel(fn.c_str(),H5PART_READ,MPI_COMM_WORLD);
#else
  H5file=H5PartOpenFile(fn.c_str(),H5PART_READ);
#endif
  
  if(!H5file) {
    ERRORMSG("File open failed:  exiting!" << endl);
    exit(0);
  }
  
  H5PartSetStep(H5file,restartStep);
  int N=(int)H5PartGetNumParticles(H5file);

  h5part_int64_t totalSteps = H5PartGetNumSteps(H5file);
  *gmsg << "total number of particles = " << N << " total steps " << totalSteps << endl;

  //TODO: do a more sophisticated distribution of particles?
  //my guess is that the end range index is EXCLUSIVE!
  int numberOfParticlesPerNode = (int) floor((double) N / Ippl::getNodes());
  long long starti = Ippl::myNode() * numberOfParticlesPerNode;
  long long endi = 0;
  // ensure that we dont miss any particle in the end
  if(Ippl::myNode() == Ippl::getNodes() - 1)
    endi = -1;
  else
    endi = starti + numberOfParticlesPerNode;

  // cout << "Node " << Ippl::myNode() << " gets particles " << starti << " - " << endi << endl;

  H5PartSetView(H5file,starti,endi);
  N = (int)H5PartGetNumParticles(H5file);

  // cout << "Node " << Ippl::myNode() << " has " << N << " particles" << endl;
  
  double actualT;
  H5PartReadStepAttrib(H5file,"TIME",&actualT);

  beam.setT(actualT);

  /*
   * eventually read more attributes a la IMPACT-T
   */

  void *varray = malloc(N*sizeof(double)); 
  double *farray = (double*)varray;

  beam.create(N);
  
  H5PartReadDataFloat64(H5file,"x",farray);
  for (unsigned long int n=0; n < N; ++n) {
    beam.R[n](0)=farray[n];
    beam.Bin[n] = 0; // not initialized
  }
  H5PartReadDataFloat64(H5file,"y",farray);
  for (unsigned long int n=0; n < N; ++n)
    beam.R[n](1)=farray[n];

  H5PartReadDataFloat64(H5file,"z",farray);
  for (unsigned long int n=0; n < N; ++n)
    beam.R[n](2)=farray[n];
  
  H5PartReadDataFloat64(H5file,"px",farray);
  for (unsigned long int n=0; n < N; ++n)
    beam.P[n](0)=farray[n];

  H5PartReadDataFloat64(H5file,"py",farray);
  for (unsigned long int n=0; n < N; ++n)
    beam.P[n](1)=farray[n];

  H5PartReadDataFloat64(H5file,"pz",farray);
  for (unsigned long int n=0; n < N; ++n)
    beam.P[n](2)=farray[n];

  if(farray)
    free(farray);

  Ippl::Comm->barrier();
  H5PartCloseFile(H5file);
  beam.boundp();
  beam.LastSection=0;
}

void Distribution::doRestart_cycl(PartBunch &beam, size_t Np, size_t restartStep)
{

  H5PartFile *H5file;
  string fn = OPAL.getInputFn(); 
  int pos=fn.find(string("."),0);  
  fn.erase(pos,fn.size()-pos);

  fn += string(".h5");

#ifdef PARALLEL_IO
  H5file=H5PartOpenFileParallel(fn.c_str(),H5PART_READ,MPI_COMM_WORLD);
#else
  H5file=H5PartOpenFile(fn.c_str(),H5PART_READ);
#endif
  
  if(!H5file) {
    ERRORMSG("File open failed:  exiting!" << endl);
    exit(0);
  }
  
  H5PartSetStep(H5file,restartStep);
  int N=(int)H5PartGetNumParticles(H5file);

  h5part_int64_t totalSteps = H5PartGetNumSteps(H5file);
  *gmsg << "total number of particles = " << N << " total steps " << totalSteps << endl;

  //TODO: do a more sophisticated distribution of particles?
  //my guess is that the end range index is EXCLUSIVE!
  int numberOfParticlesPerNode = (int) floor((double) N / Ippl::getNodes());
  long long starti = Ippl::myNode() * numberOfParticlesPerNode;
  long long endi = 0;
  // ensure that we dont miss any particle in the end
  if(Ippl::myNode() == Ippl::getNodes() - 1)
    endi = -1;
  else
    endi = starti + numberOfParticlesPerNode;

  // cout << "Node " << Ippl::myNode() << " gets particles " << starti << " - " << endi << endl;

  H5PartSetView(H5file,starti,endi);
  N = (int)H5PartGetNumParticles(H5file);

  // cout << "Node " << Ippl::myNode() << " has " << N << " particles" << endl;
  
  double actualT;
  H5PartReadStepAttrib(H5file,"TIME",&actualT);

  beam.setT(actualT);

  double lpath;
  H5PartReadStepAttrib(H5file,"LPATH",&lpath);
  beam.setLPath(lpath);

  h5part_int64_t tstep;
  H5PartReadStepAttrib(H5file,"TrackStep",&tstep);
  beam.setTrackStep((long long)tstep);

  h5part_int64_t numBunch;
  H5PartReadStepAttrib(H5file,"NumBunch",&numBunch);
  beam.setNumBunch((int)numBunch);
    
  /*
   * eventually read more attributes a la IMPACT-T
   */

  void *varray = malloc(N*sizeof(double)); 
  double *farray = (double*)varray;

  h5part_int64_t *larray = (h5part_int64_t *)varray;
  
  beam.create(N);
  
  H5PartReadDataFloat64(H5file,"x",farray);
  for (unsigned long int n=0; n < N; ++n) {
    beam.R[n](0)=farray[n];
    beam.Bin[n] = 0; // not initialized
  }
  H5PartReadDataFloat64(H5file,"y",farray);
  for (unsigned long int n=0; n < N; ++n)
    beam.R[n](1)=farray[n];

  H5PartReadDataFloat64(H5file,"z",farray);
  for (unsigned long int n=0; n < N; ++n)
    beam.R[n](2)=farray[n];
  
  H5PartReadDataFloat64(H5file,"px",farray);
  for (unsigned long int n=0; n < N; ++n)
    beam.P[n](0)=farray[n];

  H5PartReadDataFloat64(H5file,"py",farray);
  for (unsigned long int n=0; n < N; ++n)
    beam.P[n](1)=farray[n];

  H5PartReadDataFloat64(H5file,"pz",farray);
  for (unsigned long int n=0; n < N; ++n)
    beam.P[n](2)=farray[n];

  H5PartReadDataInt64(H5file,"id",larray);   
  for (unsigned long int n=0; n < N; ++n)
    beam.ID[n]=larray[n];


  if(farray)
    free(farray);

  Ippl::Comm->barrier();
  H5PartCloseFile(H5file);
  beam.boundp();
  beam.LastSection=0;
}


Distribution *Distribution::find(const string &name)
{
  Distribution *beam = dynamic_cast<Distribution*>(OPAL.find(name));

  if (beam == 0) {
    throw OpalException("Distribution::find()", "Distribution \"" + name + "\" not found.");
  }

  return beam;
}

/*
  double Distribution::getET() const
  {
  return Attributes::getReal(itsAttr[ET]);
  }

  void Distribution::setET(double value)
  {
  Attributes::setReal(itsAttr[ET], value);
  }

*/

const PartData &Distribution::getReference() const
{
  // Cast away const, to allow logically constant Distribution to update.
  const_cast<Distribution *>(this)->update();
  return reference;
}

void Distribution::update()
{

}


void Distribution::tfsDescriptors(std::ostream &os) const
{
  os << "@ Distribution     %s  " << getOpalName() << '\n' ;
}
// Emittances added to tfsDescriptors by JMJ 4/4/2000


Inform &Distribution::print(Inform &os) const {

  os << "* ************* D I S T R I B U T I O N ******************************************** " << endl;
  if (!OPAL.inRestartRun()) {
    os << "* Distribution:\t" << getOpalName() << endl;
    os << "* Distribution type:\t" << Attributes::getString(itsAttr[DISTRIBUTION]) << endl;
    os << "* sigmax=\t" << Attributes::getReal(itsAttr[SIGMAX]);
    os << "\t sigmay=\t" << Attributes::getReal(itsAttr[SIGMAY]);
    os << "\t sigmat=\t" << Attributes::getReal(itsAttr[SIGMAT]) << endl;
    
    os << "* sigmapx=\t" << Attributes::getReal(itsAttr[SIGMAPX]);
    os << "\t sigmapy=\t" << Attributes::getReal(itsAttr[SIGMAPY]);
    os << "\t sigmapt=\t" << Attributes::getReal(itsAttr[SIGMAPT]) << endl;
    
    os << "* corr x-px=\t" << Attributes::getReal(itsAttr[CORRX]);
    os << "\t corr y=py=\t" << Attributes::getReal(itsAttr[CORRY]);
    os << "\t corr t-pt=\t" << Attributes::getReal(itsAttr[CORRT]) << endl;
    if(Attributes::getReal(itsAttr[TEMISSION]) > 0.0) {
      os << "* -------- G U N -------------------------------------------------------------------" << endl;
      os << "* temission [s] \t" << Attributes::getReal(itsAttr[TEMISSION]);
      os << "\t n bins \t" << (int) Attributes::getReal(itsAttr[NBIN]) << endl;
      os << "\t dEbin [keV] \t" << Attributes::getReal(itsAttr[DEBIN]) << endl;
      os << "* ----------------------------------------------------------------------------------" << endl;
    }
  }
  else
    os << "* Distribution from restart file" << endl;
  
  os << "* ********************************************************************************** " << endl;
}

/**
For testing the implementation of particle bins do the following:

Compile: g++ -DPartBinTest -o Bins Bins.cpp ranlib.o

Run test: ./Bins 10 10000 1.0

*/


PartBins::PartBins(int bins) :
  bins_m(bins),
  xmin_m(0.0),
  nemittedBins_m(0),
  xmax_m(0.0)
{ 

  nBin_m = new size_t[bins_m];
  xbinmin_m = new double[bins_m];
  xbinmax_m = new double[bins_m];
  binsEmitted_m = new bool[bins_m];

  for (int i=0; i<bins_m; i++) {	
    nBin_m[i] = 0;
    xbinmax_m[i] = -DBL_MAX; 
    xbinmin_m[i] = DBL_MAX;
    binsEmitted_m[i]=false;
  }
}

PartBins::~PartBins() 
{ 
  if (nBin_m) {
    delete nBin_m;
    delete xbinmax_m;
    delete xbinmin_m;
    delete binsEmitted_m;
  }
  tmppart_m.clear();
  isEmitted_m.clear();
}


/** /brief Update data which have been modified by node 0 */
void PartBins::updateBinStructure(){

  if (Ippl::myNode() != 0) {
    
    nemittedBins_m = 0;
    for (int i=0; i<bins_m; i++) {	
      binsEmitted_m[i]=false;
    }
  }
  reduce(nemittedBins_m,nemittedBins_m,OpAddAssign());
  reduce(binsEmitted_m, binsEmitted_m + bins_m, binsEmitted_m, OpBitwiseOrAssign());
}

  
#ifdef PartBinTest

void PartBins::emitBin() {
  /** push bunch out of cathode such that for the last particle R(2) = 0.0 */
  //   ofstream os; 
  //   os.open("data", ios::out);
  for (int b=0; b< tmppart_m.size(); b++) {
    for (int n=0; n< tmppart_m.size(); n++) {
      if (tmppart_m[n][6] == b) {
        vector<double> p;
        p = tmppart_m[n];
        p[2] -= xbinmin_m[b];   // shift to zero
	//         os << p[0] << "\t" << p[1] << "\t" << p[2] << "\t" << b << endl;
      }
    }
  }
  os.close();
#endif

bool PartBins::getPart(size_t n, int bin, vector<double> &p) {

  if (tmppart_m[n][6] == bin) {
    p = tmppart_m[n];
    return true;
  }
  else
    return false;
}

 void PartBins::sortArray() {
   
   /** sort the vector of particles such that position of the particles decreas with increasing index. 
       Then push the particles back by xmax_m + jifactor * bunch_length. In order that the method getBin(double x) works xmin_m has to be lowered a bit more.
   */
   double jifactor = 0.15;  
   double length;
   
   std::sort(tmppart_m.begin(), tmppart_m.end(), DescendingLocationSort(2));
   
   //   calcGlobalExtrema();
   xmax_m = tmppart_m[0][2];
   xmin_m = tmppart_m[tmppart_m.size() - 1][2];
   
   length = xmax_m - xmin_m;
   for (int n=0; n< tmppart_m.size(); n++) 
     tmppart_m[n][2] -= xmax_m + jifactor * length; /* push particles back */
   
   xmin_m -= xmax_m + 0.0001*(xmax_m - xmin_m) + jifactor * length; /* lower the limits */
   xmax_m -= xmax_m + jifactor * length;
   
   hBin_m = (fabs(xmax_m - xmin_m))/(bins_m);
   calcHBins();
   for (int n=0; n<bins_m; n++)
     if (nBin_m[n] == 0) setBinEmitted(n);
   
   calcExtrema();
 }
 
 void PartBins::calcHBins() {
   
   for (int n=0; n < tmppart_m.size(); n++)
     tmppart_m[n][6] = getBin(tmppart_m[n][2]);
  calcExtrema();
  
 }
  
size_t PartBins::getSum() {
  size_t s = 0;
  for (int n=0; n< bins_m; n++) 
    s += nBin_m[n];
  return s;
}
  
void PartBins::calcGlobalExtrema(){
  xmin_m = DBL_MAX;
  xmax_m = -DBL_MAX;
  for (int n=0; n< tmppart_m.size(); n++) {
    if(tmppart_m[n][2]<=xmin_m)
      xmin_m = tmppart_m[n][2];
    if(tmppart_m[n][2]>=xmax_m)
      xmax_m = tmppart_m[n][2];    
  }
  double xdiff = 0.01*(xmax_m - xmin_m);
  xmin_m -= xdiff;
  xmax_m += xdiff;
}

void PartBins::calcExtrema() {
  for (int n=0; n< tmppart_m.size(); n++) {
    if (xbinmin_m[(int)tmppart_m[n][6]] >= tmppart_m[n][2])
      xbinmin_m[(int)tmppart_m[n][6]] = tmppart_m[n][2];

    if (xbinmax_m[(int)tmppart_m[n][6]] <= tmppart_m[n][2])
      xbinmax_m[(int)tmppart_m[n][6]] = tmppart_m[n][2];
  }
}

Inform &PartBins::print(Inform &os) { 

  os <<"-----------------------------------------"<<endl;
  os <<"     CREATE BINNED GAUSS DISTRIBUTION DONE        "<<endl;

  os << "Bins= " << bins_m << " hBin= " << hBin_m << endl;
  os << "minVal= " << xmin_m << " [m] maxVal= " << xmax_m << " [m] Particle vector lenght " << tmppart_m.size() << endl;

  for (int i=0; i<bins_m; i++)	
    os << "Bin #" << i << " bin value " << nBin_m[i] << "\t minx= " << xbinmin_m[i]  << "\t maxx= " << xbinmax_m[i] << endl;

  if (getLastemittedBin() >= 0)
    os << "Last emitted bin is " << getLastemittedBin() << endl;
  else
    os << "No bin is emitted !" << endl;
  return os;
}

 int PartBins::getBin(double x) {
   /**
      returns the index of the bin to which the particle with z = 'x' belongs.
      If getBin returns b < 0 || b >= bins_m, then is x out of range!
   */
   int b = (int) floor(fabs(xmax_m - x)/hBin_m); 
   nBin_m[b]++;
   return b;
 }

#ifdef PartBinTest
int main(int argc, char *argv[])
{ 
  vector< vector<double> > xparr;

  /* initialize random seed: */
  // srand ( time(NULL) );

  RANLIB_class *rGen = new RANLIB_class(265314159,4);

  /*                         n-bins */
  
  PartBins *p = new PartBins(atoi(argv[1])); 

  for (int n=0; n< atoi(argv[2]); n++) {
    /* generate random numbers */

    double x = rGen->gauss(0.0,1.0); //  2.*atof(argv[3]) * (0.5 - rand() / ((double)RAND_MAX + 1)) ;
    double y = rGen->gauss(0.0,1.0); //  2.*atof(argv[3]) * (0.5 - rand() / ((double)RAND_MAX + 1)) ;
    double z = rGen->gauss(0.0,1.0); // 2.*atof(argv[3]) * (0.5 - rand() / ((double)RAND_MAX + 1)) ;
    vector<double> tmp;
    tmp.push_back(x);
    tmp.push_back(y);
    tmp.push_back(z);
    tmp.push_back(x);
    tmp.push_back(y);
    tmp.push_back(z);
    tmp.push_back(0.0); // bin
    p->fill(tmp);
  }

  p->sortArray();

  p->calcHBins();


  p->emitBin();
  
  cout << *p << endl;
  
  cout << "Value sum= " << p->getSum() << endl;

  if (p)
    delete p;
  return 0;
}
#endif


 
