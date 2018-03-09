// ------------------------------------------------------------------------
// $RCSfile: GTDistribution.cc,v $
// $Header: /afs/psi.ch/user/a/adelmann/private/cvsroot/GenTrackE/GTDistribution/GTDistribution.cc,v 1.7 2004/04/06 13:23:19 adelmann Exp $
// ------------------------------------------------------------------------
// $Revision: 1.7 $
// $State : $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
// Description:
//
// ------------------------------------------------------------------------
// Class category: 
// ------------------------------------------------------------------------
//
// $Date: 2004/04/06 13:23:19 $
// $Author: adelmann $
// $Log: GTDistribution.cc,v $
// Revision 1.7  2004/04/06 13:23:19  adelmann
// *** empty log message ***
//
// Revision 1.6  2004/04/05 13:07:51  adelmann
// Many changes to use new IPPL instead of POOMA
//
// Revision 1.5  2004/04/03 12:52:18  adelmann
// New gauss distribution based on cernlib
//
// Revision 1.4  2003/05/02 13:57:21  adelmann
// First electron cloud simulation 3d without SC and probable still with the wrong
// drive beam integration scheme. Otherwise :=) the program runs
// sureprisingly stable and the first inspection of the data looks promissing
//
// Revision 1.3  2003/04/17 14:19:52  adelmann
// *** empty log message ***
//
// Revision 1.2  2003/01/29 13:50:46  adelmann
// Fix READFROMFILE works now
//
// Revision 1.1.1.1  2003/01/23 09:13:58  adelmann
// GenTrackE
//
// ------------------------------------------------------------------------

#include <cmath>
#include "hammersley.h"

#ifdef USECMEE
extern "C" {
#include "cmee.h"
}
#endif

Distribution::~Distribution()
{
  
}

Distribution::Distribution(DistrData *distrData, double len) :
  distrData_m(distrData)
    
{

}

Distribution::Distribution()
{
}


Distribution::Distribution(DistrData *distrData, string fn) :
  distrData_m(distrData)
{
  Inform msg("Read in distribution: ");
  double arg1, arg2, arg3, arg4, arg5, arg6;
  
  ifstream ifstr;

  unsigned j = 0;
  if ( beam_m->isRoot()) {
    ifstr.open(fn.c_str(), ios::in );
    if (!(ifstr.bad())) { 
      msg << "Read particle positions for tracking from file: " << fn << endl ;
      while(ifstr >> arg1 >> arg2 >> arg3 >> arg4 >> arg5 >> arg6) {
	beam_m->create(1);
	beam_m->R[j] = Vector_t(arg1,arg3,arg5);
	beam_m->P[j] = Vector_t(arg2,arg4,arg6);
	j++;
      }
      msg << "Done reading " << j << " particles " << endl;
    }
    ifstr.close();
  }
}


#ifdef TTRACK
/** 
 * 
 * 
 * @param beam the particle container to fill
 * @param particles number of total particles
 * @param xmean 
 * @param pmean 
 * @param xstddev 
 * @param pstddev 
 * @param angle 
 */
void Distribution::generateGaussianTT(ChargedParticles<playout_t> *beam,
				    unsigned long int particles,
				    Vector_t xmean,
				    Vector_t pmean,
				    Vector_t xstddev,
				    Vector_t pstddev,
				    Vector_t angle) 
{
  RANLIB_class *rGen = new RANLIB_class(265314159,4);
  unsigned int p = 0;
  unsigned int count = 0; 

  for(int i=0; i<particles; i++) {        
    // generate independent Gaussians:
    double xi   = rGen->gauss(xmean(0),xstddev(0));
    double pxi  = rGen->gauss(pmean(0),pstddev(0));
    double yi   = rGen->gauss(xmean(1),xstddev(1));
    double pyi  = rGen->gauss(pmean(1),pstddev(1));
    double psii = rGen->gauss(xmean(2),xstddev(2));
    double deli = rGen->gauss(pmean(2),pstddev(2));
    
    double X0   = xi*cos(angle(0))-pxi*sin(angle(0));
    double Px0  = xi*sin(angle(0))+pxi*cos(angle(0));
    double Y0   = yi * cos(angle(1)) - pyi * sin(angle(1));
    double Py0  = yi * sin(angle(1)) + pyi * cos(angle(1));
    double Psi0 = psii * cos(angle(2)) - deli * sin(angle(2));
    double Del0 = psii * sin(angle(2)) + deli * cos(angle(2));
      
    if (p == Ippl::myNode()) {
      beam->create(1);
      beam->R[count] = Vector_t(X0,Y0,Psi0) + xmean;
      beam->P[count] = Vector_t(Px0,Py0,Del0) + pmean;
      beam->Q[count] = beam->getQ();
      count++;
    }
    p++;
    if (p == Ippl::getNodes())
      p = 0;
  }
  beam->boundp();
}


/** 
 * 
 * 
 * @param beam the particle container to fill
 * @param particles number of total particles
 * @param xmean 
 * @param pmean 
 * @param xstddev 
 * @param pstddev 
 * @param angle 
 */
void Distribution::generateGaussianFlatTop(ChargedParticles<playout_t> *beam,
				    unsigned long int particles,
				    Vector_t xmean,
				    Vector_t pmean,
				    Vector_t xstddev,
				    Vector_t pstddev,
				    Vector_t angle) 
{
  RANLIB_class *rGen = new RANLIB_class(265314159,4);
  unsigned int p = 0;
  unsigned int count = 0; 

  for(int i=0; i<particles; i++) {        
    // generate independent Gaussians:
    double xi   = rGen->gauss(xmean(0),xstddev(0));
    double pxi  = rGen->gauss(pmean(0),pstddev(0));
    double yi   = rGen->gauss(xmean(1),xstddev(1));
    double pyi  = rGen->gauss(pmean(1),pstddev(1));
    double psii = (2.0*IpplRandom() - 1.0) * xstddev(2);
    double deli = rGen->gauss(pmean(2),pstddev(2));
    
    double X0   = xi*cos(angle(0))-pxi*sin(angle(0));
    double Px0  = xi*sin(angle(0))+pxi*cos(angle(0));
    double Y0   = yi * cos(angle(1)) - pyi * sin(angle(1));
    double Py0  = yi * sin(angle(1)) + pyi * cos(angle(1));
    double Psi0 = psii * cos(angle(2)) - deli * sin(angle(2));
    double Del0 = psii * sin(angle(2)) + deli * cos(angle(2));
      
    if (p == Ippl::myNode()) {
      beam->create(1);
      beam->R[count] = Vector_t(X0,Y0,Psi0) + xmean;
      beam->P[count] = Vector_t(Px0,Py0,Del0) + pmean;
      beam->Q[count] = beam->getQ();
      count++;
    }
    p++;
    if (p == Ippl::getNodes())
      p = 0;
  }
  beam->boundp();
}

/** 
 * 
 * 
 * @param beam the particle container to fill
 * @param particles number of total particles
 * @param xmean 
 * @param pmean 
 * @param xstddev 
 * @param pstddev 
 * @param angle 
 */
void Distribution::generateGaussianFlatTopHam(ChargedParticles<playout_t> *beam,
				    unsigned long int particles,
				    Vector_t xmean,
				    Vector_t pmean,
				    Vector_t xstddev,
				    Vector_t pstddev,
				    Vector_t angle) 
{
  // RANLIB_class *rGen = new RANLIB_class(265314159,4);
  unsigned int p = 0;
  unsigned int count = 0; 

  unsigned int DIM = 13;

  double r[DIM];
  int seed[DIM], leap[DIM], base[DIM];
  for(unsigned i=0; i < DIM;i++) seed[i] = 0;
  for(unsigned i=0; i < DIM;i++) leap[i] = 1;

  base[0] = - particles;
  base[1] = 2; base[2] = 3; base[3] = 5; base[4] = 7; base[5] = 11;
  base[6] = 13; base[7] = 17; base[8] = 19; base[9] = 23; base[10] = 29;
  base[11] = 31; base[12] = 37;

  for(int i=0; i<((double)particles) / Ippl::getNodes(); i++) {        
    i_to_hammersley_sequence(DIM,1,1 + (int)(((double)particles) * Ippl::myNode() / Ippl::getNodes() + i * DIM),    
    			     seed, leap, base, r);

    // generate independent Gaussians:
    double xi   = gauss(xmean(0),xstddev(0),r[12],r[1]);
    double pxi  = gauss(pmean(0),pstddev(0),r[2],r[3]);
    double yi   = gauss(xmean(1),xstddev(1),r[4],r[5]);
    double pyi  = gauss(pmean(1),pstddev(1),r[6],r[7]);
    double psii = (2.0*r[8] - 1.0) * xstddev(2);
    double deli = gauss(pmean(2),pstddev(2),r[10],r[11]);
    
    double X0   = xi*cos(angle(0))-pxi*sin(angle(0));
    double Px0  = xi*sin(angle(0))+pxi*cos(angle(0));
    double Y0   = yi * cos(angle(1)) - pyi * sin(angle(1));
    double Py0  = yi * sin(angle(1)) + pyi * cos(angle(1));
    double Psi0 = psii * cos(angle(2)) - deli * sin(angle(2));
    double Del0 = psii * sin(angle(2)) + deli * cos(angle(2));
      
    
    beam->create(1);
    beam->R[count] = Vector_t(X0,Y0,Psi0) + xmean;
    beam->P[count] = Vector_t(Px0,Py0,Del0) + pmean;
    beam->Q[count] = beam->getQ();
    count++;
  }
  beam->boundp();
}
#else
void Distribution::generateGaussian(ChargedParticles<playout_t> *beam,
				    unsigned long int particles,
				    Vector_t xmean,
				    Vector_t pmean,
				    Vector_t xstddev,
				    Vector_t pstddev,
				    Vector_t angle) 
{
  RANLIB_class *rGen = new RANLIB_class(265314159,4);
  unsigned int p = 0;
  unsigned int count = 0; 

  for(int i=0; i<particles; i++) {        
    // generate independent Gaussians:
    double xi   = rGen->gauss(xmean(0),xstddev(0));
    double pxi  = rGen->gauss(pmean(0),pstddev(0));
    double yi   = rGen->gauss(xmean(1),xstddev(1));
    double pyi  = rGen->gauss(pmean(1),pstddev(1));
    double psii = rGen->gauss(xmean(2),xstddev(2));
    double deli = rGen->gauss(pmean(2),pstddev(2));
    
    double X0   = xi*cos(angle(0))-pxi*sin(angle(0));
    double Px0  = xi*sin(angle(0))+pxi*cos(angle(0));
    double Y0   = yi * cos(angle(1)) - pyi * sin(angle(1));
    double Py0  = yi * sin(angle(1)) + pyi * cos(angle(1));
    double Psi0 = psii * cos(angle(2)) - deli * sin(angle(2));
    double Del0 = psii * sin(angle(2)) + deli * cos(angle(2));
      
    Vector_t x = Vector_t(X0,Y0,Psi0) + xmean;
      
    if (p == Ippl::myNode()) {
      beam->create(1);
      beam->R[count] = x; 
      beam->P[count] = Vector_t(0.0,0.0,Physics::c*beam->getBeta()); 
      beam->ndt[count] = 0;       // the individual time step of the particle [s]
      beam->dtleft[count] = 0;    // the time left in the integration scheme  [s]
      beam->elemNo[count] = 0;    // marks the location of the particle in the lattice
      beam->stat[count] = ALIVE;  // status of particle
      count++;
    }
    p++;
    if (p == Ippl::getNodes())
      p = 0;
  }
  beam->boundp();
}



void Distribution::createUniformElectronBackground(ChargedParticles<playout_t> *beam,
						   unsigned long totalP, double rCh, double lCh, double maxE) {
  
  /*
    totalP  :  total number of electrons
    rCh     :  radius of the circular vaccum chamber
    lCh     :  lenght of the circular vaccum chamber
    maxE    :  maximum kinetic energy of the electrons [eV]
  */

  const double cl = Physics::c;



  double r,rx,ry;
  double px,py,pz;
  unsigned int p = 0;
  unsigned k = 0;
  
  int ns[1];
  double bn[10];
  double bt[10];
  double bz[10];
  
  /* The material the electron hits is described
   * by an integer (1=copper, 2=stainless steel)
   */
  int mat_number = 1;

  for (int n=0; n<totalP; n++) {
    do {
      rx = 2.0*IpplRandom() - 1.0;    // -1.0 .... 1.0
      ry = 2.0*IpplRandom() - 1.0;    // -1.0 .... 1.0
      py = IpplRandom();              //  0.0 .... 1.0
      r = rx*rx + ry*ry;
    } while (r>1.0);
    

    if (p == Ippl::myNode()) {
      beam->create(1);
      beam->R[k](0) = rx*rCh;
      beam->R[k](1) = ry*rCh;
      beam->R[k](2) = lCh*IpplRandom();  //  0.0 .... l

      px = 1.0 + maxE*IpplRandom()*1E-9/Physics::m_e;   // component in eV
      py = 1.0 + maxE*IpplRandom()*1E-9/Physics::m_e;  // component in eV
      pz = 1.0 + (maxE-(0.5*(px+py)))*1E-9/Physics::m_e;

      if ( IpplRandom()>0.5)
	beam->P[k] = Vector_t(-cl*sqrt(1.0-pow(px,-2.0)),
			      -cl*sqrt(1.0-pow(py,-2.0)),
			      cl*sqrt(1.0-pow(pz,-2.0)));
      else
	beam->P[k] = Vector_t(cl*sqrt(1.0-pow(px,-2.0)),
			      cl*sqrt(1.0-pow(py,-2.0)),
			      cl*sqrt(1.0-pow(pz,-2.0)));
      
      beam->ndt[k] = 0;       // the individual time step of the particle [s]
      beam->dtleft[k] = 0;    // the time left in the integration scheme  [s]
      beam->elemNo[k] = 0;    // marks the location of the particle in the lattice
      beam->stat[k] = ALIVE;  // status of particle
      k++;
    }
    p++;
    if (p == Ippl::getNodes())
      p=0;
  }
  beam->boundp();
}


void Distribution::createBinom(unsigned long int particles,
				      Vector_t emit,
				      Vector_t alpha,
				      Vector_t beta,
				      Vector_t bincoef)
{ 
  using Physics::two_pi; 
   
  Vector_t gamma;

  gamma[0]= (1.0 + sqr(alpha[0]))/beta[0];
  gamma[1]= (1.0 + sqr(alpha[1]))/beta[1];
  gamma[2]= (1.0 + sqr(alpha[2]))/beta[2];
 
  beta *= 4.0;
  gamma *= 4.0;

  double XM=sqrt(emit[0]*beta[0]);
  double XPM=sqrt(emit[0]*gamma[0]);
  double COSCHIX=sqrt(1.0/(1.0+alpha[0]*alpha[0]));

  double SINCHIX=-alpha[0]*COSCHIX;
  double CHIX=atan2(SINCHIX,COSCHIX);
  double AMIX=1.0/bincoef[0];
  double XL=sqrt((bincoef[0]+1.0)/2.0)*XM;
  double XPL=sqrt((bincoef[0]+1.0)/2.0)*XPM;

   // Y and YP
  double YM=sqrt(emit[1]*beta[1]);
  double YPM=sqrt(emit[1]*gamma[1]);
  double COSCHIY=sqrt(1.0/(1.0+alpha[1]*alpha[1]));
  double SINCHIY=-alpha[1]*COSCHIY;
  double CHIY=atan2(SINCHIY,COSCHIY);
  double AMIY=1.0/bincoef[1];
  double YL=sqrt((bincoef[1]+1.0)/2.0)*YM;
  double YPL=sqrt((bincoef[1]+1.0)/2.0)*YPM;

  // T and TP
  double TM=sqrt(emit[2]*beta[2]);
  double TPM=sqrt(emit[2]*gamma[2]);
  double COSCHIT=sqrt(1.0/(1.0+alpha[2]*alpha[2]));
  double SINCHIT=-alpha[2]*COSCHIT;
  double CHIT=atan2(SINCHIT,COSCHIT);
  double AMIT=1.0/bincoef[2];
  double PHIT=sqrt((bincoef[2]+1.0)/2.0)*TM;
  double PHIPT=sqrt((bincoef[2]+1.0)/2.0)*TPM;
  

  unsigned int p = 0; 
  unsigned int k = 0;
  
  for (int n=0; n<2*(3-1); n++)
    IpplRandom();
  
  // what follows is arbritrary ??? and defines the distribution is
  // gaussian, if bincoef[0] is bigger than 1000. Have to check if this is 
  // reasonable !
  if (bincoef[0] < 1000) {
    for (unsigned long int n=0; n < particles; ++n) {
      double S1=IpplRandom();
      double S2=IpplRandom();
      double S3=IpplRandom();
      double S4=IpplRandom();
      double S5=IpplRandom();
      double S6=IpplRandom();
      if (p == Ippl::myNode()) {
	beam_m->create(1);

	double A=sqrt(1.0 - pow(S1,AMIX));
	double AL=two_pi*S2;
	double U=A*cos(AL);
	double V=A*sin(AL);
	beam_m->R[k](0)=XL*U;                     // X[space or momentum][part](x,y,t)
	beam_m->P[k](0)=XPL*(U*SINCHIX+V*COSCHIX);

	A=sqrt(1.0 - pow(S3,AMIY));
	AL=two_pi*S4;
	U=A*cos(AL);
	V=A*sin(AL);
	beam_m->R[k](1)=YL*U;
	beam_m->P[k](1)=YPL*(U*SINCHIY+V*COSCHIY);
	
        A=sqrt(1.0 - pow(S5,AMIY));
	AL=two_pi*S6;
	U=A*cos(AL);
	V=A*sin(AL);
	beam_m->R[k](2)=PHIT*U;
	beam_m->P[k](2)=PHIPT*(U*SINCHIT+V*COSCHIT); 

	k++;
      }
      p++;
      if (p == Ippl::getNodes())
	p=0;
    }
  }
  else {
    for (unsigned long int n=0; n < particles; ++n) {
      double S1=IpplRandom();
      double S2=IpplRandom();
      double S3=IpplRandom();
      double S4=IpplRandom();
      double S5=IpplRandom();
      double S6=IpplRandom();
       if (p == Ippl::myNode()) {
	 beam_m->create(1);
	 if (S1 == 0.0)
	   S1 += 0.00000000001;
	 double A=sqrt(2.0)/2.0 * sqrt(-log(S1));
	 double AL=two_pi*S2;
	 double U=A*cos(AL);
	 double V=A*sin(AL);
	 beam_m->R[k](0)=XM*U;
	 beam_m->P[k](0)=XPM*(U*SINCHIX+V*COSCHIX);
      
	 if (S3 == 0.0)
	   S3 += 0.00000000001;
	 A=sqrt(2.0)/2.0 * sqrt(-log(S3));
	 AL=two_pi*S4;
	 U=A*cos(AL);
	 V=A*sin(AL);
	 beam_m->R[k](1)=YM*U;
	 beam_m->P[k](1)=YPM*(U*SINCHIY+V*COSCHIY);
	
	 if (S5 == 0.0)
	   S5 += 0.00000000001;
	 A=sqrt(2.0)/2.0 * sqrt(-log(S5));
	 AL=two_pi*S6;
	 U=A*cos(AL);
	 V=A*sin(AL);
	 beam_m->R[k](2)=TM*U;
	 beam_m->P[k](2)=TPM*(U*SINCHIT+V*COSCHIT);

	 k++;
       }
       p++;
       if (p == Ippl::getNodes())
	p=0;
    }
  }
}
void Distribution::regenerateBinom(unsigned long int particles,
				      Vector_t emit,
				      Vector_t alpha,
				      Vector_t beta,
				      Vector_t bincoef)
{ 
  using Physics::two_pi; 
   
  Vector_t gamma;

  gamma[0]= (1.0 + sqr(alpha[0]))/beta[0];
  gamma[1]= (1.0 + sqr(alpha[1]))/beta[1];
  gamma[2]= (1.0 + sqr(alpha[2]))/beta[2];
 
  beta *= 4.0;
  gamma *= 4.0;

  double XM=sqrt(emit[0]*beta[0]);
  double XPM=sqrt(emit[0]*gamma[0]);
  double COSCHIX=sqrt(1.0/(1.0+alpha[0]*alpha[0]));
  double SINCHIX=-alpha[0]*COSCHIX;
  double CHIX=atan2(SINCHIX,COSCHIX);
  double AMIX=1.0/bincoef[0];
  double XL=sqrt((bincoef[0]+1.0)/2.0)*XM;
  double XPL=sqrt((bincoef[0]+1.0)/2.0)*XPM;

   // Y and YP
  double YM=sqrt(emit[1]*beta[1]);
  double YPM=sqrt(emit[1]*gamma[1]);
  double COSCHIY=sqrt(1.0/(1.0+alpha[1]*alpha[1]));
  double SINCHIY=-alpha[1]*COSCHIY;
  double CHIY=atan2(SINCHIY,COSCHIY);
  double AMIY=1.0/bincoef[1];
  double YL=sqrt((bincoef[1]+1.0)/2.0)*YM;
  double YPL=sqrt((bincoef[1]+1.0)/2.0)*YPM;

  // T and TP
  double TM=sqrt(emit[2]*beta[2]);
  double TPM=sqrt(emit[2]*gamma[2]);
  double COSCHIT=sqrt(1.0/(1.0+alpha[2]*alpha[2]));
  double SINCHIT=-alpha[2]*COSCHIT;
  double CHIT=atan2(SINCHIT,COSCHIT);
  double AMIT=1.0/bincoef[2];
  double PHIT=sqrt((bincoef[2]+1.0)/2.0)*TM;
  double PHIPT=sqrt((bincoef[2]+1.0)/2.0)*TPM;
  

  unsigned int p = 0; 
  unsigned int k = 0;
  
  for (int n=0; n<2*(3-1); n++)
    IpplRandom();
  
  // what follows is arbritrary ??? and defines the distribution is
  // gaussian, if bincoef[0] is bigger than 1000. Have to check if this is 
  // reasonable !
  if (bincoef[0] < 1000) {
    for (unsigned long int n=0; n < particles; ++n) {
      double S1=IpplRandom();
      double S2=IpplRandom();
      double S3=IpplRandom();
      double S4=IpplRandom();
      double S5=IpplRandom();
      double S6=IpplRandom();
      if (p == Ippl::myNode()) {
	double A=sqrt(1.0 - pow(S1,AMIX));
	double AL=two_pi*S2;
	double U=A*cos(AL);
	double V=A*sin(AL);
	beam_m->R[k](0)=XL*U;                     // X[space or momentum][part](x,y,t)
	beam_m->P[k](0)=XPL*(U*SINCHIX+V*COSCHIX);

	A=sqrt(1.0 - pow(S3,AMIY));
	AL=two_pi*S4;
	U=A*cos(AL);
	V=A*sin(AL);
	beam_m->R[k](1)=YL*U;
	beam_m->P[k](1)=YPL*(U*SINCHIY+V*COSCHIY);
	
        A=sqrt(1.0 - pow(S5,AMIY));
	AL=two_pi*S6;
	U=A*cos(AL);
	V=A*sin(AL);
	beam_m->R[k](2)=PHIT*U;
	beam_m->P[k](2)=PHIPT*(U*SINCHIT+V*COSCHIT); 

	k++;
      }
      p++;
      if (p == Ippl::getNodes())
	p=0;
    }
  }
  else {
    for (unsigned long int n=0; n < particles; ++n) {
      double S1=IpplRandom();
      double S2=IpplRandom();
      double S3=IpplRandom();
      double S4=IpplRandom();
      double S5=IpplRandom();
      double S6=IpplRandom();
       if (p == Ippl::myNode()) {
	 if (S1 == 0.0)
	   S1 += 0.00000000001;
	 double A=sqrt(2.0)/2.0 * sqrt(-log(S1));
	 double AL=two_pi*S2;
	 double U=A*cos(AL);
	 double V=A*sin(AL);
	 beam_m->R[k](0)=XM*U;
	 beam_m->P[k](0)=XPM*(U*SINCHIX+V*COSCHIX);
      
	 if (S3 == 0.0)
	   S3 += 0.00000000001;
	 A=sqrt(2.0)/2.0 * sqrt(-log(S3));
	 AL=two_pi*S4;
	 U=A*cos(AL);
	 V=A*sin(AL);
	 beam_m->R[k](1)=YM*U;
	 beam_m->P[k](1)=YPM*(U*SINCHIY+V*COSCHIY);
	
	 if (S5 == 0.0)
	   S5 += 0.00000000001;
	 A=sqrt(2.0)/2.0 * sqrt(-log(S5));
	 AL=two_pi*S6;
	 U=A*cos(AL);
	 V=A*sin(AL);
	 beam_m->R[k](2)=TM*U;
	 beam_m->P[k](2)=TPM*(U*SINCHIT+V*COSCHIT);

	 k++;
       }
       p++;
       if (p == Ippl::getNodes())
	p=0;
    }
  }
}

void Distribution::createUniformDistrEllipsoid(double dtime, unsigned long totalP, 
						      Vector_t semiaxi, Vector_t center, Vector_t vel ) {
    
    double r,rx,ry,rz;
    unsigned int p = 0;
    unsigned k = 0; // ada beam_m->getTopDriveBeam();
 
    for (int n=0; n<totalP; n++) {
      do {
	rx = 2.0*IpplRandom() - 1.0;
	ry = 2.0*IpplRandom() - 1.0;
	rz = 2.0*IpplRandom() - 1.0;
	r = 0.0;
	if (p == Ippl::myNode())
	  r = rx*rx + ry*ry + rz*rz;
      } while (r>1.0);
      if (p == Ippl::myNode()) {
	beam_m->create(1);
	beam_m->R[k](0) = rx*semiaxi(0);
	beam_m->R[k](1) = ry*semiaxi(1);
	beam_m->R[k](2) = rz*semiaxi(2);
	beam_m->P[k](0) = 0.0; //vel(0) + semiaxi(0)/dtime*1e-4*(2*IpplRandom()-1);
	beam_m->P[k](1) = 0.0; // vel(1) + semiaxi(0)/dtime*1e-4*(2*IpplRandom()-1);
	beam_m->P[k](2) = 0.0; //vel(2) + semiaxi(0)/dtime*1e-4*(2*IpplRandom()-1);
	k++;
      }
      p++;
      if (p == Ippl::getNodes())
	p=0;
    }
}

void Distribution::createTwoStreamInstabDistribution(unsigned long Np, double protonMacroPartRatio,
						     unsigned long Ne, double electronMacroPartRatio,
						     double rp,  // transverse radius
						     double lp,  // longitudinal lengt
						     double dlp, // longitudinal displacement
						     double re,  // transverse radius
						     double le,  // longitudinal lengt
						     double dle) // longitudinal displacement 
{

  // ada NEW    double localMin[3],double localMax[3],double grid_stretch) {
  
  Inform msg2all("createTwoStreamInstabDistribution ",INFORM_ALL_NODES);  

  double dt = 1.0;

  unsigned long n=0;    
    
  double vz = Physics::c*beam_m->getBeta();
  Vektor<double,3> p0(0.0,0.0,vz);
  
  unsigned long np = Np / Ippl::getNodes();

  msg2all << "Ne = " << Ne << "\tNp= " << Np << "\t re= " << re << "\t rp= " << rp << endl;
    
  while (np <= Np) {
    np = (Np - beam_m->getTotalNum()) / Ippl::getNodes() + 1;
    for (unsigned i=0;i<np;i++) {  
      Vektor<double,3> x = Vektor<double,3>(2.0*rp*(0.5 - IpplRandom()),
					    2.0*rp*(0.5 - IpplRandom()),
					    0.0);
      if (dot(x,x) <= rp*rp) {
	x = x + Vektor<double,3>(0.0,0.0,dlp+2.0*lp*(0.5 - IpplRandom()));
	beam_m->create(1);
	beam_m->R[n] = x;
	beam_m->P[n] = p0;
	beam_m->E[n] = Vektor<double,3>(0.0,0.0,0.0);
	beam_m->B[n] = Vektor<double,3>(0.0,0.0,0.0);
	beam_m->ndt[n] = 0;       // the individual time step of the particle [s]
	beam_m->dtleft[n] = 0;    // the time left in the integration scheme  [s]
	beam_m->elemNo[n] = 0;    // marks the location of the particle in the lattice
	beam_m->stat[n] = ALIVE;  // status of particle
	n++; 
      }
    }
    beam_m->update();
    np = beam_m->getTotalNum();
  }
  msg2all << "Nptot= " << np << " Nploc= " << beam_m->getLocalNum() << " Qp= " << np*beam_m->getQ() << endl;
  
  /* --------------------------------------------------------------------------------------------------
       Electron stuff
       --------------------------------------------------------------------------------------------------*/
    unsigned long ne = 0;
    
    if (Ne>0) {
      while (ne <= Ne) {
	ne = (Ne+np-beam_m->getTotalNum()) / Ippl::getNodes() + 1;
	for (unsigned i=0; i<ne;i++) {  
	  Vektor<double,3> x = Vektor<double,3>(IpplRandom(),
						IpplRandom(),
						IpplRandom());
	  beam_m->create(1);
	  beam_m->R[n] = x;
	  beam_m->P[n] = Vektor<double,3>(0.0,0.0,0.0);
	  beam_m->E[n] = Vektor<double,3>(0.0,0.0,0.0);
	  beam_m->B[n] = Vektor<double,3>(0.0,0.0,0.0);
	  beam_m->ndt[n] = 0;       // the individual time step of the particle [s]
	  beam_m->dtleft[n] = 0;    // the time left in the integration scheme  [s]
	  beam_m->elemNo[n] = 0;    // marks the location of the particle in the lattice
	  beam_m->stat[n] = ALIVE;  // status of particle
	  n++; 
	}
	beam_m->update();
	ne = beam_m->getTotalNum() - np;
      }
      msg2all << "Netot= " << ne << " Qe= " << -Physics::q_e*ne << endl;
    }
}

#endif




void Distribution::calcError()
{
  double pcnt=  static_cast<double>(beam_m->getTotalNum());

  Vector_t rm   = sum(beam_m->R);
  Vector_t pm   = sum(beam_m->P);
  Vector_t r2   = sum(beam_m->R*beam_m->R);
  Vector_t p2   = sum(beam_m->P*beam_m->P);
  Vector_t rp   = sum(beam_m->R*beam_m->P);
  
  Vector_t r4   = sum(beam_m->R*beam_m->R*beam_m->R*beam_m->R);
  Vector_t p4   = sum(beam_m->P*beam_m->P*beam_m->P*beam_m->P);
  Vector_t r3p  = sum(beam_m->R*beam_m->R*beam_m->R)*sum(beam_m->P);
  Vector_t rp3  = sum(beam_m->P*beam_m->P*beam_m->P)*sum(beam_m->R);
  Vector_t r2p2 = sum(beam_m->P*beam_m->P)*sum(beam_m->R*beam_m->R);

  // calculate emittances
  Vector_t emit2 = (r2*p2 - rp*rp)/pcnt;
  
  double EmitH = sqrt(emit2[0]); 
  double EmitV = sqrt(emit2[1]); 
  double EmitT = sqrt(emit2[2]);
  
  
  /*
   * Calculate all Covariances
   * ==========================
   */    
  double CxxH = r2(0)/pcnt; double CxxV = r2(1)/pcnt;
  double CxyH = rp(0)/pcnt; double CxyV = rp(1)/pcnt;
  double CyyH = p2(0)/pcnt; double CyyV = p2(1)/pcnt;
  
  double CxxT = r2(2)/pcnt; 
  double CxyT = rp(2)/pcnt; 
  double CyyT = p2(2)/pcnt; 

  /*
   * Calculate THE error matrix
   * ==========================
   */ 

  // horizontal:
  double CxxCxxH =(r4[0]/pcnt    - CxxH*CxxH) / pcnt;
  double CxyCxyH =(r2p2[0]/pcnt - CxyH*CxyH) / pcnt;
  double CyyCyyH =(p4[0]/pcnt   - CyyH*CyyH) / pcnt;
    
  double CxxCxyH =(r3p[0]/pcnt  - CxxH*CxyH) / pcnt;
  double CxxCyyH =(r2p2[0]/pcnt - CxxH*CyyH) / pcnt;
  double CxyCyyH =(rp3[0]/pcnt  - CxyH*CyyH) / pcnt;
    
  double EmatH[3][3] = {{CxxCxxH, CxxCxyH, CxxCyyH},
			{CxxCxyH, CxyCxyH, CxyCyyH},
			{CxxCyyH, CxyCyyH, CyyCyyH}};
  
  // Vertical:
  double CxxCxxV =(r4[1]/pcnt    - CxxV*CxxV) / pcnt;
  double CxyCxyV =(r2p2[1]/pcnt - CxyV*CxyV) / pcnt;
  double CyyCyyV =(p4[1]/pcnt   - CyyV*CyyV) / pcnt;
  
  double CxxCxyV =(r3p[1]/pcnt  - CxxV*CxyV) / pcnt;
  double CxxCyyV =(r2p2[1]/pcnt - CxxV*CyyV) / pcnt;
  double CxyCyyV =(rp3[1]/pcnt  - CxyV*CyyV) / pcnt;
  
  double EmatV[3][3] = {{CxxCxxV, CxxCxyV, CxxCyyV},
			{CxxCxyV, CxyCxyV, CxyCyyV},
			{CxxCyyV, CxyCyyV, CyyCyyV}};
  
  // Longitudinal:
  double CxxCxxT =(r4[2]/pcnt    - CxxT*CxxT) / pcnt;
  double CxyCxyT =(r2p2[2]/pcnt - CxyT*CxyT) / pcnt;
  double CyyCyyT =(p4[2]/pcnt   - CyyT*CyyT) / pcnt;
  
  double CxxCxyT =(r3p[2]/pcnt  - CxxT*CxyT) / pcnt;
  double CxxCyyT =(r2p2[2]/pcnt - CxxT*CyyT) / pcnt;
  double CxyCyyT =(rp3[2]/pcnt  - CxyT*CyyT) / pcnt;
  
  double EmatT[3][3] = {{CxxCxxT, CxxCxyT, CxxCyyT},
			{CxxCxyT, CxyCxyT, CxyCyyT},
			{CxxCyyT, CxyCyyT, CyyCyyT}};
  
 /*
  * Calculate Twiss parameters (TRANSVERSE ONLY)
  * ============================================
  */
    
  double AlphaH = -1.0 * rp[0] / (pcnt * EmitH); // they are identical!
  double BetaH = r2[0] / (pcnt * EmitH);
  double GammaH = p2[0] / (pcnt * EmitH);
    
  double AlphaV = -1.0 * rp[1] / (pcnt * EmitV); // they are identical!
  double BetaV = r2[1] / (pcnt * EmitV);
  double GammaV = p2[1] / (pcnt * EmitV);

  
  /*
   * Calculate RMS quantities
   * ========================
   */
  
    double xrms   = sqrt(r2[0]/pcnt-sqr(rm[0]/pcnt)); 
    double dxrms   = xrms/sqrt(pcnt);
    double xprms  = sqrt(p2[0]/pcnt-sqr(pm[0]/pcnt)); 
    double dxprms  = xprms/sqrt(pcnt);
 
    
    double yrms   = sqrt(r2[1]/pcnt-sqr(rm[1]/pcnt)); 
    double dyrms   = xrms/sqrt(pcnt);
    double yprms  = sqrt(p2[1]/pcnt-sqr(pm[1]/pcnt)); 
    double dyprms  = xprms/sqrt(pcnt);
     
    double trms   = sqrt(r2[2]/pcnt-sqr(rm[2]/pcnt)); 
    double dtrms   = xrms/sqrt(pcnt);
    double tprms  = sqrt(p2[2]/pcnt-sqr(pm[2]/pcnt)); 
    double dtprms  = xprms/sqrt(pcnt);
 

    double dxxprms,dyyprms,dttprms;
    RmsUncertainty(CxxH,CxyH,CyyH,EmatH,&dxrms,&dxprms,&dxxprms);
    RmsUncertainty(CxxV,CxyV,CyyV,EmatV,&dyrms,&dyprms,&dyyprms);
    RmsUncertainty(CxxT,CxyT,CyyT,EmatT,&dtrms,&dtprms,&dttprms);
  

    double dAlphaH = 0.0;
    double dBetaH  = 0.0;
    double dGammaH = 0.0;
    double dAlphaV = 0.0;
    double dBetaV  = 0.0;
    double dGammaV = 0.0;
  
    TwissUncertainty(CxxH,CxyH,CyyH, EmatH, &dAlphaH, &dBetaH, &dGammaH);
    TwissUncertainty(CxxV,CxyV,CyyV, EmatV, &dAlphaV, &dBetaV, &dGammaV);
    
    double dEmittH = EmittanceUncertainty(CxxH, CxyH, CyyH, EmatH);
    double dEmittV = EmittanceUncertainty(CxxV, CxyV, CyyV, EmatV);

    INFOMSG("AlphaH " << AlphaH << " +- " << dAlphaH << endl);
    INFOMSG("BetaH  " << BetaH << "  +- " << dBetaH  << endl);
    INFOMSG("GammaH " << GammaH << " +- " << GammaH << endl << endl);

    INFOMSG("AlphaV " << AlphaV << " +- " << dAlphaV << endl);
    INFOMSG("BetaV  " << BetaV << "  +- " << dBetaV  << endl);
    INFOMSG("GammaV " << GammaV << " +- " << GammaV << endl << endl);
    
    INFOMSG("EmitH " << EmitH << " +- " << dEmittH << endl);
    INFOMSG("EmitV " << EmitV << " +- " << dEmittV << endl);
}

// Added by H.Krueger
double Distribution::gauss(double av, double sd, double ran, double ran2) 
       {return av+sd*sqrt(-2*log(ran))*cos(6.283185307179586*ran2);}						     







