// ------------------------------------------------------------------------
// $RCSfile: Distribution.cpp,v $
// $Header: /afs/psi.ch/user/a/adelmann/private/cvsroot/pcyclint/prog/Distribution/Distribution.cpp,v 1.3 2004/10/01 20:33:00 adelmann Exp $
// ------------------------------------------------------------------------
// $Revision: 1.3 $
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
// $Date: 2004/10/01 20:33:00 $
// $Author: adelmann $
// $Log: Distribution.cpp,v $
// Revision 1.3  2004/10/01 20:33:00  adelmann
// add GAUSS and UNIFORM distribution
//
// Revision 1.2  2004/09/24 19:19:05  adelmann
// - Add stuff for Gaussian Distribution
//
// Revision 1.1.1.1  2004/09/22 12:10:45  adelmann
// Imported sources pcyclubt
//
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

Distribution::~Distribution()
{
  
}

Distribution::Distribution(ChargedParticles<double,3> *beam, DistrData *distrData, double len) :
    beam_m(beam),
    distrData_m(distrData)
  
{

    Inform msg("Create distribution  ");

    msg << "Create distribution not implemented yet" << endl;
    
}

Distribution::Distribution(ChargedParticles<double,3> *beam) :
  beam_m(beam)  
{
/*
  Vector_t emit(distrData_m->emitx,distrData_m->emity,distrData_m->emitt);
  Vector_t alpha;
  Vector_t beta;
  Vector_t bincoef = Vector_t(9999999.,9999999.,9999999.); // gauss
    
  beta[0] =  sqr(beam->actSigma[0])/beam->emit[0];
  beta[1] =  sqr(beam->actSigma[3])/beam->emit[1];
  beta[2] =  sqr(beam->actSigma[6])/beam->emit[2];
	
  alpha[0] = -beam->actSigma[2];
  alpha[1] = -beam->actSigma[5];
  alpha[2] = -beam->actSigma[8];
	
  regenerateBinom(beam->getTotalNum(),emit,alpha,beta,bincoef);
  
  beam_m->rescaleAndUpdate();
  beam_m->qm = 1.0;
  beam_m->mass = beam_m->getMassKg();
  beam_m->elemNo = 0; // start in first element
*/
}

Distribution::Distribution(ChargedParticles<double,3> *beam, DistrData *distrData, string fn) :
  beam_m(beam),
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

Vector_t Distribution::generateGaussian(
    double qTot, 
    double mTot, 
    unsigned long int particles,
    Vector_t xmean, 
    Vector_t pmean, 
    Vector_t xstddev, 
    Vector_t pstddev, 
    Vector_t angle,
    double ekin, 
    double r, 
    double alpha, 
    double vvertical,
    double PR)
{
    Inform msg2all("generateGaussianDistr ", INFORM_ALL_NODES);
    Inform msg("generateGaussianDistr ");

    unsigned int p = 0;
    unsigned int k = 0;      
    double q = qTot/particles;
    double m = mTot/particles;
    
    RANLIB_class *rGen = new RANLIB_class(265314159,4);

    Vector_t dr,dv;  // deviations from r and v
    Vector_t rv;  //  
        
    double ptot= sqrt(pow((ekin/MPMEV+1.0),2) -1.0); //TOTAL MOMENTUM IN CU
    double pr = PR /1e3;                              //FIXPO MCU to CU
    double pt = sqrt(ptot*ptot-pr*pr);
    double alphav = 0.0;
        
    rv = Vector_t(r*cos(alpha),r*sin(alpha),0.0);


    if (beam_m->isRoot()) {
        for(unsigned long int i=0; i<particles; i++) {   
            double dx   = rGen->gauss(xmean(0),sqrt(xstddev(0)));
            double dpx  = rGen->gauss(pmean(0),sqrt(pstddev(0)));
            double dy   = rGen->gauss(xmean(1),sqrt(xstddev(1)));
            double dpy  = rGen->gauss(pmean(1),sqrt(pstddev(1)));
            double dz   = rGen->gauss(xmean(2),sqrt(xstddev(2)));
            double dpz  = rGen->gauss(pmean(2),sqrt(pstddev(2)));
            
            double x0   = dx * cos(angle(0)) - dpx * sin(angle(0));
            double px0  = dx * sin(angle(0)) + dpx * cos(angle(0));
            
            double y0   = dy * cos(angle(1)) - dpy * sin(angle(1));
            double py0  = dy * sin(angle(1)) + dpy * cos(angle(1));
            
            double z0   = dz * cos(angle(2)) - dpz * sin(angle(2));
            double pz0  = dz * sin(angle(2)) + dpz * cos(angle(2));
            
            beam_m->create(1);
          
            beam_m->R[k] = rv + Vector_t(x0,y0,z0)*xstddev;
            double npr = pr + px0*pstddev[0];
            double npt = sqrt(ptot*ptot-npr*npr);
            beam_m->P[k] = Vector_t(
                beam_m->cu2mks(npr*cos(alpha)-npt*sin(alpha)),
                beam_m->cu2mks(npr*sin(alpha)+npt*cos(alpha)),
                vvertical*pz0
                );
            beam_m->mass_m[k] = m;
            beam_m->q_m[k] = q;
            beam_m->bunchNo[k]= 0;
            k++;
        }
    }
    msg2all << "Particles created: " << beam_m->getLocalNum() << endl;
    beam_m->boundp();
    msg2all << "Particles after boundp: " << beam_m->getLocalNum() << endl;
    Vector_t cm = sum(beam_m->R);
    cm /= beam_m->getTotalNum();
    return cm;
}
    
Vector_t  Distribution::generateUniform(
    double qTot, 
    double mTot, 
    unsigned long totalP, 
    Vector_t sigmar, 
    Vector_t sigmav, 
    double ekin, 
    double r, 
    double alpha, 
    double vvertical,
    double PR) {
        
    Inform msg2all("generateUniform ", INFORM_ALL_NODES);
    Inform msg("generateUniform ");

    unsigned k = 0;
    
    double q = qTot/totalP;
    double m = mTot/totalP;
    
    Vector_t dr,dv;  // deviations from r and v
    Vector_t rv;  //  
        
    double ptot= sqrt(pow((ekin/MPMEV+1.0),2) -1.0); //TOTAL MOMENTUM IN CU
    double pr = PR /1e3;                              //FIXPO MCU to CU
    double pt = sqrt(ptot*ptot-pr*pr);

    rv = Vector_t(r*cos(alpha),r*sin(alpha),0.0);
    
    if (beam_m->isRoot()) {
        for (int n=0; n<totalP; n++) {
            do {
                dr = Vector_t(2.0*IpplRandom() - 1.0,2.0*IpplRandom() - 1.0,2.0*IpplRandom() - 1.0);
            } while (dot(dr,dr)>1.0);
            
            do {
                dv = Vector_t(2.0*IpplRandom() - 1.0,2.0*IpplRandom() - 1.0,2.0*IpplRandom() - 1.0);
            } while (dot(dv,dv)>1.0);
            
            beam_m->create(1);                
            
            beam_m->R[k] = rv + dr*sigmar;   // rv +- delta
            double npr = pr + dv[0]*sigmav[0];
            double npt = sqrt(ptot*ptot-npr*npr);
            beam_m->P[k] = Vector_t(
                beam_m->cu2mks(npr*cos(alpha)-npt*sin(alpha)),
                beam_m->cu2mks(npr*sin(alpha)+npt*cos(alpha)),
                vvertical*sigmav[2]);
            beam_m->mass_m[k] = m;
            beam_m->q_m[k] = q;
            beam_m->bunchNo[k]= 0;
            k++;
        }
    }
    msg2all << "Particles created: " << beam_m->getLocalNum() << endl;
    beam_m->boundp();
    msg2all << "Particles after boundp: " << beam_m->getLocalNum() << endl;
    Vector_t cm = sum(beam_m->R);
    cm /= beam_m->getTotalNum();
    return cm;
} 

Vector_t Distribution::readInputDistribution(string fn) {
    double  arg1, arg2, arg3, arg4, arg5, arg6;

    Inform msg2all("readInputDistribution ", INFORM_ALL_NODES);
    Inform msg("readInputDistribution ");
    
    ifstream *ifstr = new ifstream;
    
    if (beam_m->isRoot()) {
        ifstr->open(fn.c_str(), ios::in );
        unsigned int count = 0;
    
        *ifstr >> arg1 >> arg2 >> arg3 >> arg4 >> arg5 >> arg6;
        
        while(*ifstr) {
            beam_m->create(1);
            beam_m->R[count](0) = arg1;
            beam_m->R[count](1) = arg3;
            beam_m->R[count](2) = arg5;
            beam_m->P[count](0) = arg2;
            beam_m->P[count](1) = arg4;
            beam_m->P[count](2) = arg6;
            beam_m->bunchNo[count]= 0;
            count++;
            *ifstr >> arg1 >> arg2 >> arg3 >> arg4 >> arg5 >> arg6;
        }
    }
    msg2all << "Particles created: " << beam_m->getLocalNum() << endl;
    beam_m->boundp();
    msg2all << "Particles after boundp: " << beam_m->getLocalNum() << endl;
    Vector_t cm = sum(beam_m->R);
    cm /= beam_m->getTotalNum();
    return cm;
}
    
 
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


void Distribution::createBinom(unsigned long int particles,
				      Vector_t emit,
				      Vector_t alpha,
				      Vector_t beta,
				      Vector_t bincoef)
{ 

  const double two_pi = 8.0*atan(1.0);
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

    const double two_pi = 8.0*atan(1.0);
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















