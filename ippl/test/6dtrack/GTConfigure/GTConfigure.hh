#ifndef CONFIGURE_HH
#define CONFIGURE_HH
// ------------------------------------------------------------------------
// $RCSfile: GTConfigure.hh,v $
// $Header: /afs/psi.ch/user/a/adelmann/private/cvsroot/GenTrackE/GTConfigure/GTConfigure.hh,v 1.7 2004/04/28 08:04:57 adelmann Exp $
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
// $Date: 2004/04/28 08:04:57 $
// $Author: adelmann $
// $Log: GTConfigure.hh,v $
// Revision 1.7  2004/04/28 08:04:57  adelmann
// Defines MAPS and cut out many map related stuff in order
// to make SEABORG happy
//
// Revision 1.6  2004/04/06 13:21:45  adelmann
// *** empty log message ***
//
// Revision 1.5  2004/04/05 13:07:06  adelmann
// Many changes to use new IPPL instead of POOMA
//
// Revision 1.4  2003/05/02 13:57:39  adelmann
// First electron cloud simulation 3d without SC and probable still with the wrong
// drive beam integration scheme. Otherwise :=) the program runs
// sureprisingly stable and the first inspection of the data looks promissing
//
// Revision 1.3  2003/04/17 14:18:43  adelmann
// *** empty log message ***
//
// Revision 1.2  2003/01/24 15:23:42  adelmann
// Have to define DBL_EPSILON by hand. On the pc2836
// I do not have to to that ??
//
// Revision 1.1.1.1  2003/01/23 09:13:52  adelmann
// GenTrackE
//
// ------------------------------------------------------------------------
#include "GTConst.hh"
#include "GTElemData.hh"

#ifdef MAPS
#include "GTSbendMap.hh"
#include "GTDriftMap.hh"
#include "GTMultipoleMap.hh"
#include "GTRecquadMap.hh"
//#include "GTSymplectifyLinear.hh"
#include "GTRFCavityMap.hh"
#include "GTIBS.hh"
#endif

#include "Algorithms/PartData.h"

#include "MadParser/MadParser.h"
#include "Parser/FileStream.h"

#include <iostream>
#include <fstream>

#include <new>
#include <exception>

#include "Ippl.h"

/* Prototyles */
ElemData getElementData(FileStream *is, Token tk);

bool isConfigureSection (FileStream *is) {

    Token token = is->readToken();
    
    if (token.isKey("CONFIGURE")) {
	token = is->readToken();
	if (token.isKey("BEGIN")) {
	    return true;
	}
	else { 
	    return false;
	}
    }
    else {
	return false;	
    }
}

bool getToken (FileStream *is, Token *token) {
  // read until END is reached
  *token = is->readToken();
  return (!token->isKey("END"));
}

bool isLatticeSection (FileStream *is) {

    Token token = is->readToken();
    if (token.isKey("LATTICE")) {
	token = is->readToken();
	if (token.isKey("BEGIN")) {
	    return true;
	}
	else { 
	    return false;
	}
    }
    else {
	return false;	
    }
}

/**
 * Read in simulation configuration from inputfile 
 */
bool FSConfigure(FileStream *is, SimCfgData *cfgData, DistrData *distrData)
{
  string bc_str;
  string interPol_str;
  string dist_str;

  cfgData->Qtot = 0.0;
  cfgData->trackSymplectic=false;
  cfgData->deltaEPS=-1.0;
  cfgData->timeinvers=false;
  cfgData->dt=-1;
  cfgData->ndt = 10;
  cfgData->distrInputFile=string("");
  cfgData->I0 = -1.0;

  cfgData->doIBS = false;
  cfgData->dipoleIBSKicks = 1;

  cfgData->doRadiation = false;

  cfgData->latticeTest = false;

  cfgData->doTouschek = false;

  cfgData->interrad = 0;
  cfgData->meshmultfactor = 2;

  cfgData->diagnostics = false;
  cfgData->doRestart = false;
  cfgData->doInbalance = false;
  cfgData->doPhaseSpace = false;
  cfgData->doStatistics = false;
  cfgData->doLostPlace = false;
  cfgData->doAllInteractions = false;

  cfgData->scattering = 1;
  cfgData->distributiontype = 1;

  if (isConfigureSection (is)) {
    Token tk;
    getToken (is, &tk);
    while (!tk.isKey("END")) {
      if (tk.isKey("DIPOLEIBSKICKS")) {
	getToken (is, &tk);
	cfgData->dipoleIBSKicks =  tk.getInteger();
      } else if (tk.isKey("RADIATION")) {
	cfgData->doRadiation = true;
	getToken (is, &tk);
	if (tk.isKey("QUANTUM"))
	  cfgData->radiationProcess = QUANTUM;
	else
	  cfgData->radiationProcess = CLASSIC;
      } else if (tk.isKey("IBS")) {
	cfgData->doIBS = true;
      } else if (tk.isKey("DT")) {
	getToken (is, &tk);
	cfgData->dt = tk.getReal();
      } else if (tk.isKey("NDT")) {
	getToken (is, &tk);
	cfgData->ndt = tk.getInteger();
      } else  if (tk.isKey("TIMEINVERS")) {
	cfgData->timeinvers = true;
      } else  if (tk.isKey("TITLE")) {
	getToken (is, &tk);
	cfgData->title = tk.getString();
      } else if (tk.isKey("I" )) {
	getToken (is, &tk);
	cfgData->I0 = tk.getReal();
      } else if (tk.isKey("EK" )) {
	getToken (is, &tk);
	cfgData->Ek = tk.getReal();
      } else if (tk.isKey("FREQU")) {
	getToken (is, &tk);
	cfgData->frequ = tk.getReal();
      } else if (tk.isKey("PERIODS" )) {
	getToken (is, &tk);
	cfgData->periods =  tk.getInteger();
      } else if (tk.isKey("ORDER")) {
	getToken (is, &tk);
	cfgData->order = tk.getInteger();
      } else if (tk.isKey("GRID")) {
	getToken (is, &tk);
	cfgData->nx = tk.getInteger();
	getToken (is, &tk);
	cfgData->ny = tk.getInteger();
	getToken (is, &tk);
	cfgData->nz = tk.getInteger();
      } else if (tk.isKey("CENTER")) {
	getToken (is, &tk);
	cfgData->centerx = tk.getReal();
	getToken (is, &tk);
	cfgData->centery = tk.getReal();
	getToken (is, &tk);
	cfgData->centert = tk.getReal();
      } else if (tk.isKey("NPINIT")) {
	getToken (is, &tk);
	cfgData->nPInit = tk.getInteger();
      } else if (tk.isKey("NEINIT")) {
	getToken (is, &tk);
	cfgData->nEInit = tk.getInteger();
      } 
      else if (tk.isKey("SUBTSTEPS")) {
	getToken (is, &tk);
	cfgData->subTSteps = tk.getInteger();
      } 
      else if (tk.isKey("PARTICLESTOSAVE")) {
	getToken (is, &tk);
	cfgData->particlesToSave = tk.getInteger();
      } else if (tk.isKey("QTOT")) {
	getToken (is, &tk);
	cfgData->Qtot = tk.getReal();
      }  else if (tk.isKey("USEQUAD")) {
	getToken (is, &tk);
	cfgData->Usequad =  tk.getBool();
      } else if (tk.isKey("LATTICETEST")) {
	getToken (is, &tk);
	cfgData->latticeTest = tk.getBool();
      } else if (tk.isKey("ELCLINPUT")) {
	getToken (is, &tk);
	cfgData->elclInputFile = tk.getString();
      } else if (tk.isKey("DIAGFILE")) {
	getToken (is, &tk);
	cfgData->diagfile = tk.getString();
      } else if (tk.isKey("DATADIR")) {
	getToken (is, &tk);
	cfgData->datadir = tk.getString();
      } else if (tk.isKey("INTERACTION")) {
	getToken (is, &tk);
	if (tk.getString()=="GRAVITY")
	  cfgData->interaction=GRAVITY;
	else if (tk.getString()=="COULOMB")
	  cfgData->interaction=COULOMB;
	else if (tk.getString()=="ELECTRONCLOUD")
	  cfgData->interaction=ELECTRONCLOUD;
	else
	  cfgData->interaction=COULOMB;
      } else if (tk.isKey("GRIDSCHEME")) {
	getToken (is, &tk);
	if (tk.getString()=="NGP") 
	  cfgData->interPol = NGP;
	else 
	  cfgData->interPol = CIC;
      }  else if (tk.isKey("SOLVER")) {
	getToken (is, &tk);
	if (tk.getString()=="P2P") {
	  cfgData->fsType = P2PSolver;
	} else if (tk.getString()=="FFTOOO") {
	  cfgData->fsType = FFTSolverOOO;
	} else if (tk.getString()=="FFTOOP") {
	  cfgData->fsType = FFTSolverOOP;
	} else if (tk.getString()=="TREE") {
	  cfgData->fsType = TreeSolver;
	} else if (tk.getString()=="MULTIGRIDFEM") {
	  cfgData->fsType = MultigridFEM;
	} else if (tk.getString()=="MGAMRSOLVER") {
	  cfgData->fsType = MGAMRSolver;
	} else 
	  cfgData->fsType = None;
      } else if (tk.isKey("BC")) {
	  getToken (is, &tk);
	  if (tk.getString()=="OOO") {
	      cfgData->bc = OOO;
	  } else if (tk.getString()=="OOP") {
	      cfgData->bc = OOP;
	  } else if (tk.getString()=="NEUMANN") {
	      cfgData->bc = NEUMANN;
	  } else if (tk.getString()=="DIRICHLET") {
	      cfgData->bc = DIRICHLET;
	  } else if (tk.getString()=="INHOMDIRICHLET") {
	      cfgData->bc = INHOMDIRICHLET;
	  } 
	  else 
	      cfgData->bc = OOO;
      } else if (tk.isKey("UNITS")) {
	getToken (is, &tk);
	if (tk.getString()=="MARYLIE"){
	  cfgData->uType=MARYLIE;
	}
	else if (tk.getString()=="ELCL"){
	  cfgData->uType=ELCL;
	}
	else {
	  cfgData->uType=MAD9;
	}
      } else if (tk.isKey("TRACKSYMPLECTIC")){
	cfgData->trackSymplectic=true;
      } else if (tk.isKey("TRACKOneTurnMap")) {
        cfgData->trackOneTurnMap=true;
      } else if (tk.isKey("DISTRTYPE")) {
	getToken (is, &tk);
	if (tk.isKey("ELLIPSOIDALUNIFORM"))
	  cfgData->distrType=ELLIPSOIDALUNIFORM;
	else if (tk.isKey("READFROMFILE")) {
	  cfgData->distrType=READFROMFILE;
	  getToken (is, &tk);
	  cfgData->distrInputFile=tk.getString();
	}
	else 
	  cfgData->distrType=BINOMINAL;
      }  else if (tk.isKey("SIGX")) {
	getToken (is, &tk);
	distrData->sigx = tk.getReal();
      } else if (tk.isKey("SIGPX")) {
	getToken (is, &tk);
	distrData->sigpx = tk.getReal();
      } else if (tk.isKey("CORXPX")) {
	getToken (is, &tk);
	distrData->corxpx = tk.getReal();
      } else if (tk.isKey("MX")) {
	getToken (is, &tk);
	distrData->mx = tk.getReal();
      } else if (tk.isKey("EMITX")) {
	getToken (is, &tk);
	distrData->emitx = tk.getReal();
      } else if (tk.isKey("SIGY")) {
	getToken (is, &tk);
	distrData->sigy = tk.getReal();
      } else if (tk.isKey("SIGPY")) {
	getToken (is, &tk);
	distrData->sigpy = tk.getReal();
      } else if (tk.isKey("CORYPY")) {
	getToken (is, &tk);
	distrData->corypy = tk.getReal();
      } else if (tk.isKey("MY")) {
	getToken (is, &tk);
	distrData->my = tk.getReal();
      } else if (tk.isKey("EMITY")) {
	getToken (is, &tk);
	distrData->emity = tk.getReal();
      } else if (tk.isKey("SIGT")) {
	getToken (is, &tk);
	distrData->sigt = tk.getReal();
      } else if (tk.isKey("SIGPT")) {
	getToken (is, &tk);
	distrData->sigpt = tk.getReal();
      } else if (tk.isKey("CORTPT")) {
	getToken (is, &tk);
	distrData->cortpt = tk.getReal();
      } else if (tk.isKey("MT")) {
	getToken (is, &tk);
	distrData->mt = tk.getReal();
      } else if (tk.isKey("EMITT")) {
	getToken (is, &tk);
	distrData->emitt = tk.getReal();
      } else if (tk.isKey("DELTAEPS")) {
	getToken (is, &tk);
#define DBL_EPSILON 2.2e-16
	cfgData->deltaEPS = tk.getReal()*DBL_EPSILON;
      } else if (tk.isKey("TIME")) {
	cfgData->indepVar = TIME;
	getToken (is, &tk);
	if (tk.isKey("VERLET"))
	  cfgData->integrator = VERLET;
	else if (tk.isKey("RK4"))
	  cfgData->integrator = RK4;
	else
	  cfgData->integrator = LEAPFROG;
      } else if (tk.isKey("ARCLENGTH")) {
          cfgData->indepVar = ARCLEN;
          cfgData->integrator = SPLIT;
      } else if (tk.isKey("PROTON")) {
          cfgData->specis = PROTONS;
      } else if (tk.isKey("ELECTRON")) {
          cfgData->specis = ELECTRONS;
      } else if (tk.isKey("HMINUS")) {
          cfgData->specis = HMINUS;
      } else if (tk.isKey("MAXMGLEVEL")) {
          getToken (is, &tk);
          cfgData->maxMgLevel = tk.getInteger();
      } else if (tk.isKey("GRIDSIZE")) {
          getToken (is, &tk);
          cfgData->hInit = tk.getReal();
      } else if (tk.isKey("RELERROR")) {
          getToken (is, &tk);
          cfgData->relErr = tk.getReal();
      } else if (tk.isKey("RELAX")) {
          getToken (is, &tk);
          cfgData->iter_relax = tk.getInteger();
      } else if (tk.isKey("GRIDSTRETCH")) {
          getToken (is, &tk);
          cfgData->grid_stretch = tk.getReal();
      } else if (tk.isKey("LAMBDAP")) {
          getToken (is, &tk);
          cfgData->lambdaP = tk.getReal();
      } else if (tk.isKey("LAMBDAE")) {
          getToken (is, &tk);
          cfgData->lambdaE = tk.getReal();
      } else if (tk.isKey("NEUTRALFAC")) {
	getToken (is, &tk);
	cfgData->neutralfac = tk.getReal();
      } else if (tk.isKey("TOUSCHEK")) {
	cfgData->doTouschek = true;
      } else if (tk.isKey("DIVISION")) {
	getToken (is, &tk);
	cfgData->Division = tk.getReal();
      } else if (tk.isKey("DISPERSION")) {
	getToken (is, &tk);
	cfgData->dispersion = tk.getReal();
      } else if (tk.isKey("RADIUS")) {
	getToken (is, &tk);
	cfgData->radius = tk.getReal();
      } else if (tk.isKey("FIELDINDEX")) {
	getToken (is, &tk);
	cfgData->fieldindex = tk.getReal();
      } else if (tk.isKey("LAMBDA")) {
	getToken (is, &tk);
	cfgData->lambda = tk.getReal();
      } else if (tk.isKey("ACCPX")) {
	getToken (is, &tk);
	cfgData->accpx = tk.getReal();
      } else if (tk.isKey("ACCPY")) {
	getToken (is, &tk);
	cfgData->accpy = tk.getReal();
      } else if (tk.isKey("ACCPS")) {
	getToken (is, &tk);
	cfgData->accps = tk.getReal();
      } else if (tk.isKey("INTERRAD")) {
	getToken (is, &tk);
	cfgData->interrad = tk.getReal();
      } else if (tk.isKey("MESHMULTFACTOR")) {
	getToken (is, &tk);
	cfgData->meshmultfactor = tk.getInteger();
      } else if (tk.isKey("DIAGNOSTICS")) {
	cfgData->diagnostics = true;
      } else if (tk.isKey("DORESTART")) {
	cfgData->doRestart = true;
      } else if (tk.isKey("DOSTATISTICS")) {
	cfgData->doStatistics = true;
      } else if (tk.isKey("DOINBALANCE")) {
	cfgData->doInbalance = true;
      } else if (tk.isKey("DOPHASESPACE")) {
	cfgData->doPhaseSpace = true;
      } else if (tk.isKey("SCATTERING")) {
	getToken(is, &tk);
	cfgData->scattering = tk.getInteger();
      } else if (tk.isKey("DISTRIBUTIONTYPE")) {
	getToken(is, &tk);
	cfgData->distributiontype = tk.getInteger();
      } else if (tk.isKey("DOLOSTPLACE")) {
	cfgData->doLostPlace = true;
      } else if (tk.isKey("DOALLINTERACTIONS")) {
	cfgData->doAllInteractions = true;
      } else if (tk.isKey("DOCROSSSECTION")) {
	cfgData->doCrosssection = true;
      } 
	getToken (is, &tk);
    }     
    // update distrData
    distrData->N=cfgData->nPInit;
    distrData->distrType = cfgData->distrType;
    return tk.isKey("END");
  }
  else {
    return false;
  }
}

#ifdef MAPS
void simplectify (FMatrix<double,2*DIMENSION,2*DIMENSION> &m) {
  /*
    Based on F. Neris ML implementation
  */
  unsigned int d = 2*DIMENSION;
  
  double  qq,pq,qp,pp;
  unsigned int kp,kq,lp,lq,jp,jq,i;
  
  for (kp=1; kp<d;kp+=2) {                                  //do 100 kp=2,2*n,2
    kq = kp-1;
    for (lp=1; kp-1;lp+=2) {                                //do 200 lp=2,kp-2,2
      lq = lp-1;
      qq = 0.0;
      pq = 0.0;
      qp = 0.0;
      pp = 0.0;
      for (jp=1; jp<d;jp+=2) {   	                    //do 300 jp=2,2*n,2
	jq = jp-1;
	qq = qq + m(lq,jq)*m(kq,jp) - m(lq,jp)*m(kq,jq);
	pq = pq + m(lp,jq)*m(kq,jp) - m(lp,jp)*m(kq,jq);
	qp = qp + m(lq,jq)*m(kp,jp) - m(lq,jp)*m(kp,jq);
	pp = pp + m(lp,jq)*m(kp,jp) - m(lp,jp)*m(kp,jq);
      }                                                     //300 continue
      for (i=0; i<d;i++) {                                  //do 400 i=1,2*n
	m(kq,i) = m(kq,i) - qq*m(lp,i) + pq*m(lq,i);
	m(kp,i) = m(kp,i) - qp*m(lp,i) + pp*m(lq,i);
      }                                                     //400 continue
    }                                                       //200 continue
    qp = 0.0;
    for (jp=1; jp<d;jp+=2) {                                //do 500 jp=2,2*n,2
      jq = jp-1;
      qp = qp + m(kq,jq)*m(kp,jp) - m(kq,jp)*m(kp,jq);
    }                                                       //500 continue
    for (i=0; i<d;i++)                                      //do 600 i=1,2*n
      m(kp,i) = m(kp,i)/qp;                                 //600     continue
  }                                                         //100  continue
}
#endif


double buildLattice (FileStream *is, vector< ElMap > &lattice,
		     SimCfgData &simCfg, PartData &pdata)
{  
  Inform msg("buildLattice ");  
  Token tk;
  
  double eps = 1.0e-10;
  
  double actualS = 0.0;
  double entrance = 0.0;
  double exit = 0.0;
  
  int elemNo = -1;  // in order to start with zero
  
  if (isLatticeSection (is)) {
    getToken (is, &tk);
    while (!tk.isKey("END")) { 
      
      ElMap newElement;
      ElemData elemData = getElementData(is,tk);
      
      /*
	number all electromagnetic elements
	needed for indirectionlist
      */
      if ((elemData.elType != STATMARKER) &&
	  (elemData.elType != PARTMARKER)) {
	elemNo++;
	elemData.elemNo = elemNo;
      }
      

      elemData.order = simCfg.order;
      elemData.uType = simCfg.uType;
      newElement.elemData = elemData;
      
#ifdef MAPS      
      FMatrix<double,2*DIMENSION,2*DIMENSION> m;  // linear map
      FVpsT M,C;                                  // Taylor series
      //      DragtFinnMap<DIMENSION> dfMap;
#endif      

      double L;
      
      L = elemData.length;
      newElement.elemData.length = L;

      if (elemData.elType == CONSTFOC) {
	actualS += L;
	newElement.actualS=actualS;
	/*
	  used in the t code 
	*/
	entrance = exit;
	exit = entrance + L;
	newElement.elemData.entrance = entrance;
	
	lattice.push_back(newElement);
      }

      /*
        Build a DRIFT
      */
      if (elemData.elType == DRIFT) {
	
	actualS += L;
	newElement.actualS=actualS;
        /*
          used in the t code
	*/
        entrance = exit;
        exit = entrance + L;
        newElement.elemData.entrance = entrance;
#ifdef MAPS      	
	DriftMap *mDrift = new DriftMap(newElement.elemData,pdata);
	M = mDrift->getTaylorSeries(M);
	m = M.linearTerms();
	newElement.Msymp.push_back(m);
	newElement.taylor.push_back(M);
	delete(mDrift);
#endif
        lattice.push_back(newElement);
      }
      
      /*
        Build a DIPOLE
      */
      if (elemData.elType == DIPOLE) {
        actualS += L;
        newElement.actualS=actualS;
        /*
          used in the t code
	*/
        entrance = exit;
        exit = entrance + L;
        newElement.elemData.entrance = entrance;
	
        lattice.push_back(newElement);
      }
      
      /*
        Build a QUADRUPOLE
      */
      if (elemData.elType == QUADRUPOLE) {

        actualS += L;
        newElement.actualS=actualS;
        /*
          used in the t code
	*/
        entrance = exit;
        exit = entrance + L;
        newElement.elemData.entrance = entrance;
#ifdef MAPS      		
	MultiPoleMap  *mQuad  = new MultiPoleMap(newElement.elemData,pdata);
	M = mQuad->getTaylorSeries(M);
	m = M.linearTerms();
	newElement.Msymp.push_back(m);
	newElement.taylor.push_back(M);
        delete(mQuad);
#endif	
        lattice.push_back(newElement);
      }
      getToken (is, &tk);
    }   
  }
  
  /*
    update lastElIdx;
  */
  vector< ElMap >::iterator latit = lattice.begin();
  for( ; latit != lattice.end(); latit++) 
    (*latit).lastElIdx = elemNo;
  
  double len;
  len = actualS;
  return len;
}
  
  
  void get2Token (FileStream *is, Token *tk)
 {
    getToken (is, tk);
    getToken (is, tk);
}

ElemData getElementData(FileStream *is, Token tk)
{
  Inform msg("getElementData ");  
    /*
      Assume token is set to valid elementtype or END 
     */ 
    
    ElemData newElement;
      
    if (tk.isKey("DRIFT")) {
      newElement.elType = DRIFT;
      get2Token (is, &tk);
      newElement.length=tk.getReal();
      get2Token (is, &tk);
      newElement.a=tk.getReal();
      get2Token (is, &tk);
      newElement.b=tk.getReal();
      get2Token (is, &tk);
      newElement.missX=tk.getReal();
      get2Token (is, &tk);
      newElement.missY=tk.getReal();
      
    }
    else if (tk.isKey("DIPOLE")) {
      newElement.elType = DIPOLE;
      get2Token (is, &tk);
      newElement.length=tk.getReal();
      get2Token (is, &tk);
      newElement.K0=tk.getReal();
      get2Token (is, &tk);
      newElement.a=tk.getReal();
      get2Token (is, &tk);
      newElement.b=tk.getReal();
      get2Token (is, &tk);
      newElement.missX=tk.getReal();
      get2Token (is, &tk);
      newElement.missY=tk.getReal();
    }
    else if (tk.isKey("SBEND")) {
      newElement.elType = SBEND;
      get2Token (is, &tk);            // -length
      newElement.length = tk.getReal();
      get2Token (is, &tk);            // e1
      newElement.e1 = tk.getReal();
      get2Token (is, &tk);            // e2
      newElement.e2 = tk.getReal();
      get2Token (is, &tk);            // angle
      newElement.angle = tk.getReal();
      get2Token (is, &tk);            // hgap
      newElement.hgap = tk.getReal();
      get2Token (is, &tk);            // fint
      newElement.fint = tk.getReal();
      get2Token (is, &tk);            // K0
      newElement.K0 = tk.getReal();
      get2Token (is, &tk);            // K0S
      newElement.K0S = tk.getReal();
      get2Token (is, &tk);            // K1
      newElement.K1 = tk.getReal();
      get2Token (is, &tk);            // K1S
      newElement.K1S = tk.getReal();
      get2Token (is, &tk);            // K2
      newElement.K2 = tk.getReal();
      get2Token (is, &tk);            // K2S
      newElement.K2S = tk.getReal();
      get2Token (is, &tk);            // K3
      newElement.K3 = tk.getReal();
      get2Token (is, &tk);            // K3S
      newElement.K3S = tk.getReal();
      get2Token (is, &tk);            // a
      newElement.a=tk.getReal();
      get2Token (is, &tk);            // b
      newElement.b=tk.getReal();
      get2Token (is, &tk);
      newElement.missX=tk.getReal();
      get2Token (is, &tk);
      newElement.missY=tk.getReal();
    } 
    else if (tk.isKey("VSBEND")) {
      newElement.elType = VSBEND;
      get2Token (is, &tk);            // -length
      newElement.length = tk.getReal();
      get2Token (is, &tk);            // e1
      newElement.e1 = tk.getReal();
      get2Token (is, &tk);            // e2
      newElement.e2 = tk.getReal();
      get2Token (is, &tk);            // angle
      newElement.angle = tk.getReal();
      get2Token (is, &tk);            // hgap
      newElement.hgap = tk.getReal();
      get2Token (is, &tk);            // fint
      newElement.fint = tk.getReal();
      get2Token (is, &tk);            // K0
      newElement.K0 = tk.getReal();
      get2Token (is, &tk);            // K0S
      newElement.K0S = tk.getReal();
      get2Token (is, &tk);            // K0
      newElement.K1 = tk.getReal();
      get2Token (is, &tk);            // K0S
      newElement.K1S = tk.getReal();
      get2Token (is, &tk);            // K0
      newElement.K2 = tk.getReal();
      get2Token (is, &tk);            // K0S
      newElement.K2S = tk.getReal();
      get2Token (is, &tk);            // K0
      newElement.K3 = tk.getReal();
      get2Token (is, &tk);            // K0S
      newElement.K3S = tk.getReal();
      get2Token (is, &tk);
      newElement.a=tk.getReal();
      get2Token (is, &tk);
      newElement.b=tk.getReal();
      get2Token (is, &tk);
      newElement.missX=tk.getReal();
      get2Token (is, &tk);
      newElement.missY=tk.getReal();
    } 
    else if  (tk.isKey("VSBEND")) {
      newElement.elType = VSBEND;
      get2Token (is, &tk);            // -length
      newElement.length = tk.getReal();
      get2Token (is, &tk);            // e1
      newElement.e1 = tk.getReal();
      get2Token (is, &tk);            // e2
      newElement.e2 = tk.getReal();
      get2Token (is, &tk);            // angle
      newElement.angle = tk.getReal();
      get2Token (is, &tk);            // hgap
      newElement.hgap = tk.getReal();
      get2Token (is, &tk);            // fint
      newElement.fint = tk.getReal();
      get2Token (is, &tk);            // K0
      newElement.K0 = tk.getReal();
      get2Token (is, &tk);            // K0S
      newElement.K0S = tk.getReal();
      get2Token (is, &tk);            // K0
      newElement.K1 = tk.getReal();
      get2Token (is, &tk);            // K0S
      newElement.K1S = tk.getReal();
      get2Token (is, &tk);            // K0
      newElement.K2 = tk.getReal();
      get2Token (is, &tk);            // K0S
      newElement.K2S = tk.getReal();
      get2Token (is, &tk);            // K0
      newElement.K3 = tk.getReal();
      get2Token (is, &tk);            // K0S
      newElement.K3S = tk.getReal();
      get2Token (is, &tk);
      newElement.a=tk.getReal();
      get2Token (is, &tk);
      newElement.b=tk.getReal();
      get2Token (is, &tk);
      newElement.missX=tk.getReal();
      get2Token (is, &tk);
      newElement.missY=tk.getReal();
    } 
    else if (tk.isKey("VRBEND")) {
      newElement.elType = VRBEND;
      get2Token (is, &tk);            // -length
      newElement.length = tk.getReal();
      get2Token (is, &tk);            // e1
      newElement.e1 = tk.getReal();
      get2Token (is, &tk);            // e2
      newElement.e2 = tk.getReal();
      get2Token (is, &tk);            // angle
      newElement.angle = tk.getReal();
      get2Token (is, &tk);            // hgap
      newElement.hgap = tk.getReal();
      get2Token (is, &tk);            // fint
      newElement.fint = tk.getReal();
      get2Token (is, &tk);            // K0
      newElement.K0 = tk.getReal();
      get2Token (is, &tk);            // K0S
      newElement.K0S = tk.getReal();
      get2Token (is, &tk);            // K0
      newElement.K1 = tk.getReal();
      get2Token (is, &tk);            // K0S
      newElement.K1S = tk.getReal();
      get2Token (is, &tk);            // K0
      newElement.K2 = tk.getReal();
      get2Token (is, &tk);            // K0S
      newElement.K2S = tk.getReal();
      get2Token (is, &tk);            // K0
      newElement.K3 = tk.getReal();
      get2Token (is, &tk);            // K0S
      newElement.K3S = tk.getReal();
      get2Token (is, &tk);
      newElement.a=tk.getReal();
      get2Token (is, &tk);
      newElement.b=tk.getReal();
      get2Token (is, &tk);
      newElement.missX=tk.getReal();
      get2Token (is, &tk);
      newElement.missY=tk.getReal();
    } 
    else if  ( (tk.isKey("QUADRUPOLE")) || (tk.isKey("SEXTUPOLE")) ) {
      if (tk.isKey("QUADRUPOLE"))
	newElement.elType = QUADRUPOLE;
      else
	newElement.elType = SEXTUPOLE;
      get2Token (is, &tk);            // length
      newElement.length = tk.getReal();
      get2Token (is, &tk);            // K0
      newElement.K0 = tk.getReal();
      get2Token (is, &tk);            // K0S
      newElement.K0S = tk.getReal();
      get2Token (is, &tk);            // K1
      newElement.K1 = tk.getReal();
      get2Token (is, &tk);            // K1S
      newElement.K1S = tk.getReal();
      get2Token (is, &tk);            // K2
      newElement.K2 = tk.getReal();
      get2Token (is, &tk);            // K2S
      newElement.K2S = tk.getReal();
      get2Token (is, &tk);            // K3
      newElement.K3 = tk.getReal();
      get2Token (is, &tk);            // K3S
      newElement.K3S = tk.getReal();
      get2Token (is, &tk);
      newElement.a=tk.getReal();
      get2Token (is, &tk);
      newElement.b=tk.getReal();
      get2Token (is, &tk);
      newElement.missX=tk.getReal();
      get2Token (is, &tk);
      newElement.missY=tk.getReal();
    } 
    else if (tk.isKey("RECQUADRUPOLE")) {
      newElement.elType = RECQUADRUPOLE;
      get2Token (is, &tk);            // -length
      newElement.length = tk.getReal();
      get2Token (is, &tk);            // K0
      newElement.r1 = tk.getReal();
      get2Token (is, &tk);            // K0S
      newElement.r2 = tk.getReal();
      get2Token (is, &tk);              // a
      newElement.a=tk.getReal();
      get2Token (is, &tk);              // b
      newElement.b=tk.getReal();
      get2Token (is, &tk);
      newElement.missX=tk.getReal();
      get2Token (is, &tk);
      newElement.missY=tk.getReal();
    } else if (tk.isKey("RFCAVITY")) {
      newElement.elType = RFCAVITY;
      get2Token (is, &tk);              // L             
      newElement.length=tk.getReal();   
      get2Token (is, &tk);              // V             
      newElement.v=tk.getReal();   
      get2Token (is, &tk);              // PHI             
      newElement.phi=tk.getReal();   
      get2Token (is, &tk);              // FREQ
      newElement.freq=tk.getReal();        
      get2Token (is, &tk);              // a
      newElement.a=tk.getReal();
      get2Token (is, &tk);              // b
      newElement.b=tk.getReal();
      get2Token (is, &tk);
      newElement.missX=tk.getReal();
      get2Token (is, &tk);
      newElement.missY=tk.getReal();
    } 
    else if (tk.isKey("STATMARKER")) {
      newElement.elType = STATMARKER;
      newElement.length = 0.0;
    } 
    else if (tk.isKey("PARTMARKER")) {
      newElement.elType = PARTMARKER;
      newElement.length = 0.0;
    } 
    else if (tk.isKey("SCMAP")) {
      newElement.elType = SCMAP;
      get2Token (is, &tk);
      newElement.length = tk.getReal();
      get2Token (is, &tk);
      newElement.a=tk.getReal();
      get2Token (is, &tk);
      newElement.b=tk.getReal();
      get2Token (is, &tk);              // a
      newElement.missX=tk.getReal();
      get2Token (is, &tk);
      newElement.missY=tk.getReal();
    }
    else if  ( (tk.isKey("CONSTFOC")) ) {
      newElement.elType = CONSTFOC;
      get2Token (is, &tk);            // is circumference 
      newElement.length = tk.getReal();
      get2Token (is, &tk);            // K0 x 
      newElement.K0 = tk.getReal();
      get2Token (is, &tk);            // K1 y
      newElement.K1 = tk.getReal();
      get2Token (is, &tk);            // K2 z
      newElement.K2 = tk.getReal();
      get2Token (is, &tk);              // a
      newElement.a=tk.getReal();
      get2Token (is, &tk);              // b
      newElement.b=tk.getReal();
    } 
    newElement.isMisaligned = (newElement.missX != 0.0 || newElement.missY);
    return newElement; 
}
#endif
