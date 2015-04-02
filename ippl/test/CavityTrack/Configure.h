#ifndef CONFIGURE_H
#define CONFIGURE_H

class DistrData {
 public:
    DistT distrType;
    unsigned long N;
    
    double sigx;
    double sigpx;
    double corxpx;
    double emitx;
    double mx;

    double sigy;
    double sigpy;
    double corypy;
    double emity;
    double my;

    double sigz;
    double sigpz;
    double corzpz;
    double emitz;
    double mz;
};

bool cavityTrackConfigure(int argc, char *argv[], 
			  double *Emin, 
			  double *Emax,
			  string *eFieldFn,
			  string *bFieldFn,
			  double *I0, double *nuRf, 
			  InterPolT *interPol, GreenT *greenType,
			  BCT *bc, DistT *dist, IntegratorT *integr,
			  Vector_t &semiaxi,
			  Vector_t &center,
			  Vector_t &vel,
			  double threshold[],
			  double deflection[],
			  double *dT, double *Tol,
			  unsigned long *totalP,
			  unsigned int *nx, unsigned int *ny, unsigned int *nz,
			  unsigned int *nadj,
			  FsT *fsType, string *distInpFn,
			  string *colDefFile, string *datadir,
			  string *phaseSpaceFn, string *statisticFn,
			  double *initialRadius,
			  double *initialPhase,
			  double *deltaZ,
			  double *initvelang,
			  double *phaseofcavity,
			  string *probeFn,
			  bool *accelstatus,
			  bool *writeRestart,
			  bool *readRestart,
			  string *restartFn,
			  bool *tune,
			  int *turns,
			  double *dphiStat,
			  double *dphiPart,
			  int *maxInbalance,
			  int *saveTurn, 
			  double *saveAngle, 
			  unsigned int *numOfNeighbors,
			  DistrData &distrData,
			  double *PCU, double *PRCU,
			  string *flatTopFn,
			  string *meshfilename, 
			  double *cavitylength,
			  double *xtrans,
			  double *ytrans,
			  double *ztrans,
			  double *xrot,
			  double *yrot,
			  double *zrot,
			  double *scaling,
			  double *peakvoltage)
{    
    *writeRestart = *readRestart = false;
    *tune = false;
    *turns = 51;
    *integr = RK78;
    *interPol = CIC;
    *nadj=0;
    *fsType=None;
    *bc=OOO;
    *dist=ELLIPSOIDALUniform;        

    *dphiStat = 1.0;     // output in 1.0 degree steps statistical data        
    *dphiPart = 10.0;    // output in 10.0 degree steps phase space data        
    *maxInbalance = 20;  // we allow maximal 20% unbalance in the particle distribution among the processors

    *saveTurn = -1;
    *saveAngle = 0.0;
    *numOfNeighbors = 0;

    distrData.distrType=NOTKNOWN;

    *PCU = 0.0;
    *PRCU = 0.0;

    *flatTopFn = string("");
    
    for (int i=1; i < argc; ++i) {
        string s(argv[i]);
        if (s == "-FFTSolver")
            *fsType = FFTSolver;
        if (s == "-grid") {
            *nx = atoi(argv[++i]);
            *ny = atoi(argv[++i]);
            *nz = atoi(argv[++i]);
        } else if (s == "-neighbors") {
            *saveTurn = atoi(argv[++i]);
            *saveAngle = atof(argv[++i]);
            *numOfNeighbors = atoi(argv[++i]);
        } else if (s == "-inbalance" ) {
            *maxInbalance = atoi(argv[++i]);
        } else if (s == "-initR" ) {
            *initialRadius = atof(argv[++i]);
        } else if (s == "-dphiStat" ) {
            *dphiStat = atof(argv[++i]);
        } else if (s == "-dphiPart" ) {
            *dphiPart = atof(argv[++i]);
        } else if (s == "-initPhi" ) {
            *initialPhase = atof(argv[++i]);
        } else if (s == "-initdZ" ) {
            *deltaZ = atof(argv[++i]); 
        } else if (s == "-initialangleofvelocity" ) {
            *initvelang = atof(argv[++i]);
        } else if (s == "-phaseofcavity" ) {
            *phaseofcavity = atof(argv[++i]);
        } else if (s == "-accelerationstatus" ) {
	  *accelstatus = (atoi(argv[++i])==0) ? false: true ;
        } else if (s == "-Threshold") {
            threshold[0] = atof(argv[++i]);
            threshold[1] = atof(argv[++i]);
            threshold[2] = atof(argv[++i]);
            threshold[3] = atof(argv[++i]);
            threshold[4] = atof(argv[++i]);
            threshold[5] = atof(argv[++i]);
        }  else if (s == "-deflection") {
            deflection[0] = atof(argv[++i]);
            deflection[1] = atof(argv[++i]);
            deflection[2] = atof(argv[++i]);
        } else if (s == "-dT" ) {
            *dT = atof(argv[++i]);
        } else if (s == "-Tol" ) {
            *Tol = atof(argv[++i]);
        } else if (s == "-Emin" ) {
            *Emin = atof(argv[++i]);
        } else if (s == "-Emax" ) {
            *Emax = atof(argv[++i]);
        } else if (s == "-EfieldFn") {
            *eFieldFn = string(argv[++i]);
        } else if (s == "-probeFn") {
            *probeFn = string(argv[++i]);
        } else if (s == "-BfieldFn") {
            *bFieldFn = string(argv[++i]);
        } else if (s == "-FlatTopfieldFn") {
            *flatTopFn   = string(argv[++i]);
        }
        else if (s == "-I0") {
            *I0 = atof(argv[++i]);
        } else if (s == "-nuRf") {
            *nuRf = atof(argv[++i]);
        }  else if (s == "-semiaxi") {
            semiaxi[0] = atof(argv[++i]);
            semiaxi[1] = atof(argv[++i]);
            semiaxi[2] = atof(argv[++i]);
        } else if (s == "-center") {
            center[0] = atof(argv[++i]);
            center[1] = atof(argv[++i]);
            center[2] = atof(argv[++i]);
        } else if (s == "-vel") {
            vel[0] = atof(argv[++i]);
            vel[1] = atof(argv[++i]);
            vel[2] = atof(argv[++i]);
        } else if (s == "-N") {
            *totalP = static_cast<unsigned int>(atof(argv[++i]));
        } else if (s == "-distrInpFn") {
            *distInpFn = string(argv[++i]);
	    distrData.distrType = READFROMFILE;
        } else if (s == "-collimatorFn") {
            *colDefFile = string(argv[++i]);
        } else if (s == "-dataDir") {
            *datadir = string(argv[++i]);
        } else if (s == "-phaseSpaceData") {
            *phaseSpaceFn = string(argv[++i]);
        } else if (s == "-statisticData") {
            *statisticFn = string(argv[++i]);
        } else if (s == "-writeRestart") {
            *writeRestart = true;
            *restartFn = string("restart.dat");
        } else if (s == "-readRestart") {
            *readRestart = true;
            *restartFn = string(argv[++i]);
        } else if (s == "-Integrator") {
	  string intStr = string(argv[++i]);
	  if (intStr == string("RK45"))
	    *integr = RK45;
          else if (intStr == string("RK4"))
	    *integr = RK4;
          else if (intStr == string("RK4A"))
	    *integr = RK4A;
	  else if (intStr == string("RK78"))
	    *integr = RK78;
	  else if (intStr == string("LeapFrog"))
	    *integr = LEAPFROG;
	  else if (intStr == string("VelVerlet"))
	    *integr = VELVERLET;
        } else if (s == "-tune") {
            *tune = true;
        } else if (s == "-turns") {
            *turns = atoi(argv[++i]);
        } else if (s == "-distr") {
	  if ("GAUSS" == string(argv[++i]))
	      distrData.distrType = GAUSS;
	  if ("sigx" == string(argv[++i]))
	    distrData.sigx =  atof(argv[++i]);
	  if ("sigy" == string(argv[++i]))
	    distrData.sigy =  atof(argv[++i]);
	  if ("sigz" == string(argv[++i]))
	    distrData.sigz =  atof(argv[++i]);
	  if ("sigpx" == string(argv[++i]))
	    distrData.sigpx =  atof(argv[++i]);
	  if ("sigpy" == string(argv[++i]))
	    distrData.sigpy =  atof(argv[++i]);
	  if ("sigpz" == string(argv[++i]))
	    distrData.sigpz =  atof(argv[++i]); 
	  if ("xpx" == string(argv[++i]))
	    distrData.corxpx =  atof(argv[++i]);
	  if ("ypy" == string(argv[++i]))
	    distrData.corypy =  atof(argv[++i]);
	  if ("zpz" == string(argv[++i]))
	      distrData.corzpz =  atof(argv[++i]);
        } else if (s == "-PCU") {
            *PCU = atof(argv[++i]);
        } else if (s == "-PRCU") {
            *PRCU = atof(argv[++i]);
	} else if (s == "-cavitylength") {       //this is new for the length of drift-cavity-drift
	    *cavitylength = atof(argv[++i]);
	} else if (s == "-cavitymeshfile") {     // and this is new for cavity meshfile
	    *meshfilename = string(argv[++i]);    
	} else if (s == "-xtrans") {             // next 7 ifs are for cavity coordinate system parameters
	    *xtrans = atof(argv[++i]);     
	} else if (s == "-ytrans") {             
	    *ytrans = atof(argv[++i]);
	} else if (s == "-ztrans") {             
	    *ztrans = atof(argv[++i]);
	} else if (s == "-xrot") {             
	    *xrot = atof(argv[++i]);
	} else if (s == "-yrot") {             
	    *yrot = atof(argv[++i]);
	} else if (s == "-zrot") {             
	    *zrot = atof(argv[++i]);
	} else if (s == "-scaling") {             
	    *scaling = atof(argv[++i]);
	} else if (s == "-peakvoltage") {             
	    *peakvoltage = atof(argv[++i]);
	}
    }
    return true; 
}
#endif
