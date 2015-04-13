/*

cavitytrack.cpp

Aim: Integrate particles through a box cavity using a second order Leap
Frog algorithm. No space charge and no external forces are considered.
     
Initial particle positions are read in from file.

Units:

x, y, z  : m        
time     : s
vx xy vz : m/s
B-field  : tesla
E-field  : volt/m

Usage:

*/

const unsigned int Dim = 3;

#include <ctime>
#include <iostream>
#include <fstream>
#include <mpi.h>
#include "CavitySolver.h"
#include "Ippl.h"
#ifdef USE_PBE
#include "pbe.h"
#endif
#include "Const.hh"
#include "ChargedParticles.hh"
typedef ChargedParticles<double,Dim>::Mesh_t Mesh_t;
typedef ChargedParticles<double,Dim>::Layout_t Layout_t;
typedef ChargedParticles<double,Dim>::Vector_t Vector_t;  
typedef ChargedParticles<double,Dim>::FieldLayout_t FieldLayout_t;
#include "cavitytrack.h"
#include "Distribution/Distribution.hh"
#include "Configure.h"

using namespace std; 

int main(int argc, char** argv) 
{
#ifdef __GNUC__
    std::set_terminate (__gnu_cxx::__verbose_terminate_handler);
#endif

    MPI_Init(&argc, &argv);    

    int turn = 0; 
    int i;
        
    // Configure the simulation
    // ************************
        
    double Emin;                // initial total energy of the reference particle in MeV
    double Emax;                // final energy, determines the numer of turns
    double PMCU;
    double PRMCU;
    double PTMCU;
    
    string eFieldFn;            // where the electric field data is stored
    string bFieldFn;            // where the B-Field data is stored
    string flatTopFn;           // where the electric field data for flat top cavities is stored
  
    double I0;                  // initial current in Ampere where Q is then: Q = I0/nuRf 
    double nuRf;                // RF frequency in Hz
    InterPolT interPol;         // CIC or NGP
    GreenT greenType;            
    
    BCT bc;                     // OOO or OOP
    DistT dist; 
    IntegratorT integrator;     // RK4 
    
    Vector_t semiaxi;           // sigma of initial distribution in meter
    Vector_t center;            // of axis definition in meter
    Vector_t vel;               // sigma of initial velocity spread 
    double thres[NEQ];
    double deflection[3];
    

    unsigned long totalP;       // number of macroparticles
    unsigned int nx,ny,nz;      // meshsizes
    unsigned int nadj;          // not used
    
    FsT fsType; 
    string distrInputFile;      // if we read in a distribution  
    string diagFile;            // not used
    string dataDir;             // directory where the outputdata gets stored
    string phaseSpaceData;      // basefilename of the phasespace fieles
    string statisticData;       // file name of statistical informations
    string collimatorFn("");    // file where the collimator are defined
    string probeFn("");         // file where the probes (Sonden) are defined
    
    double rad=0.0;
    double ang=0.0;
    double dZ=0.0;
    double initvelang=0.0;
    double phase=0.0;
    double dT=0.0;
    double TOL=0.0;
    bool accelstatus=false;

    double dphiStat;
    double dphiPart;
            
    bool writeRestart;
    bool readRestart;
    string restartFn;

    bool calcTune;
    int turns;

    const double pi = 4.0*atan(1.0);

    double r,phi,ddy[3],halfy[6]; 
    double rCenter,phiCenter;
        
    double phiAct = 0.0;

    DistrData distrData;

    double cavitylength= 1;                                 // the total length of drift-cavity-drift
    string meshfilename("");                           // the name of the cavity's meshfile
    double xtrans=0;
    double ytrans=0;
    double ztrans=0;
    double xrot=0;                                          // Cavity coordinate system parameters 
    double yrot=0;
    double zrot=0;
    double scaling=1;
    double peakvoltage = 0;
    
    CavitySolver *cavity_solver;

    ChargedParticles<double,Dim>* Bunch;

    double t = 0.0;

    ofstream of;    
    
    Ippl ippl(argc, argv);
    Inform msg("cavitytrack ");  
    Inform msg2all("cavitytrack ",INFORM_ALL_NODES); 
  
#ifdef USE_PBE
    pbe_init("PCYCLINT",10,0,NULL);
#endif
#ifdef USE_PBE
    pbe_start(1,"totalprog");
#endif
    int maxInbalance = 0;
    int cloneTurn  = -1;
    double cloneAngle = 0.0;
    unsigned int numOfNeighbors;
    
    bool res = cavityTrackConfigure(argc, argv, &Emin, &Emax, &eFieldFn, &bFieldFn,
				    &I0,&nuRf, &interPol,&greenType,&bc,&dist,&integrator,
				    semiaxi,center,vel,thres,deflection,&dT,&TOL,
				    &totalP,&nx,&ny,&nz,&nadj,&fsType, 
				    &distrInputFile,
				    &collimatorFn,&dataDir,&phaseSpaceData,&statisticData,&rad,&ang,&dZ,&initvelang,
				    &phase,&probeFn,&accelstatus,&writeRestart,&readRestart,&restartFn,
				    &calcTune,&turns,&dphiStat,&dphiPart,&maxInbalance,
				    &cloneTurn,&cloneAngle,&numOfNeighbors,distrData,&PMCU,&PRMCU,
				    &flatTopFn, 
				    &meshfilename,
				    &cavitylength,
				    &xtrans,&ytrans,&ztrans,&xrot,&yrot,&zrot,&scaling,&peakvoltage); 
    
     
    // ---------- Initial and derived values ----------------------
    
    double gamma= 1.0 + Emin/(m_p*1.0e3); 
    double beta =sqrt(1.0-(1.0/(gamma*gamma))); 
    double gammaFinal = 1.0 + Emax/(m_p*1.0e3);
    double Qtot = I0/nuRf; 
    double Mtot = Qtot*MASS/CHARGE;  

    //     << "P [MCU]              \t= " << PMCU << "    \t PR [MCU] = " << PRMCU << " \t PT [MCU] " << PTMCU << endl 
    if ((PMCU!=0.0) && (PRMCU!=0.0)) {
        PTMCU = sqrt((PMCU*PMCU)-(PRMCU*PRMCU));
    }
    else {
        PMCU = 1e3*sqrt(pow((Emin/MPMEV + 1.0),2) - 1.0);
        PTMCU = sqrt((PMCU*PMCU)-(PRMCU*PRMCU));
    }

    if (Ippl::myNode() == 0) {
        of.open(string(dataDir + "/diag.dat").c_str() ,ios::out);
        of << "# t, x, px, Excav, Bxcav, y, py, Eycav, Bycav, z, pz, Ezcav, Bzcav, gamma, ID ::: meshfile" << meshfilename << endl;
        of.close();
    }
       
    Index I(nx), J(ny), K(nz);
    NDIndex<Dim> domain1;
    domain1[0] = nx+1;
    domain1[1] = ny+1;
    domain1[2] = nz+1;  
  
    // For efficiency in the FFT's, we will use a parallel decomposition
    // which is serial in the first dimension.
    e_dim_tag decomp[Dim];
    decomp[0] = SERIAL; decomp[1] = SERIAL; decomp[2] = PARALLEL; 
    
    // create mesh and layout objects for this problem domain
    Mesh_t *mesh;
    FieldLayout_t *FL;
    mesh = new Mesh_t(domain1);
    FL = new FieldLayout_t(*mesh, decomp);
    Layout_t * PL = new Layout_t(*FL, *mesh); 
   
    // create an empty ChargedParticles object and set the BC's
    Bunch = new ChargedParticles<double,Dim>(PL,bc);    
    Bunch->setCoupling(1.0/(4.0*pi*EPS0));
    Bunch->setCurrent(I0);

    msg << "Bunch object created " << endl;
    
    Distribution *d = new Distribution(Bunch);
    dist = READFROMFILE;
    if (dist == READFROMFILE) 
        center = d->readInputDistribution(distrInputFile);
    else {
        ERRORMSG("Distribution not known" << endl);
        exit(1);
    }    
    Bunch->setN0();

    
    /*
      All particles have the same energy
    */
    for (unsigned long k=0; k<Bunch->getLocalNum(); k++) {
        Bunch->P[k](2) = beta*CLIGHT;
        Bunch->Ef[k] = Vector_t(0.0);
    }
    Bunch->updateBunches();

    IpplInfo::Comm->barrier();	
    
    msg << "/////////////////////////////////////////////////////////////////////////" << endl 
        << "Simulation parameters, using " << Ippl::getNodes() << " nodes(s)" 
        << endl;
    if (distrData.distrType==READFROMFILE) {
        msg << "Read particle distributon from file: " << distrInputFile << endl;
        dist = READFROMFILE;
    }
    else
        msg << "not known"  << endl;
    
    msg << "Integrator ";
    integrator = LEAPFROG;
    switch (integrator) {  
    case LEAPFROG: msg << "Leap Frog " << endl; break;
    default: msg << "not known I will use Leap Frog " << endl;
    }
    msg << "No Field solver selected " << endl;
    msg << "Diagnostics in " << string(dataDir + "/diag.dat") << endl;

    Vector_t p = Bunch->getMeanP() ;
    
    msg << "------------------------------------------------------------------" << endl
        << "E[MeV]= "  << Emin 
	<< " gamma= " << Bunch->getGamma() << " beta= " << beta << " p[mCU]= " << Bunch->mks2mcu(sqrt(dot(p,p))) << endl 
	<< "px[mCU]= " << Bunch->mks2mcu(p[0]) << " py[mCU]= " << Bunch->mks2mcu(p[1]) << " pz[mCU]= " << Bunch->mks2mcu(p[2]) 
        << endl
        << "RFfrequency[Hz]= " << nuRf 
        << "\t phase of cavity= " << phase
        << endl;
    msg << "dT= " << dT << endl;
    msg << "Cavity meshfile: " << meshfilename 
	<< " length of cavity [m] = " << cavitylength << endl;
    msg << "xtrans = " << xtrans 
	<< " ytrans = " << ytrans
	<< " ztrans = " << ztrans << endl;
    msg << "xrot = " << xrot 
	<< " yrot = " << yrot 
	<< " zrot = " << zrot << endl;
    msg << "scaling = " << scaling 
        << " phase = " << phase 
        << " peakvoltage [V] = " << peakvoltage << endl;
    msg << "/////////////////////////////////////////////////////////////////////////" << endl ;
        
    msg << "Initialize cavity_solver...";
    cavity_solver = new CavitySolver(meshfilename,cavitylength,xtrans,ytrans,ztrans,xrot,yrot,zrot,scaling,peakvoltage,phase);
    msg  << "....cavity_solver initialized." << endl << endl;
        
    msg << "Going to run the femaxx solver......" << endl;
    cavity_solver->runSolver();
    msg << "cavity_solver->runSolver successfully run. Going to start the time integration." << endl;
    
    
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // ********************* Time steping  ***********************
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    Bunch->h = t;    
        
    bool notFinish = true;        

    Vector_t cm(0.0);

    unsigned long eigenvector_n = 0; 
    double lambda = cavity_solver->getLambda(eigenvector_n);  // just needed for output here
    double omega = sqrt(lambda * pow(CLIGHT,2));   // just needed for output here 
    double length =  cavity_solver->getlength();  // just needed for output here 
    
    msg << "lambda= " << lambda << ", omega= " << omega << ", cl [m] =  " << length << endl;
    
    double curPos[3]  = {0.0,0.0,0.0};
    double Ecavity[3] = {0.0,0.0,0.0};
    double Bcavity[3] = {0.0,0.0,0.0};
    double E[3]       = {0.0,0.0,0.0};
    double B[3]       = {0.0,0.0,0.0};
              
    cm = Bunch->getMeanR(0);    

    for (unsigned int i=0; i<Dim; i++)
        curPos[i] = Bunch->R[3](i);
    
    double mult = 1.0;
    cavity_solver->getECavity(&curPos[0],&Ecavity[0],eigenvector_n,0);
    if (Ecavity[2]>0.0)
	mult=-1.0;

    beta = sqrt(dot(Bunch->P[0],Bunch->P[0]))/CLIGHT;
    gamma = 1.0/(sqrt(1.0-pow(beta,2)));
    
    msg << "t = " << t << " beta= " << beta << " Ek[Mev]= " <<  (gamma-1.0)*(m_p*1.0e3) << endl;    

    for (unsigned long k=0; k<Bunch->getLocalNum(); k++) {	
	for (unsigned int i=0; i<Dim; i++)
	    curPos[i] = Bunch->R[k](i);
	beta = sqrt(dot(Bunch->P[k],Bunch->P[k]))/CLIGHT;
	gamma = 1.0/(sqrt(1.0-pow(beta,2)));
	// here get the cavity fields 
	cavity_solver->getECavity(&curPos[0], &Ecavity[0],eigenvector_n, t);
	// B-cavity 
	cavity_solver->getBCavity(&curPos[0], &Bcavity[0],eigenvector_n, t);
	
	if (Bunch->ID[k]==3) {
	    of <<  Bunch->Ti[k] << ' ';
	    for (int jj = 0 ; jj < Dim ; ++jj) {
		of << setw(PWIC) << Bunch->R[k](jj) << ' ';
		of << setw(PWIC) << Bunch->P[k](jj) << ' ';
		of << setw(PWIC) << Ecavity[jj] << ' ';
		of << setw(PWIC) << Bcavity[jj] << ' ';
	    }
	    of << setw(PWIC) << gamma << ' ';
	    of << setw(PWIC) << Bunch->ID[k] << endl;
	}
    }
    
    while (notFinish) {
        for (unsigned long k=0; k<Bunch->getLocalNum(); k++) {
            double mult = CHARGE / MASS * dT / 2.0;
            for (unsigned int i=0; i<Dim; i++)
                curPos[i] = Bunch->R[k](i);
		
            // here get the cavity fields 
            cavity_solver->getECavity(&curPos[0], &Ecavity[0],eigenvector_n, Bunch->Ti[k]);
            // B-cavity 
            cavity_solver->getBCavity(&curPos[0], &Bcavity[0],eigenvector_n, Bunch->Ti[k]);
		
            // now do the time- integration
            Bunch->P[k](0) += (Ecavity[0] + (-1 * Bunch->P[k](2)) * Bcavity[1] +       Bunch->P[k](1)  * Bcavity[2]) * mult;
            Bunch->P[k](1) += (Ecavity[1] +       Bunch->P[k](2)  * Bcavity[0] + (-1 * Bunch->P[k](0)) * Bcavity[2]) * mult;
            Bunch->P[k](2) += (Ecavity[2] + (-1 * Bunch->P[k](1)) * Bcavity[0] +       Bunch->P[k](0)  * Bcavity[1]) * mult;
		
            Bunch->R[k](0) += Bunch->P[k](0) * dT;
            Bunch->R[k](1) += Bunch->P[k](1) * dT;
            Bunch->R[k](2) += Bunch->P[k](2) * dT;
		
            for (unsigned int i=0; i<Dim; i++)
                curPos[i] = Bunch-> R[k](i);
		
            // here get the cavity fields 
            cavity_solver->getECavity(&curPos[0], &Ecavity[0],eigenvector_n, Bunch->Ti[k]+dT);
            // B-cavity 
            cavity_solver->getBCavity(&curPos[0], &Bcavity[0],eigenvector_n, Bunch->Ti[k]+dT);
		
            Bunch->P[k](0) += (Ecavity[0] + (-1 * Bunch->P[k](2)) * Bcavity[1] +       Bunch->P[k](1)  * Bcavity[2]) * mult;
            Bunch->P[k](1) += (Ecavity[1] +       Bunch->P[k](2)  * Bcavity[0] + (-1 * Bunch->P[k](0)) * Bcavity[2]) * mult;
            Bunch->P[k](2) += (Ecavity[2] + (-1 * Bunch->P[k](1)) * Bcavity[0] +       Bunch->P[k](0)  * Bcavity[1]) * mult;
	  
            Bunch->h[k]   = dT;
            Bunch->Ti[k] += dT;
	    
            beta = sqrt(dot(Bunch->P[k],Bunch->P[k]))/CLIGHT;
            gamma = 1.0/(sqrt(1.0-pow(beta,2)));
            if (Bunch->ID[k]==3) {
                of.open(string(dataDir + "/diag.dat").c_str() ,ios::app);
                of <<  Bunch->Ti[k] << ' ';
                for (int jj = 0 ; jj < Dim ; ++jj) {
                    of << setw(PWIC) << Bunch->R[k](jj) << ' ';
                    of << setw(PWIC) << Bunch->P[k](jj) << ' ';
                    of << setw(PWIC) << Ecavity[jj] << ' ';
                    of << setw(PWIC) << Bcavity[jj] << ' ';
                }
                of << setw(PWIC) << gamma << ' ';
                of << setw(PWIC) << Bunch->ID[k] << endl;
                of.close();
            }
	}
        Bunch->updateBunches();
        t+=dT; 
	cm = Bunch->getMeanR(0);
        notFinish = (abs(cm[2]) < cavitylength);
        IpplInfo::Comm->barrier();	
	msg << "t = " << t << " beta= " << beta << " Ek[Mev]= " <<  (gamma-1.0)*(m_p*1.0e3) << endl;
    }
#ifdef USE_PBE
    pbe_stop(1);
    pbe_dump(); 
    pbe_finalize(1);
#endif

    // FIXME: Some deletes are missing.
    delete cavity_solver;

    MPI_Finalize();
    return 0;
}
