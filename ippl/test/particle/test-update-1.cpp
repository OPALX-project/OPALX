// -*- C++ -*-
/**************************************************************************************************************************************
 *
 * The IPPL Framework
 *
 * This program was prepared by PSI.
 * All rights in the program are reserved by PSI.
 * Neither PSI nor the author(s)
 * Neither PSI nor the author(s)
 * makes any warranty, express or implied, or assumes any liability or
 * responsibility for the use of this software
 *

 This test program sets up a simple sine-wave electric field in 3D,
   creates a population of particles with random q/m values (charge-to-mass
   ratio) and velocities, and then tracks their motions in the static
   electric field using nearest-grid-point interpolation.

Usage:                                                                         Serial dimension
                                               Time steps                      |
                                        N-Part |                               |     Movement test [PUSH|DRIFT|EDGE|NONE]
                            grid size   |      |                               |      |
                            |           |      |                               |      | 
 mpirun -np 2 test-update-1 128 128 128 10000 10 NGP [OOO|PPP] [GUARDCELLS|no] 2 5.0 PUSH --commlib mpi --info 0
                                                 |                               |
                                                 | [NGP|CIC]                     | Defines the maximum particle inbalance allowed (%)

Example:
mpirun -np 32 test-update-1 32 32 32 300 400 NGP OOO GUARDCELLS 3 5.0 --commlib mpi --info 0 | tee test-update-1.out

*************************************************************************************************************************************/

#include "Ippl.h"
#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include <set>
#include <limits>
#include <cmath>

#include "mpi.h"

#define noVTKOUT 1
#define GUARDCELL 1


// dimension of our positions
const unsigned Dim = 3;

// some typedefs
typedef ParticleSpatialLayout<double,Dim>::SingleParticlePos_t Vector_t;
typedef ParticleSpatialLayout<double,Dim> playout_t;
typedef UniformCartesian<Dim,double> Mesh_t;
typedef Cell                                       Center_t;
typedef CenteredFieldLayout<Dim, Mesh_t, Center_t> FieldLayout_t;
typedef Field<double, Dim, Mesh_t, Center_t>       Field_t;
typedef Field<Vector_t, Dim, Mesh_t, Center_t>     VField_t;
typedef IntCIC IntrplCIC_t;
typedef IntNGP IntrplNGP_t;

typedef Field<dcomplex, Dim, Mesh_t, Center_t>     CxField_t;
typedef FFT<RCTransform, Dim, double>              FFT_t;

enum BC_t {OOO,OOP,PPP};
enum InterPol_t {NGP,CIC};
enum Movement_t {NONE, PUSH, DRIFT, EDGE, WALK};

const double pi = acos(-1.0);
const double qmmax = 1.0;       // maximum value for particle q/m
const double dt = 1.0;          // size of timestep


void dumpVTK(Field<Vektor<double,3>,3> &EFD, NDIndex<3> lDom, int nx, int ny, int nz, int iteration,
             double dx, double dy, double dz) {

    ofstream vtkout;
    vtkout.precision(10);
    vtkout.setf(ios::scientific, ios::floatfield);

    std::stringstream fname;
    fname << "data/ef_";
    fname << setw(4) << setfill('0') << iteration;
    fname << ".vtk";

    //SERIAL at the moment
    //if (Ippl::myNode() == 0) {

    // open a new data file for this iteration
    // and start with header
    vtkout.open(fname.str().c_str(), ios::out);
    vtkout << "# vtk DataFile Version 2.0" << endl;
    vtkout << "pic3d" << endl;
    vtkout << "ASCII" << endl;
    vtkout << "DATASET STRUCTURED_POINTS" << endl;
    vtkout << "DIMENSIONS " << nx << " " << ny << " " << nz << endl;
    vtkout << "ORIGIN 0 0 0" << endl;
    vtkout << "SPACING " << dx << " " << dy << " " << dz << endl;
    vtkout << "POINT_DATA " << nx*ny*nz << endl;

    vtkout << "VECTORS E-Field float" << endl;
    for(int z=lDom[2].first(); z<lDom[2].last(); z++) {
        for(int y=lDom[1].first(); y<lDom[1].last(); y++) {
            for(int x=lDom[0].first(); x<lDom[0].last(); x++) {
                Vektor<double, 3> tmp = EFD[x][y][z].get();
                vtkout << tmp(0) << "\t"
                       << tmp(1) << "\t"
                       << tmp(2) << endl;
            }
        }
    }

    // close the output file for this iteration:
    vtkout.close();
}


void dumpVTK(Field<double,3> &EFD, NDIndex<3> lDom, int nx, int ny, int nz, int iteration,
             double dx, double dy, double dz) {

    ofstream vtkout;
    vtkout.precision(10);
    vtkout.setf(ios::scientific, ios::floatfield);

    std::stringstream fname;
    fname << "data/scalar_";
    fname << setw(4) << setfill('0') << iteration;
    fname << ".vtk";

    //SERIAL at the moment
    //if (Ippl::myNode() == 0) {

    // open a new data file for this iteration
    // and start with header
    vtkout.open(fname.str().c_str(), ios::out);
    vtkout << "# vtk DataFile Version 2.0" << endl;
    vtkout << "toyfdtd" << endl;
    vtkout << "ASCII" << endl;
    vtkout << "DATASET STRUCTURED_POINTS" << endl;
    vtkout << "DIMENSIONS " << nx << " " << ny << " " << nz << endl;
    vtkout << "ORIGIN 0 0 0" << endl;
    vtkout << "SPACING " << dx << " " << dy << " " << dz << endl;
    vtkout << "POINT_DATA " << nx*ny*nz << endl;

    vtkout << "SCALARS E-Field float" << endl;
    vtkout << "LOOKUP_TABLE default" << endl;
    for(int z=lDom[2].first(); z<=lDom[2].last(); z++) {
        for(int y=lDom[1].first(); y<=lDom[1].last(); y++) {
            for(int x=lDom[0].first(); x<=lDom[0].last(); x++) {
                vtkout << EFD[x][y][z].get() << endl;
            }
        }
    }

    // close the output file for this iteration:
    vtkout.close();
}




template<class PL>
class ChargedParticles : public IpplParticleBase<PL> {
public:
    ParticleAttrib<double>     qm; // charge-to-mass ratio
    typename PL::ParticlePos_t E;  // electric field at particle position

    ChargedParticles(PL* pl, BC_t bc, InterPol_t interpol, e_dim_tag decomp[Dim], bool gCells) :
        IpplParticleBase<PL>(pl),
        bco_m(bc),
        interpol_m(interpol),
        fieldNotInitialized_m(true),
        withGuardCells_m(gCells)
    {
        // register the particle attributes
        this->addAttribute(qm);
        this->addAttribute(E);
        setupBCs();
        for(int i=0; i<Dim; i++)
            decomp_m[i]=decomp[i];

	partPerNode_m = new double[Ippl::getNodes()];
	globalPartPerNode_m = new double[Ippl::getNodes()];

	gradienTimer_m = IpplTimings::getTimer("gradienTimer");
	divergeTimer_m = IpplTimings::getTimer("divergeTimer");
	loadbalTimer_m = IpplTimings::getTimer("loadbalTimer");
	parpushTimer_m = IpplTimings::getTimer("parpushTimer");
	myupdatTimer_m = IpplTimings::getTimer("myupdatTimer");
	ngpgathTimer_m = IpplTimings::getTimer("ngpgathTimer");
	cicgathTimer_m = IpplTimings::getTimer("cicgathTimer");
	ngpscatTimer_m = IpplTimings::getTimer("ngpscatTimer");
	cicscatTimer_m = IpplTimings::getTimer("cicscatTimer");
    pardrifTimer_m = IpplTimings::getTimer("pardrifTimer");
    paredgeTimer_m = IpplTimings::getTimer("paredgeTimer");
	parwalkTimer_m = IpplTimings::getTimer("parwalkTimer");

	/*
	  Create the fft object
	*/
	// NDIndex<Dim> ndiStandard0h = getFieldLayout().getDomain();
	// ndiStandard0h[0] = Index(nr_m[0]/2+1);
	// fft_m = new FFT_t(getFieldLayout().getDomain(), ndiStandard0h[0], true);
    } 
    
    /*
      In case we have OOP or PPP boundary conditions
      we must define the domain, i.e can not be deduced from the
      particles as in the OOO case.
    */

    ChargedParticles(PL* pl, BC_t bc, InterPol_t interpol, Vector_t hr, Vector_t rmin, Vector_t rmax, e_dim_tag decomp[Dim], bool gCells) :
        IpplParticleBase<PL>(pl),
        bco_m(bc),
        interpol_m(interpol),
        hr_m(hr),
        rmin_m(rmin),
        rmax_m(rmax),
        fieldNotInitialized_m(true),
        withGuardCells_m(gCells)
    {
        // register the particle attributes
        this->addAttribute(qm);
	this->addAttribute(E);
     
        setupBCs();
        for(int i=0; i<Dim; i++)
            decomp_m[i]=decomp[i];

	partPerNode_m = new double[Ippl::getNodes()];
	globalPartPerNode_m = new double[Ippl::getNodes()];

	gradienTimer_m = IpplTimings::getTimer("gradienTimer");
	divergeTimer_m = IpplTimings::getTimer("divergeTimer");
	loadbalTimer_m = IpplTimings::getTimer("loadbalTimer");
	parpushTimer_m = IpplTimings::getTimer("parpushTimer");
	myupdatTimer_m = IpplTimings::getTimer("myupdatTimer");
	ngpgathTimer_m = IpplTimings::getTimer("ngpgathTimer");
	cicgathTimer_m = IpplTimings::getTimer("cicgathTimer");
	ngpscatTimer_m = IpplTimings::getTimer("ngpscatTimer");
	cicscatTimer_m = IpplTimings::getTimer("cicscatTimer");
    pardrifTimer_m = IpplTimings::getTimer("pardrifTimer");
    paredgeTimer_m = IpplTimings::getTimer("paredgeTimer");
	parwalkTimer_m = IpplTimings::getTimer("parwalkTimer");

	/*
	  Create the fft object

	  create half-size domain for RC transform along zeroth axis                                                                                                 
	*/
	NDIndex<Dim> ndiStandard0h = getFieldLayout().getDomain();
	ndiStandard0h[0] = Index(nr_m[0]/2+1);
	//	fft_m = new FFT_t(getFieldLayout().getDomain(), ndiStandard0h[0], true);
    }
    
    void setupBCs() {
      if (bco_m == OOO)
	setBCAllOpen();
      if (bco_m == PPP)
	setBCAllPeriodic();
      if (bco_m == OOP)
	setBCOOP();
    }

    inline const Mesh_t& getMesh() const { return this->getLayout().getLayout().getMesh(); }

    inline Mesh_t& getMesh() { return this->getLayout().getLayout().getMesh(); }

    inline const FieldLayout_t& getFieldLayout() const {
        return dynamic_cast<FieldLayout_t&>( this->getLayout().getLayout().getFieldLayout());
    }

    inline FieldLayout_t& getFieldLayout() {
        return dynamic_cast<FieldLayout_t&>(this->getLayout().getLayout().getFieldLayout());
    }

    void gather(int iteration) {

        if (interpol_m==CIC)
            gatherCIC();
        else
            gatherNGP();

	#ifdef VTKOUT
        NDIndex<Dim> lDom = getFieldLayout().getLocalNDIndex();
        dumpVTK(scalF_m,lDom,nr_m[0],nr_m[1],nr_m[2],iteration,hr_m[0],hr_m[1],hr_m[2]);
	#endif
    }

    double calcGridCharge() {
       Inform m ("calcGridCharge() ");
      
       /*
	Gives the same results thatn sum(scalF_m) !
       */

	NDIndex<Dim> elem;
	double Q = 0.0;
        const int me = Ippl::myNode();

        scalF_m.fillGuardCells();
	
        FieldLayout<Dim> & FL = scalF_m.getLayout();
	NDIndex<Dim> gDom = FL.getDomain();

	vector<NDIndex<Dim> > lDoms;

	lDoms.resize(Ippl::getNodes());
	lDoms[me] = FL.getLocalNDIndex();

	for (FieldLayout<Dim>::iterator_dv rdv = FL.begin_rdv(); rdv != FL.end_rdv(); ++ rdv) 
	  lDoms[((*rdv).second)->getNode()] = ((*rdv).second)->getDomain();
	
	for (int i = lDoms[me][0].first(); i <= lDoms[me][0].last(); ++ i)  {
	  elem[0] = Index(i,i);
	  for (int j = lDoms[me][1].first(); j <= lDoms[me][1].last(); ++ j)  {
	    elem[1] = Index(j,j);
	    for (int k = lDoms[me][2].first(); k <= lDoms[me][2].last(); ++ k)  {
	      elem[2] = Index(k,k);
	      Q += scalF_m.localElement(elem); 
	    }
	  }
	}
        reduce(Q,Q,OpAddAssign());
	return Q;
    }


    void scatter() {
        Inform m ("scatter: " );
	m.precision(10);
	m.setf(ios::fixed,ios::floatfield);

	scalF_m = 0.0;
	
        if (interpol_m==CIC)
	  scatterCIC();
        else
	  scatterNGP();
	
	double Q = calcGridCharge();
	
	m << "Qgrid= " << Q << " sum(qe)= " << sum(qm) << " deltaQ= " << (sum(qm)-Q) << endl;
    }
    
    bool hasZeroParticlesOnNode() {
      /*
	getLoadInbalance must be called before!
       */
      	bool low = false;
	for (int i=0; i<Ippl::getNodes(); i++) {
	  low = (globalPartPerNode_m[i] == 0);
	  if (low)
	    break;
	}
	return low;
    }

    pair<double,double> getMaxMinLoadInbalance() {
      /*
	getLoadInbalance must be called before!
       */
      size_t totpart  = this->getTotalNum();
      double idealdis = totpart/Ippl::getNodes();
      double maxnpart = (double)std::numeric_limits<size_t>::min();
      double minnpart = (double)std::numeric_limits<size_t>::max();
      for (int i=0; i<Ippl::getNodes(); i++) {
	if (maxnpart < globalPartPerNode_m[i])
	  maxnpart = globalPartPerNode_m[i];
	if (minnpart > globalPartPerNode_m[i])
	  minnpart = globalPartPerNode_m[i];
      }
      // now calculate the difference to the ideal ditribution in %
      maxnpart = abs(maxnpart-idealdis)/totpart*100.0;
      minnpart = abs(minnpart-idealdis)/totpart*100.0;
      return pair<double,double>(maxnpart,minnpart);
    }

    void getLoadInbalace() {
      for (int i=0; i<Ippl::getNodes(); i++)
	partPerNode_m[i] = globalPartPerNode_m[i] = 0.0;
      partPerNode_m[Ippl::myNode()] = this->getLocalNum();
      reduce(partPerNode_m,partPerNode_m+Ippl::getNodes(),globalPartPerNode_m,OpAddAssign());
    }
 
    void checkLoadBalance(double threshold) {
        Inform m("checkLoadBalance ");
	IpplTimings::startTimer(this->loadbalTimer_m);
        getLoadInbalace();
	pair<double,double> maxmin = getMaxMinLoadInbalance();
	if ((std::max(maxmin.first,maxmin.second) > threshold) || hasZeroParticlesOnNode()) {
	  m << "Load imbalance " << std::max(maxmin.first,maxmin.second); 
	  BinaryRepartition(*this);
	  getLoadInbalace();
	  maxmin = getMaxMinLoadInbalance();
	  m << " (%) AFTER REPARTITIONING: load imbalance " << std::max(maxmin.first,maxmin.second) << " (%) " << endl;
	}
      	IpplTimings::stopTimer(this->loadbalTimer_m);
    }

    void initFields() {

        NDIndex<Dim> domain = getFieldLayout().getDomain();

        for(int i=0; i<Dim; i++)
            nr_m[i] = domain[i].length();

        int nx = nr_m[0];
        int ny = nr_m[1];
        int nz = nr_m[2];

        double phi0 = 0.1*nx;

        Index I(nx), J(ny), K(nz);

        assign(vectF_m[I][J][K](0), -2.0*pi*phi0/nx * cos(2.0*pi*(I+0.5)/nx) * cos(4.0*pi*(J+0.5)/ny) * cos(pi*(K+0.5)/nz));

        assign(vectF_m[I][J][K](1),  4.0*pi*phi0/ny * sin(2.0*pi*(I+0.5)/nx) * sin(4.0*pi*(J+0.5)/ny));

        assign(vectF_m[I][J][K](2),  4.0*pi*phi0/ny * sin(2.0*pi*(I+0.5)/nx) * sin(4.0*pi*(J+0.5)/ny));

        assign(scalF_m[I][J][K],
               vectF_m[I][J][K](0) * vectF_m[I][J][K](0) +
               vectF_m[I][J][K](1) * vectF_m[I][J][K](1) +
               vectF_m[I][J][K](2) * vectF_m[I][J][K](2));
    }

    Vector_t getRMin() { return rmin_m;}
    Vector_t getRMax() { return rmax_m;}
    Vector_t getHr() { return hr_m;}

    void setRMin(Vector_t x) { rmin_m = x; }
    void setHr(Vector_t x) { hr_m = x; }

    void savePhaseSpace(string fn, int idx) {

        int tag = Ippl::Comm->next_tag(IPPL_APP_TAG4, IPPL_APP_CYCLE);
        vector<double> tmp;
        tmp.clear();

        // every node checks if he has to dump particles
        for(unsigned i=0; i<this->getLocalNum(); i++) {
            tmp.push_back(this->ID[i]);
            tmp.push_back(this->R[i](0));
            tmp.push_back(this->R[i](1));
            tmp.push_back(this->R[i](2));
            tmp.push_back(this->P[i](0));
            tmp.push_back(this->P[i](1));
            tmp.push_back(this->P[i](2));
        }

        if(Ippl::myNode() == 0) {
            ofstream ostr;
            string Fn;
            char numbuf[6];
            sprintf(numbuf, "%05d", idx);
            Fn = fn  + string(numbuf) + ".dat";
            ostr.open(Fn.c_str(), ios::out );
            ostr.precision(15);
            ostr.setf(ios::scientific, ios::floatfield);

            ostr << " x, px, y, py t, pt, id, node" << endl;

            unsigned int dataBlocks=0;
            double x,y,z,px,py,pz,id;
            unsigned  vn;

            for (unsigned i=0; i < tmp.size(); i+=7)
                ostr << tmp[i+1] << " " << tmp[i+4] << " " << tmp[i+2]  << " " \
                     << tmp[i+5] << " " << tmp[i+3] << " " << tmp[i+6]  << " " \
                     << tmp[i]   << " " << Ippl::myNode() << endl;

            int notReceived =  Ippl::getNodes() - 1;
            while (notReceived > 0) {
                int node = COMM_ANY_NODE;
                Message* rmsg =  Ippl::Comm->receive_block(node, tag);
                if (rmsg == 0)
                    ERRORMSG("Could not receive from client nodes in main." << endl);
                notReceived--;
                rmsg->get(&dataBlocks);
                rmsg->get(&vn);
                for (unsigned i=0; i < dataBlocks; i+=7) {
                    rmsg->get(&id);
                    rmsg->get(&x);
                    rmsg->get(&y);
                    rmsg->get(&z);
                    rmsg->get(&px);
                    rmsg->get(&py);
                    rmsg->get(&pz);
                    ostr << x  << "\t " << px  << "\t " << y  << "\t " \
                         << py << "\t " << z << "\t " << pz << "\t "   \
                         << id   << "\t " << vn << endl;
                }
                delete rmsg;
            }
            ostr.close();
        }
        else {
            unsigned dataBlocks = 0;
            Message* smsg = new Message();
            dataBlocks = tmp.size();
            smsg->put(dataBlocks);
            smsg->put(Ippl::myNode());
            for (unsigned i=0; i < tmp.size(); i++)
                smsg->put(tmp[i]);
            bool res = Ippl::Comm->send(smsg, 0, tag);
            if (! res)
                ERRORMSG("Ippl::Comm->send(smsg, 0, tag) failed " << endl;);
        }
    }

    inline void gradTest() {
      vectF_m = Vector_t(0.0);
      IpplTimings::startTimer(this->gradienTimer_m);
      vectF_m = -Grad(scalF_m,vectF_m);
      IpplTimings::stopTimer(this->gradienTimer_m);
    }

    inline void divTest() {
      IpplTimings::startTimer(this->divergeTimer_m);
      scalF_m = Div(vectF_m,scalF_m);
      IpplTimings::stopTimer(this->divergeTimer_m);
    }
  

    void myUpdate() {
        IpplTimings::startTimer(this->myupdatTimer_m);
	double hz;
	double zmin;
	double zmax;

        if (bco_m != PPP) {

	  // in case of OOP, we have already set the size of the parallel dimension.	  
	  if (bco_m == OOP) {
	    hz   = hr_m[2];
	    zmin = rmin_m[2];
	    zmax = rmax_m[2];
	  }
	  
	  bounds(this->R, rmin_m, rmax_m);
	
	  double stretch = 0.000;
	  Vector_t diff = rmax_m - rmin_m;
	  rmax_m += stretch*diff;
	  rmin_m -= stretch*diff;
	
	  NDIndex<Dim> domain = this->getFieldLayout().getDomain();

	  for(int i=0; i<Dim; i++)
	    nr_m[i] = domain[i].length();

	  for(int i=0; i<Dim; i++)
	    hr_m[i] = (rmax_m[i] - rmin_m[i]) / (nr_m[i] - 1.0);

	  if (bco_m == OOP) {  
	    rmin_m[2] = zmin;
	    rmax_m[2] = zmax;
	    hr_m[2] = hz;
	  }
	  getMesh().set_meshSpacing(&(hr_m[0]));
	  getMesh().set_origin(rmin_m);

	  if(withGuardCells_m) {
	    vectF_m.initialize(getMesh(), getFieldLayout(), GuardCellSizes<Dim>(1), vbc_m);
	    scalF_m.initialize(getMesh(), getFieldLayout(), GuardCellSizes<Dim>(1), bc_m);
	  }
	  else {
	    vectF_m.initialize(getMesh(), getFieldLayout(), vbc_m);
	    scalF_m.initialize(getMesh(), getFieldLayout(), bc_m);
	  }
        }
        else { // case PPP
	  if(fieldNotInitialized_m) {
	    fieldNotInitialized_m=false;
	    getMesh().set_meshSpacing(&(hr_m[0]));
	    getMesh().set_origin(rmin_m);
	    if(withGuardCells_m) {
	      vectF_m.initialize(getMesh(), getFieldLayout(), GuardCellSizes<Dim>(1), vbc_m);
	      scalF_m.initialize(getMesh(), getFieldLayout(), GuardCellSizes<Dim>(1), bc_m);
	    }
	    else {
	      vectF_m.initialize(getMesh(), getFieldLayout(), vbc_m);
	      scalF_m.initialize(getMesh(), getFieldLayout(), bc_m);
	    }
	  }
        }
        this->update();
	IpplTimings::stopTimer(this->myupdatTimer_m);
    }

    inline void particlePushTest() {
      /* 
	 push particles randomly 
      */
      IpplTimings::startTimer(this->parpushTimer_m);
      for (size_t i = 0; i< this->getLocalNum(); i++) {
	  for (int d = 0; d<3; d++)
		this->R[i](d) =  IpplRandom() * nr_m[d];
      }
      IpplTimings::stopTimer(this->parpushTimer_m);
    }
    
    inline void particleWalkTest() {
      /* 
	 push particles randomly 
      */
      IpplTimings::startTimer(this->parwalkTimer_m);
      for (size_t i = 0; i< this->getLocalNum(); i++) {
	  for (int d = 0; d<3; d++)
		this->R[i](d) +=  0.01 * IpplRandom() * nr_m[d];
      }
      IpplTimings::stopTimer(this->parwalkTimer_m);
    }

    inline void particleDriftTest() {
      /* 
	 push particles in direction (1,1,1)
      */
      
      //Inform dmsg("DriftTest", INFORM_ALL_NODES);
      //dmsg << "Vnode" << Ippl::myNode() << ' ' << this->getLocalNum() << " particles" << endl;
      IpplTimings::startTimer(this->pardrifTimer_m);
      
      if(Ippl::myNode() == 0)
		for (size_t i = 0; i< this->getLocalNum(); i++)
		{
				for (int d = 0; d<3; d++)
					this->R[i](d) +=  0.5;
		}
      	IpplTimings::stopTimer(this->pardrifTimer_m);

    }
    
    inline void particleOnEdgeTest()
    {
		/*
		 * one particle each on the "sides" of the region 
		 */
     IpplTimings::startTimer(this->paredgeTimer_m);
 
		NDRegion<double, Dim> localRegion =
		  ((*(getLayout().getLayout().begin_iv())).second)->getDomain();
		
		for (int i = 0; i < 6; i++)
		{
		  R[i](i%3) = (i%2) ? localRegion[i%3].min() : localRegion[i%3].max();
		}
		
		IpplTimings::stopTimer(this->paredgeTimer_m);

	}

private:

    inline void setBCAllOpen() {
        for (int i=0; i < 2*Dim; i++) {
            bc_m[i]  = new ZeroFace<double  ,Dim,Mesh_t,Center_t>(i);
            vbc_m[i] = new ZeroFace<Vector_t,Dim,Mesh_t,Center_t>(i);
            this->getBConds()[i] = ParticleNoBCond;
        }
    }

    inline void setBCAllPeriodic() {
        for (int i=0; i < 2*Dim; i++) {
            this->getBConds()[i] = ParticlePeriodicBCond;
            bc_m[i]  = new PeriodicFace<double  ,Dim,Mesh_t,Center_t>(i);
            vbc_m[i] = new PeriodicFace<Vector_t,Dim,Mesh_t,Center_t>(i);
        }
    }

    inline void setBCOOP() {
        for (int i=0; i < 2*Dim - 2; i++) {
            bc_m[i]  = new ZeroFace<double  ,Dim,Mesh_t,Center_t>(i);
            vbc_m[i] = new ZeroFace<Vector_t,Dim,Mesh_t,Center_t>(i);
            this->getBConds()[i] = ParticleNoBCond;
        }
        for (int i= 2*Dim - 2; i < 2*Dim; i++) {
            bc_m[i]  = new PeriodicFace<double  ,Dim,Mesh_t,Center_t>(i);
            vbc_m[i] = new PeriodicFace<Vector_t,Dim,Mesh_t,Center_t>(i);
            this->getBConds()[i] = ParticlePeriodicBCond;
        }
    }

    inline void gatherNGP() {
        // create interpolater object (nearest-grid-point method)
        IpplTimings::startTimer(ngpgathTimer_m); 
        E.gather(vectF_m, this->R, IntrplNGP_t());
	IpplTimings::stopTimer(ngpgathTimer_m); 
    }

    inline void gatherCIC() {
        // create interpolater object (cloud in cell method)
        IpplTimings::startTimer(cicgathTimer_m); 
        E.gather(vectF_m, this->R, IntrplCIC_t());
        IpplTimings::stopTimer(ngpgathTimer_m); 
    }

    inline void scatterCIC() {
        // create interpolater object (cloud in cell method)
        IpplTimings::startTimer(cicscatTimer_m); 
	qm.scatter(scalF_m, this->R, IntrplCIC_t());
        IpplTimings::stopTimer(cicscatTimer_m); 
    }

    inline void scatterNGP() {
        // create interpolater object (cloud in cell method)
        IpplTimings::startTimer(ngpscatTimer_m); 
	qm.scatter(scalF_m, this->R, IntrplNGP_t());
	IpplTimings::stopTimer(ngpscatTimer_m); 
    }

    /*  
	A Dim dimensional vector and boundary condition object field used to test operator
    */
    Field<Vektor<double,Dim>,Dim> vectF_m;
    BConds<Vector_t,Dim,Mesh_t,Center_t> vbc_m;

    /*  
	A Dim dimensional scalar field and boundary condition object used to test operator
    */
    Field<double,Dim> scalF_m;
    BConds<double,Dim,Mesh_t,Center_t> bc_m;

    /*  
	A Dim dimensional complex  field and boundary condition object used to test operator
    */
    Field<dcomplex,Dim> complF_m;

    /* 
       The FFT object
    */
    FFT_t *fft_m;

    /*
        Defines index space for the global domain
    */
    Vektor<int,Dim> nr_m;

    /*
        Defines mesh space for the global domain
    */
    Vector_t hr_m;

    /*
        Maximun and minimum extend of  the global domain
    */
    Vector_t rmax_m;
    Vector_t rmin_m;

    BC_t bco_m;
    InterPol_t interpol_m;
    bool fieldNotInitialized_m;
    bool withGuardCells_m;
    e_dim_tag decomp_m[Dim];

  
    double *partPerNode_m;
    double *globalPartPerNode_m;

    IpplTimings::TimerRef gradienTimer_m;
    IpplTimings::TimerRef divergeTimer_m;
    IpplTimings::TimerRef loadbalTimer_m;
    IpplTimings::TimerRef parpushTimer_m;
    IpplTimings::TimerRef myupdatTimer_m;
    IpplTimings::TimerRef ngpgathTimer_m; 
    IpplTimings::TimerRef cicgathTimer_m ;
    IpplTimings::TimerRef ngpscatTimer_m; 
    IpplTimings::TimerRef cicscatTimer_m ;
    IpplTimings::TimerRef pardrifTimer_m ;
    IpplTimings::TimerRef paredgeTimer_m ;
    IpplTimings::TimerRef parwalkTimer_m ;

};

int main(int argc, char *argv[]){
    Ippl ippl(argc, argv);
    Inform msg(argv[0]);
    Inform msg2all(argv[0],INFORM_ALL_NODES);

    static IpplTimings::TimerRef mainprgTimer = IpplTimings::getTimer("mainprgTimer");
    IpplTimings::startTimer(mainprgTimer);

    /************************************************************************************************************************************
	READ IN STUFF FROM CMD-LINE
    */

    Vektor<int,Dim> nr(atoi(argv[1]),atoi(argv[2]),atoi(argv[3]));

    const unsigned int totalP = atoi(argv[4]);
    const int nt              = atoi(argv[5]);

    InterPol_t myInterpol;
    if (string(argv[6])==string("CIC"))
        myInterpol = CIC;
    else
        myInterpol = NGP;

    BC_t myBC;
    if (string(argv[7])==string("OOO")) {
      myBC = OOO; // open boundary
    }
    if (string(argv[7])==string("OOP")) {
      myBC = OOP; // open boundary in x and y, periodic in z
    }
    if (string(argv[7])==string("PPP")) {
      myBC = PPP; // all periodic
    }

    bool gCells;
    gCells =  (string(argv[8])==string("GUARDCELLS"));

    e_dim_tag decomp[Dim];
    int serialDim = atoi(argv[9]);

    double threshold = atof(argv[10]); // threshold for load balancing in %

	Movement_t movementTest = PUSH;
	if(argc>11)
	{
		if (string(argv[11])==string("PUSH")) {
			movementTest = PUSH;
		}
		else if (string(argv[11])==string("DRIFT")) {
			movementTest = DRIFT;
		}
		else if (string(argv[11])==string("EDGE")) {
			movementTest = EDGE;
		}
		else if (string(argv[11])==string("WALK")) {
			movementTest = WALK;
		}

		else if (string(argv[11])==string("NONE")) {
			movementTest = NONE;
		}

    }

    /** 
	END READ IN STUFF FROM CMD-LINE
     ************************************************************************************************************************************/ 



    /************************************************************************************************************************************ 
	PRINTOUT CONFIGURATION
    */

    msg << "nt " << nt << " Np= " << totalP << " grid = " << nr;

    if (myInterpol==CIC)
      msg << " Cloud in cell (CIC) interpolation selected";
    else
      msg << " Nearest grid point (NGP) interpolation selected";

    if (gCells)
        msg << " Using guard cells" << endl;
    else
        msg << " Not using guard cells" << endl;

    if (myBC == OOO) {
      msg << "BC == OOO";
    }
    if (myBC == OOP) { 
      msg << "BC == OOP";
    }
    if (myBC == PPP) { 
      msg << "BC == PPP";
    }
    
    if (serialDim >= 3)
      msg << " 3D Domain decomposition " << endl;
    else 
      msg << " 2D Domain decomposition, serial dimension " << serialDim  << endl;
    
    msg << "Threshold for load balancing is " << threshold << " (%) based on the total number of particles" << endl;

    /** 
	END PRINTOUT CONFIGURATION
    ************************************************************************************************************************************/

    Mesh_t *mesh;
    FieldLayout_t *FL;
    ChargedParticles<playout_t>  *P;

    NDIndex<Dim> domain;
    if (gCells) {
        for(int i=0; i<Dim; i++)
            domain[i] = domain[i] = Index(nr[i] + 1);
    }
    else {
        for(int i=0; i<Dim; i++)
            domain[i] = domain[i] = Index(nr[i]);
    }

    for (int d=0; d < Dim; ++d)
        decomp[d] = (d == serialDim) ? SERIAL : PARALLEL;

    // create mesh and layout objects for this problem domain
    mesh          = new Mesh_t(domain);
    FL            = new FieldLayout_t(*mesh, decomp);
    playout_t* PL = new playout_t(*FL, *mesh);

    if (myBC==OOO)
        P = new ChargedParticles<playout_t>(PL,myBC,myInterpol,decomp,gCells);
    else {
        /*
          In case of periodic BC's define the domain with hr and rmin fix for all time steps
        */
        Vector_t hr(1.0);
        Vector_t rmin(0.0);
        Vector_t rmax(nr);
        P = new ChargedParticles<playout_t>(PL,myBC,myInterpol,hr,rmin,rmax,decomp,gCells);
    }

    // initialize the particle object: do all initialization on one node, and distribute to others
    double q = 1.0/totalP;
    unsigned long int nloc = totalP / Ippl::getNodes();
    double dummy = IpplRandom();

    P->create(nloc);
    for (unsigned long int i = 0; i< nloc; i++) {
      for (int d = 0; d<3; d++) 
	P->R[i](d) =  IpplRandom() * nr[d];
    }
    
    P->qm = q;

    /* 
	redistribute particles based on spatial layout and update all the fields (layout etc.)
    */

    P->myUpdate();

    INFOMSG(P->getMesh() << endl);
    INFOMSG(P->getFieldLayout() << endl);

    P->scatter();

    msg << endl << endl;

    //    P->initFields(); // Note: not used for the moment

    for (unsigned int it=0; it<nt; it++) {

      P->checkLoadBalance(threshold);
	
      /*
	After this we have set the scalar field
      */

      P->scatter();

      P->gradTest();

      P->divTest();

	  switch(movementTest)
		{
			case PUSH:
				P->particlePushTest();
				break;
      
			case DRIFT:
				P->particleDriftTest();
				break;
		
			case EDGE:
				P->particleOnEdgeTest();
				break;
		
			case WALK:
				P->particleWalkTest();
				break;
			
			default:
				break;
		}
		
      P->myUpdate();

      P->gather(it);
      
      if (!(it%10))
	msg << "Iter: " << it << " - min/max r and h " << P->getRMin()
	    << P->getRMax() << P->getHr() << endl;
	    
	
    }
    Ippl::Comm->barrier();
    
 
    msg << "Particle test ended." << endl;

    IpplTimings::stopTimer(mainprgTimer);
    IpplTimings::print();
    IpplTimings::print(string(argv[0])+string(".timing"));


	
	
    return 0;
}
/***************************************************************************
 * $RCSfile: addheaderfooter,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:17 $
 * IPPL_VERSION_ID: $Id: addheaderfooter,v 1.1.1.1 2003/01/23 07:40:17 adelmann Exp $
 ***************************************************************************/

