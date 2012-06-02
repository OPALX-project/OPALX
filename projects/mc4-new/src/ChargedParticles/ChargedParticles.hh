#ifndef CHARGED_PARTICLES_HH
#define CHARGED_PARTICLES_HH

#include "Ippl.h"
#include "FieldSolvers/FieldSolver.hh"

#include "Parser/FileStream.h"
#include "SimData.hh"
#include "PwrSpec/PwrSpec.hh"

template<class T, unsigned int Dim> class ChargedParticles;

#include "DataSink/DataSink.h"

#include <iostream>
#include <fstream>
#include <strstream>
#include <map>
#include "ListElem.h"

//corr stuff
#include <vector>
#include <ctime>
#include <cstdlib>
#include <cmath>
#include <algorithm>

static IpplTimings::TimerRef TIUpdate_m = IpplTimings::getTimer("Initial Update");



template <class T, unsigned int Dim>
class ChargedParticles : public ParticleBase< ParticleCashedLayout<T,Dim> > {
//class ChargedParticles : public ParticleBase< ParticleSpatialLayout<T,Dim> > {

private:

#ifdef IPPL_XT3
  size_t fragments_i;
  unsigned long total_free_i, largest_free_i, total_used_i;

  size_t fragments_f;
  unsigned long total_free_f, largest_free_f, total_used_f;
#endif

  ///mapping from species ID's to start index in particle container
  std::map<unsigned int, unsigned int> SpeciesIDs;
  ///current species index
  unsigned int SpeciesIdx;



public:

    typedef IntCIC IntrplCIC_t; 

    //    typedef typename ParticleSpatialLayout<T,Dim>::ParticlePos_t Ppos_t;
    // typedef typename ParticleSpatialLayout<T,Dim>::ParticleIndex_t PID_t;

    typedef UniformCartesian<Dim,T> Mesh_t;
    typedef typename ParticleCashedLayout<T,Dim>::SingleParticlePos_t Vector_t;
    //typedef typename ParticleSpatialLayout<T,Dim>::SingleParticlePos_t Vector_t;

    //switching to layout handling cashed particles around domain
    //needed for parallel halo finder
    typedef ParticleCashedLayout<T, Dim, Mesh_t>      Layout_t;
    //typedef ParticleSpatialLayout<T, Dim, Mesh_t>      Layout_t;

    typedef Vert                                       Center_t;

    typedef CenteredFieldLayout<Dim, Mesh_t, Center_t> FieldLayout_t;
    typedef Field<T, Dim, Mesh_t, Center_t>            Field_t;
    typedef Field<Vector_t, Dim, Mesh_t, Center_t>     VField_t;

    typedef Field<T, Dim, Mesh_t, Center_t>       RxField_t;

    typedef FFT<CCTransform, Dim, T>                   FFT_t;

    ParticleInteractAttrib<Vector_t>   V;      // particle velocity
    ParticleAttrib<Vector_t>   F;      // force field at particle position
    ParticleAttrib<T>          M;      // mass

    T ti_m;                            // actual integration time
    T dt_m;                            // timestep

    Vector_t rmin_m, rmax_m;

    Vector_t hr_m;

    Vektor<int,Dim> nr_m;

    BConds<T,Dim,Mesh_t,Center_t> bc_m;

    FieldSolver<T,Dim> *fs_m;

    bool fieldNotInitialized_m;

    SimData<T,Dim> simData_m;

    Vector_t boxSize_m;
    DataSink *itsDataSink_m;

    //constructor to create univ0
    ChargedParticles(Layout_t *playout, SimData<T,Dim> simData) : 
        ParticleBase< ParticleCashedLayout<T,Dim> >(playout),
        simData_m(simData)
    {
        SpeciesIdx = 0;

        // register the particle attributes
        addAttribute(M);
        addAttribute(V);
        addAttribute(F); 

        setBCAllPeriodic();

        // random number check
        IpplRandom.SetSeed(simData_m.iseed);
        T dummy = IpplRandom();

        fieldNotInitialized_m = true;
    }


    //constructor to create univ
    ChargedParticles(Layout_t *playout, SimData<T,Dim> simData, ParticleAttrib<Vector_t> &r ,
            ParticleAttrib<Vector_t>  &v,
            ParticleAttrib<T>    &m, size_t npLocal,
		     Vector_t boxSize, size_t chunksize, int iterations, size_t rest) : 
        ParticleBase< ParticleCashedLayout<T,Dim> >(playout),
        simData_m(simData)
    {
        SpeciesIdx = 0;
        
        Inform msg ("Create particles ");

        // register the particle attributes
        addAttribute(M);
        addAttribute(V);
        addAttribute(F); 

        setBoxSize(boxSize);

        setBCAllPeriodic();
        // random number check
        IpplRandom.SetSeed(simData_m.iseed);
        T dummy = IpplRandom();

        fieldNotInitialized_m = true;
	
	bool iCreate = (npLocal > 0);

	if (chunksize > 0) {    
            for (int i=0; i<iterations;i++) {
                if (iCreate) {
                    size_t offset = this->getLocalNum();
                    this->create(chunksize); 
                    for (size_t k = 0; k < chunksize; k++) {
		      if (isnan(r[k](0)) ||  isnan(r[k](1)) || isnan(r[k](2))) {
			; // skip 
		      }
		      else{ 
                        this->R[offset+k] = r[k];
                        this->V[offset+k] = v[k];
                        this->M[offset+k] = m[k];
		      }
                    }
                    r.destroy(chunksize,0,true);
                    v.destroy(chunksize,0,true);
                    m.destroy(chunksize,0,true);
		    
                }
                if (i==0)
                    rescaleAll();
                else
                    this->update();
		
                msg << "Step " << i << " np= " << scientific << (double) this->getTotalNum() << endl;
            }
            if (iCreate) {
                size_t offset = this->getLocalNum();
		rest  = r.size();
		msg << "Last Step rest= " << rest << endl;
		if (rest > 0) {
		  this->create(rest); 
		  for (size_t k = 0; k < rest; k++) {
                    this->R[offset+k] = r[k];
                    this->V[offset+k] = v[k];
                    this->M[offset+k] = m[k];
		  }
		  r.destroy(rest,0,true);
		  v.destroy(rest,0,true);
		  m.destroy(rest,0,true);
		}
            }

            if (this->getTotalNum() == 0)
	      rescaleAll();
            else
	      this->update();
	    
            msg << "Last Step (rest) done np= " << this->getTotalNum() << endl;
        }
        else {
            if (iCreate) {
                this->create(npLocal);
                this->R = r;
                this->V = v;
                this->M = m;
            }
            rescaleAll();
            r.destroy(npLocal,0,true);
            v.destroy(npLocal,0,true);
            m.destroy(npLocal,0,true);
        }
	itsDataSink_m = NULL; 
    }
    
    ChargedParticles(Layout_t *playout, SimData<T,Dim> simData, Vector_t boxSize) : 
        ParticleBase< ParticleCashedLayout<T,Dim> >(playout),
        simData_m(simData)
    {
        SpeciesIdx = 0;

        // register the particle attributes
        addAttribute(M);
        addAttribute(V);
        addAttribute(F); 

        setBoxSize(boxSize);

        setBCAllPeriodic();
        // random number check
        IpplRandom.SetSeed(simData_m.iseed);
        T dummy = IpplRandom();

        fieldNotInitialized_m = true;

        itsDataSink_m = NULL; 
    }


    ~ChargedParticles() {  

	if (itsDataSink_m)
	    delete itsDataSink_m;
    } 


    size_t getNeighborInformation(double bb) {
	
        return 0;

    }

    /// fortrans nint function
    inline T nint(T x)
    {
        return ceil(x + 0.5) - (fmod(x*0.5 + 0.25, 1.0) != 0);
    }

    /*
     *  scatter() 
     *
     *
     *  Created by Andreas Adelmann on 30.04.08.
     *  Copyright 2008 __MyCompanyName__. All rights reserved.
     *
     */
    void scatter() {
        fs_m->CICforwardAndCopy(this);
    }


    void addSpecies(ParticleAttrib<Vector_t> &r ,
		  ParticleAttrib<Vector_t>  &v,
		  ParticleAttrib<T>    &m) 
    {
       
        //assume we cannot add particles to previous species
        SpeciesIdx++;
        SpeciesIDs.insert( pair<unsigned int, unsigned int>(SpeciesIdx, this->getTotalNum()) );
        this->add(r, v, m);
	
    }

    void add(ParticleAttrib<Vector_t> &r ,
		  ParticleAttrib<Vector_t>  &v,
		  ParticleAttrib<T>    &m) 
    {
        	
        Inform msg ("add ");
        size_t curSize = this->getLocalNum();
        this->create(r.size()); 
        msg << "adding " << r.size() << " to local " << curSize << endl;
        for(int i=0; i<r.size(); i++) {
            this->R[curSize+i] = r[i];
            this->V[curSize+i] = v[i];
            this->M[curSize+i] = m[i];
        }
        rescaleAll(); 
        r.destroy(r.size(),0,true);
        v.destroy(v.size(),0,true);
        m.destroy(m.size(),0,true);
    }

    void tStep(DataSink *ds) {

        Inform msg ("tStep ");
        Inform msg2all ("tStep ",INFORM_ALL_NODES); 

        T omegab   = simData_m.deut / simData_m.hubble / simData_m.hubble;
        T omegatot = simData_m.omegadm + omegab;

        msg << "Begin timestepping simData_m.irun= " << simData_m.irun << endl;

        /* Rescaling of mesh needs to be taken care of in mc4.cpp */

        T ain  = 1.0 /(1+simData_m.zin);
        T afin = 1.0 /(1+simData_m.zfin);
        T apr  = 1.0 /(1+simData_m.zpr);
        T ainv = 1.0 /(simData_m.alpha);
        T pp   = std::pow(ain ,simData_m.alpha);
        T pfin = std::pow(afin,simData_m.alpha);

        // Energy test not implemented.
        if (simData_m.entest==0)
            msg << "Cosmic energy test not implemented!" << endl;

        if (simData_m.irun==0) {
            //----This is INTEG----------------------------------------
            T ppr  = std::pow(apr,simData_m.alpha);

            msg << "Time: a(t) = "<< std::pow(pp,ainv) << " z = " << 1.0/std::pow(pp,ainv)-1.0 << endl;

            T tau = (pfin-pp)/(simData_m.nsteps); //eps in integ - timestep

            //2nd order integrator: 
            if (simData_m.norder==2) {
                //norder locked to 2 for now, so this is not necessary, but..
                T tau1 = tau;
                T tau2 = 0.5*tau;
		
                msg << "Timesteping 2'd order  simData_m.nsteps= " << simData_m.nsteps << endl;
                for (int i = 0; i < simData_m.nsteps; i++) {
                    map1(tau2, omegatot, simData_m.alpha, pp);
                    pp=pp+tau2;  //first time update
                    map2(tau1, omegatot, simData_m.alpha, pp);
                    pp=pp+tau2; //second time update
                    map1(tau2, omegatot, simData_m.alpha, pp);
		
                    //Print positions and velocities, if required.
                    if (((i+1) % simData_m.iprint)==0) { 
                        //dump to text file if on one core with few particle 
                        if ((this->getTotalNum() < 260000) && (Ippl::getNodes()==1)) {
                            string fn1("mc4-out/step-");
                            ostringstream oss;
                            oss << fn1 << i+1;
			    
                            msg << "Time: a(t) = "<<std::pow(pp,ainv) << " z = " << 1.0/std::pow(pp,ainv)-1.0
                                << "\tdumping particle to file " << oss.str() << endl;
                            dumpParticles(string(oss.str()));      
                        }
                        if (ds) ds->writePhaseSpace(std::pow(pp, ainv), 1.0/std::pow(pp,ainv)-1.0, i);
                    } // end if iprint.
                    else
                        msg << "Time: a(t) = "<<std::pow(pp,ainv) << " z = " << 1.0/std::pow(pp,ainv)-1.0 << endl;
                }  //end for
            }  //end if norder==2
            //----End IEWTEG-------------------------------------------
        } //end if irun==0
        else { //use zprint to determine print intervals.
            //----This is INTEGZ-------------------------------------
            //Read zprnt file
            T zprnt[simData_m.irun];
            T pprnt[simData_m.irun+2];
            if (isRoot()) {
                ifstream ifstr;
                ifstr.open(string("zprnt").c_str(),ios::in);
                ifstr.precision(8);
                for (int i=0; i<simData_m.irun; i++ ) 
                    ifstr >> zprnt[i];
                ifstr.close();
            } //end zprnt file read 

            // Work out print intervals.
            T prange = pfin - pp;
            pprnt[0]=pp;
            pprnt[simData_m.irun+1]=pfin;
            for (int i=0; i<simData_m.irun; i++) 
                pprnt[i+1]=std::pow((T)(1.0+zprnt[i]),(T)(-ainv));

            int nsum = 0;
            T ptime = pp;
            for (int i=1; i<simData_m.irun+1; i++) { //<---F90 indices corrected here.
                T delp = pprnt[i]-pprnt[i-1];
                int nsloc = nint((delp/prange)*simData_m.nsteps);
                T epsloc = delp/nsloc;
                nsum += nsloc;
                //2nd order integrator - input locked to norder = 2 atm.
                if (simData_m.norder == 2) {
                    T tau1 = epsloc ;
                    T tau2 = 0.5*epsloc ;
                    for (int j=0; j<nsloc; j++) {
                        map1(tau2, omegatot, simData_m.alpha, pp);
                        pp=pp+tau2;  //first time update
                        map2(tau1, omegatot, simData_m.alpha, pp);
                        pp=pp+tau2; //second time update
                        map1(tau2, omegatot, simData_m.alpha, pp);
                    } // end for j
                    /* TO DO:
                    //Particle print out at specified zs.
                    // Print out 2D density field at specified zs.
                    */
                } //end if norder == 2
            } //end for i
            //----End INTEGZ------------------------------------------
        } //end else (i.e. irun!=0)
    } //end tStep 


    void genLocalUniformICiter(T scale) {

        Inform msg ("genUniformIC ");

        size_t count = 0;
        //size_t numParticles = (size_t)(simData_m.np*simData_m.np*simData_m.np)/Ippl::getNodes();
        size_t numParticles = simData_m.np;
        numParticles *= simData_m.np;
        numParticles *= simData_m.np;
        numParticles /= Ippl::getNodes();
        NDIndex<3> lDomain = getFieldLayout().getLocalNDIndex();
        NDIndex<Dim> domain = getFieldLayout().getDomain();

        msg << "Start generating uniform initial conditions scale=" <<  scale << endl;

        IpplRandom.SetSeed(simData_m.iseed);

        T tmp = IpplRandom();

        msg << "start loop creating " << numParticles << " particles" << endl;
        for(count=0; count < numParticles; count++) {
            T r1 = IpplRandom();
            T r2 = IpplRandom();
            T r3 = IpplRandom();
            T r4 = IpplRandom();
            T r5 = IpplRandom();
            T r6 = IpplRandom();

            this->create(1);
            this->R[count](0) = scale*(lDomain[0].first() + r1*(lDomain[0].length()))/domain[0].length();
            this->R[count](1) = scale*(lDomain[1].first() + r2*(lDomain[1].length()))/domain[1].length();
            this->R[count](2) = scale*(lDomain[2].first() + r3*(lDomain[2].length()))/domain[2].length();
            this->V[count](0) = r4;
            this->V[count](1) = r5;
            this->V[count](2) = r6;
            this->M[count] = 1.0;

        }
        msg << "after loop" << endl;
	IpplTimings::startTimer(TIUpdate_m);
        this->update();
	IpplTimings::stopTimer(TIUpdate_m);
        msg << "Particles " << this->getTotalNum() << " created" << endl;
    }

    inline const Mesh_t& getMesh() const { return this->getLayout().getLayout().getMesh(); }

    inline Mesh_t& getMesh() { return this->getLayout().getLayout().getMesh(); }

    inline const FieldLayout_t& getFieldLayout() const {
        return dynamic_cast<FieldLayout_t&>(this->getLayout().getLayout().getFieldLayout());
    }

    inline FieldLayout_t& getFieldLayout() {
        return dynamic_cast<FieldLayout_t&>(this->getLayout().getLayout().getFieldLayout());
    }

    inline void set_FieldSolver(FieldSolver<T,Dim> *fs) { fs_m = fs; }

    inline T getTime() {return ti_m; }
    inline void setTime(T ti) {ti_m = ti; }

    inline T getdT() {return dt_m; }
    inline void setdT(T dt) {dt_m = dt; }

    inline Vector_t getGridSize() { return nr_m; }
    inline Vector_t getMeshSpacing() { return hr_m; }
    inline Vector_t getRmin() { return rmin_m; }
    inline Vector_t getRmax() { return rmax_m; }

    inline bool isRoot() { return (Ippl::Comm->myNode() == 0); }

    inline void     setBoxSize(Vector_t bs) { boxSize_m = bs;}
    inline Vector_t getBoxSize() { return boxSize_m;}

    void rescaleAll() {  

        Inform msg ("rescaleAll ");

        rmin_m = Vector_t(0.0);
        // XXX IMPORTANT
        // Particles always have to be in the range
        //      [0,ng_comp] x [0,ng_comp] x [0,ng_comp]
        // scale accordingly!
        rmax_m = getBoxSize();

        Vector_t len;
        len = rmax_m - rmin_m;

        NDIndex<Dim> domain = getFieldLayout().getDomain(); 
        for(int i=0; i<Dim; i++) {
            nr_m[i] = domain[i].length();
            hr_m[i] = (len[i] / (nr_m[i]-1));
        }

        msg << "rmax and rmin (is boxsixe) " << rmax_m << rmin_m << " nr_m= " << nr_m << " hr_m= " << hr_m << endl;

        getMesh().set_meshSpacing(&(hr_m[0]));
        getMesh().set_origin(rmin_m);
        setBCAllPeriodic();
        Ippl::Comm->barrier();
        this->update();
    }


    inline void do_binaryRepart() {
        BinaryRepartition(*this);
        boundp();
    }

    inline void setBCAllPeriodic() {
        for (int i=0; i < 2*Dim; ++i) {
            bc_m[i] = new ParallelPeriodicFace<T,Dim,Mesh_t,Center_t>(i);
            this->getBConds()[i] = ParticlePeriodicBCond;
        } 
    }

    void readInputDistribution(string fn) {

        Inform msg2all("readInputDistribution ", INFORM_ALL_NODES);
        Inform msg ("readInputDist ");

        //if file ends with ".h5" read in h5 file
        string ending = ".h5";
        if(fn.length() > ending.length() && (0 == fn.compare (fn.length() - ending.length(), ending.length(), ending))) {
            itsDataSink_m = new DataSink(fn, this, 0);
        }
        else {
            if (isRoot()) {
                T  arg1, arg2, arg3, arg4, arg5, arg6;
                T maxx = 0.0;
                T maxy = 0.0;
                T maxz = 0.0;
                ifstream *ifstr = new ifstream;
                ifstr->open(fn.c_str(), ios::in );
                size_t count = 0;

                *ifstr >> arg1 >> arg2 >> arg3 >> arg4 >> arg5 >> arg6;

                while(*ifstr) {
                    this->create(1);
                    this->R[count](0) = arg1;
                    maxx = max(arg1, maxx);
                    this->R[count](1) = arg3;
                    maxy = max(arg3, maxy);
                    this->R[count](2) = arg5;
                    maxz = max(arg5, maxz);
                    this->V[count](0) = arg2;
                    this->V[count](1) = arg4;
                    this->V[count](2) = arg6;
                    this->M[count] = 1.0;
                    count++;
                    *ifstr >> arg1 >> arg2 >> arg3 >> arg4 >> arg5 >> arg6;
                }
                ifstr->close();
            }
        }
        rescaleAll();
        msg2all << "Particles created: " << this->getLocalNum() << endl;
    }

    void dumpParticles(string fname) {

        Inform msg("dumpParticles ");

        unsigned int pwi = 14;
        ofstream of;
        unsigned int ok = 0;
        int toSave=0;
        if (this->getTotalNum() < 260000) {
            int tag1 = IPPL_APP_TAG1;
            int tag2 = IPPL_APP_TAG2;
            stringstream fn;
            fn << fname;
            fn << "-P";
            fn << Ippl::getNodes();
            fn << ".dat";

            double grid2box = simData_m.rL/simData_m.ng_comp;

            if (this->singleInitNode()) {
                of.open(fn.str().c_str(),ios::out);
                for (size_t ii = 0 ; ii < this->getLocalNum(); ++ii) {
                    for (int jj = 0 ; jj < Dim ; ++jj) {
                        of << setw(pwi) << this->R[ii](jj)*grid2box << ' ';
                        of << setw(pwi) << this->V[ii](jj) << ' ';
                    }
                    of << endl;
                }
                of.close(); 

                for (toSave=1; toSave <= Ippl::getNodes()-1; toSave++) {
                    // notify toSave that he can save
                    Message *msg = new Message();
                    msg->put(toSave);
                    bool res = Ippl::Comm->send(msg, toSave, tag1); 
                    if (! res) 
                        ERRORMSG("Ippl::Comm->send(smsg, 0, tag1) failed " << endl;);
                    // wail til toSave reports done with save
                    Message* doneMsg;
                    doneMsg = Ippl::Comm->receive_block(toSave,tag2);
                    doneMsg->get(&ok);
                    delete doneMsg;
                    if (ok != toSave) {
                        ERRORMSG("Node0 got wrong msg from node " << toSave);
                        ERRORMSG("received ID = " << ok  << endl);
                    }
                }
            }
            else {
                // Wail til we get ready to save
                Message* readyToSaveMesg;
                int z=0;

                readyToSaveMesg = Ippl::Comm->receive_block(z, tag1);
                readyToSaveMesg->get(&ok);
                delete readyToSaveMesg;
                if (ok == Ippl::myNode()) {
                    of.open(fn.str().c_str(),ios::app);
                    for (size_t ii = 0 ; ii < this->getLocalNum(); ++ii) {
                        for (int jj = 0 ; jj < Dim ; ++jj) {
                            of << setw(pwi) << this->R[ii](jj)*grid2box << ' ';
                            of << setw(pwi) << this->V[ii](jj) << ' ';
                        }
                        of << endl;
                    }
                    of.close();
                    // Notify node0 we are done by sending the node id
                    Message *msg = new Message();
                    msg->put(ok);
                    bool res = Ippl::Comm->send(msg, 0, tag2);
                    if (! res) 
                        ERRORMSG("Ippl::Comm->send(smsg, 0, tag2) failed " << endl;);
                }
                else {
                    ERRORMSG("Node" << Ippl::myNode() << "received ID = " << ok  << endl);
                }
            }
            Ippl::Comm->barrier(); 
        }
    }


    inline void writeDiagnostics() {


    } 

    inline void map2(T tau, T omegatot, T alpha, T pp) {
        /// adot in units of H_0^-1
        Inform msg("map2: ");
        if (fs_m) {
            // Uncomment out below and fix CIC in computeForceField
            fs_m->computeForceField(this);

            T adot = sqrt((omegatot+(1.0-omegatot)*std::pow(pp,((T)3.0/alpha))) / std::pow(pp,(T)(1.0/alpha)) );
            T phiscale = 1.5*omegatot/std::pow(pp,(T)(1.0/alpha));
            T fscale   = phiscale/(alpha*adot*std::pow(pp,(T)(1.0 - 1.0/alpha)));
            INFOMSG("adot= " << adot << " phiscale= " << phiscale << " fscale= " << fscale << " tau= " << tau << endl);
            assign(V, V + fscale*tau*F);
        }
    }  

    inline void map1(T tau, T omegatot, T alpha, T pp) {
        Inform msg("map1: ");
        /// adot in units of H_0^-1
        T adot = sqrt((omegatot+(1.0-omegatot)*std::pow(pp,((T)3.0/alpha))) / std::pow(pp,(T)(1.0/alpha)) );
        T pre  = 1.0/( alpha*adot*std::pow(pp,(T)(1.0 + 1.0/alpha)));

        assign(this->R, this->R + tau * pre * this->V);
        this->update();
    }

    void corr2(int N, vector<double> x, vector<double> y, vector<double> z, vector<double> &r_weight, vector<double> &corr_f)
    {
        double h=0.7, redshift=0.0, box_size=simData_m.rL;
        double r_min=0.02, r_max=box_size/2.0*1.732051;
        double pi=3.14159265;
        int Nbins = corr_f.size() + 1;

        double rr, dr, dV, norm, d_max;
        vector<double> r(N);
        vector<int> count(Nbins-1);

        count.assign(Nbins-1, 0);
        r_weight.assign(r_weight.size(), 0.0);

        d_max = 0.0;
        dr = log10(r_max/r_min)/(Nbins-1);

        //XXX: for cross correlation i = light and j = heavy (i.e)
        for(int i=0; i<N; i++) {
            for(int j=0; j<N; j++) {
                r[j] = sqrt(pow((min(abs(x[i]-x[j]), abs(abs(x[i]-x[j])-box_size))),2.0) +
                        pow((min(abs(y[i]-y[j]), abs(abs(y[i]-y[j])-box_size))),2.0) +
                        pow((min(abs(z[i]-z[j]), abs(abs(z[i]-z[j])-box_size))),2.0));
            }

            std::sort(r.begin(), r.end());

            if(r[N-1] > d_max) 
                d_max = r[N-1];

            // bin the data
            // not accounting for pairing with itself
            int k = 1;
            for(k; k<N; k++) {
                if(r[k] >= r_min) 
                    break;
            }

            rr = r_min;
            for(int j=0; j<Nbins-2; j++) {
                for(k; k<N; k++) {
                    if(r[k] > rr*pow(10,dr)) 
                        break;
                    count[j]++;
                    r_weight[j] += r[k];
                }
                rr = rr * pow(10,dr);
            }

            //FIXME: point of synchronization
            //reduce r_weight,count (sum all)
        
        }

        cout << "largest distance between two objects :" << d_max << endl;
        cout << "maximum possible distance :" << r_max << endl;
        /*cout << "total number of pairs :" << sum(count(:))/2.0, (N**2-N)/2.0 << endl;*/
        /*cout << "difference :" << sum(count(:)) - (N**2-N) << endl;*/

        for(int j=0; j<Nbins-2; j++) {
            if(count[j] > 0)
                r_weight[j] = r_weight[j]/count[j];
            else
                r_weight[j] = (r_min*pow(pow(10,dr),j) + r_min*pow(pow(10,dr),j+1))/2.0;
            count[j] /= 2;
        }
        
        rr = r_min;
        for(int i=0; i<Nbins-2; i++) {
            dV = 4.0*pi/3.0 * (pow(rr*(pow(10,dr)),3.0)-pow(rr,3.0));
            norm = 1.0/dV * pow(box_size,3.0)/(N*(N-1)/2.0);
            corr_f[i] = count[i]*norm;
            if(r_weight[i] < rr && r_weight[i] > rr*pow(10,dr)){
                cout << "PROBLEM WITH R IN BIN " << i << ":";
                /*cout << rr, rr*10**dr, r_weight(i), corr_f(i)*/
                exit(-1);
            }
            rr = rr*pow(10,dr);
        }

        return;
    }

    void ran_cat(int N, vector<double> &ran_corr)
    {
        double box_size=256.0;
        int Nbins = ran_corr.size()+1;

        vector<double> x(N), y(N), z(N);
        vector<double> r_weight(Nbins-1);

        cout << "Building random catalog with " << N << " particles" << endl;

        //generate random number between [0,1]
        srand(time(NULL));
        for(int i=0; i<N; i++) {
            x[i] = (rand() / ((double)RAND_MAX+1))*box_size;
            y[i] = (rand() / ((double)RAND_MAX+1))*box_size;
            z[i] = (rand() / ((double)RAND_MAX+1))*box_size;
        }

        corr2(N, x, y, z, r_weight, ran_corr);

        return;
    }

    //FIXME: parallelize!
    void calcAndDumpCorrelation(int Nbins=31) {

        size_t N = this->getTotalNum();
        //FIXME: this should be 2*N or 4*N (if sample is large, maybe =N)
        int Nran=2*N;
        vector<double> x(N), y(N), z(N);
        vector<double> r_weight(Nbins-1), ran_corr(Nbins-1);
        vector<double> r_corr(Nbins-1), hh_corr(Nbins-1);

        for(int i=0; i<N; i++) {
            x[i] = this->R[i](0)*simData_m.rL/simData_m.ng_comp;
            y[i] = this->R[i](1)*simData_m.rL/simData_m.ng_comp;
            z[i] = this->R[i](2)*simData_m.rL/simData_m.ng_comp;
        }

        corr2(N, x, y, z, r_weight, hh_corr);

        ran_cat(Nran, ran_corr);

        ofstream corrout;
        corrout.open("mc4.corr", ios::out);
        for(int i=0; i<Nbins-2; i++) {
            hh_corr[i] = hh_corr[i]/ran_corr[i]-1.0;
            r_corr[i] = r_weight[i];
            corrout << r_corr[i] << "\t" << hh_corr[i] << endl;
        }
        corrout.close();
    }


    ///SERIAL at the moment
    void dumpVTK(unsigned int timestep) {

        ofstream vtkout;
        vtkout.precision(10);
        vtkout.setf(ios::scientific, ios::floatfield);

        std::stringstream fname;
        fname << "vtkdata/step_";
        fname << setw(4) << setfill('0') << timestep;
        fname << ".vtk";

        vtkout.open(fname.str().c_str(), ios::out);
        vtkout << "# vtk DataFile Version 2.0" << endl;
        vtkout << "MC4 position, mass and velocity data on unstructured grid" << endl;
        vtkout << "ASCII" << endl;
        vtkout << "DATASET UNSTRUCTURED_GRID" << endl;
        vtkout << "POINTS " << this->getTotalNum() << " FLOAT" << endl;

        size_t numParticles = this->getTotalNum();
        for (int i=0; i < numParticles; i++)
            vtkout << this->R[i](0) << "\t" << this->R[i](1) << "\t" << this->R[i](2) << endl;
        vtkout << endl;
        
        vtkout << "CELLS " << numParticles << " " << 2*numParticles << endl;
        for (int i=0; i < numParticles; i++)
            vtkout << "1 " << i << endl;
        vtkout << endl;

        vtkout << "CELL_TYPES " << numParticles << endl;
        for(int i=0; i < numParticles; i++)
            vtkout << "2" << endl;
        vtkout << endl;

        vtkout << "POINT_DATA " << numParticles << endl;
        vtkout << "SCALARS mass float 1" << endl
               << "LOOKUP_TABLE default" << endl;
        for (int i=0; i < numParticles; i++)
            vtkout << M[i] << endl;
        vtkout << endl;

        vtkout << "VECTORS V float" << endl;  
        for (int i=0; i < numParticles; i++)
            vtkout << V[i](0) << "\t" << V[i](1) << "\t" << V[i](2) << endl;

        vtkout.close();

    }

    /**
    Data Format [http://t8web.lanl.gov/people/heitmann/arxiv/data.html]
    There are two files in the archives: particles_name, which has the following content:

        x[Mpc], v_x[km/s], y[Mpc], v_y[km/s], z[Mpc], v_z[km/s], particle mass[M_sol], particle tag

    This file is written in single precision and binary format. The files were written in
    "direct access format", i.e., no extra 4 byte from FORTRAN. The byteorder is little-endian,
    which will be fine on Intels and Alphas.
    */
    void dumpCosmo(unsigned int timestep) {

        ofstream out;
        std::stringstream fname;
        double grid2box = simData_m.rL/simData_m.ng_comp;

        fname << "vtkdata/step_";
        fname << setw(4) << setfill('0') << timestep;
        fname << ".cosmo";

        out.open(fname.str().c_str(), ios::out|ios::binary);

        for(int i=0; i<this->getTotalNum(); i++) {

            //make sure we dump float values
            float x = this->R[i](0)*grid2box;
            float y = this->R[i](1)*grid2box;
            float z = this->R[i](2)*grid2box;
            float vx = this->V[i](0);
            float vy = this->V[i](1);
            float vz = this->V[i](2);
            float m = this->M[i];

            out.write((char *)&x, sizeof(float));
            out.write((char *)&vx, sizeof(float));
            out.write((char *)&y, sizeof(float));
            out.write((char *)&vy, sizeof(float));
            out.write((char *)&z, sizeof(float));
            out.write((char *)&vz,sizeof(float));
            out.write((char *)&m, sizeof(float));
            out.write((char *)&i, sizeof(int));

        }

        out.close();
        
        std::stringstream fname1;
        fname1 << "vtkdata/step_";
        fname1 << setw(4) << setfill('0') << timestep;
        fname1 << "_light.cosmo";

        out.open(fname1.str().c_str(), ios::out|ios::binary);

        for(int i=0; i<this->getTotalNum(); i+=2) {

            //make sure we dump float values
            float x = this->R[i](0)*grid2box;
            float y = this->R[i](1)*grid2box;
            float z = this->R[i](2)*grid2box;
            float vx = this->V[i](0);
            float vy = this->V[i](1);
            float vz = this->V[i](2);
            float m = this->M[i];

            out.write((char *)&x, sizeof(float));
            out.write((char *)&vx, sizeof(float));
            out.write((char *)&y, sizeof(float));
            out.write((char *)&vy, sizeof(float));
            out.write((char *)&z, sizeof(float));
            out.write((char *)&vz,sizeof(float));
            out.write((char *)&m, sizeof(float));
            out.write((char *)&i, sizeof(int));

        }

        out.close();
        
        std::stringstream fname2;
        fname2 << "vtkdata/step_";
        fname2 << setw(4) << setfill('0') << timestep;
        fname2 << "_heavy.cosmo";

        out.open(fname2.str().c_str(), ios::out|ios::binary);

        for(int i=1; i<this->getTotalNum(); i+=2) {

            //make sure we dump float values
            float x = this->R[i](0)*grid2box;
            float y = this->R[i](1)*grid2box;
            float z = this->R[i](2)*grid2box;
            float vx = this->V[i](0);
            float vy = this->V[i](1);
            float vz = this->V[i](2);
            float m = this->M[i];

            out.write((char *)&x, sizeof(float));
            out.write((char *)&vx, sizeof(float));
            out.write((char *)&y, sizeof(float));
            out.write((char *)&vy, sizeof(float));
            out.write((char *)&z, sizeof(float));
            out.write((char *)&vz,sizeof(float));
            out.write((char *)&m, sizeof(float));
            out.write((char *)&i, sizeof(int));

        }

        out.close();
    }


    void print(ostream &os) 
    {
        os << " ------ Conditions @ t= " << getTime() << " ------ " << endl;
        os << " N= " << this->getTotalNum() << " Grid= " << getGridSize() << endl;
        os << " h= " << getMeshSpacing() << " []" << endl ;
        os << " rmax=" << getRmax() << " []" << endl;
        os << " rmin=" << getRmin() << " []" << endl;
    }

};

inline ostream &operator<<(ostream &os, ChargedParticles<double,3> &data)
{
	data.print(os);
	return os;
}

#endif 

