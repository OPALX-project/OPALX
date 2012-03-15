// ------------------------------------------------------------------------
// $RCSfile: ChargedParticles.hh,v $
// $Header: /afs/psi.ch/user/a/adelmann/private/cvsroot/pcyclint/prog/ChargedParticles/ChargedParticles.hh,v 1.7 2004/10/06 15:23:32 adelmann Exp $
// ------------------------------------------------------------------------
// $Revision: 1.7 $
// $State : $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
// Diption:
//
// $Date: 2004/10/06 15:23:32 $
// $Author: adelmann $
// $Log: ChargedParticles.hh,v $
// Revision 1.7  2004/10/06 15:23:32  adelmann
// Add a function clean do delete particles which are to
// far away.
//
// Add function getMinPartPerNode()
//
// Revision 1.6  2004/10/05 12:00:47  adelmann
// add maximum inbalance to control binary repartition
// add to writeStatistics some information on computational issues (inbalance)
// and write to a sepatate file.
// Fix Lorenty transformation
//
// Revision 1.5  2004/10/01 20:36:48  adelmann
// -Remove GAUSS and UNIFORM distribution
// -fix bug in convert from/to CU and mks, the direction of the
//  momenta was list
// - local number is now correct after boundp
// - write out momenta in m/s (verlocity)
// - only one writePhaseSpace the one with the clone argument
//
// Revision 1.4  2004/09/29 18:51:59  adelmann
// Version with new distribution generation. Distribution without
// space charge goes til 72 MeV, see logfile
//
// Revision 1.3  2004/09/28 04:41:30  adelmann
// Add  many structures for cloning
//
// Revision 1.2  2004/09/24 19:28:46  adelmann
// Add bunchNo, all bunches with bunchNo <> 0 are clones!
//
// Change scaling of p in write phase space
//
// Revision 1.1.1.1  2004/09/22 12:10:45  adelmann
// Imported sources pcyclubt
//
// Revision 1.3  2003/10/23 04:44:48  adelmann
// getRMean recomputes the mean value
//
// Revision 1.2  2003/10/07 09:03:53  adelmann
// writeStatistics: add boolean calcTune to write a separate file with
// informations for LOMB tune analysis
//
// Revision 1.1.1.1  2003/10/03 11:55:47  wittberger
// Marcus Wittbergers Summerwork 2003
//
// ------------------------------------------------------------------------
#ifndef CHARGED_PARTICLES_HH
#define CHARGED_PARTICLES_HH

#include "Ippl.h"

#ifdef USE_PBE
#include "pbe.h"
#endif

#include <iostream>
#include <fstream>
#include <strstream>
#include "Const.hh"


template <class T, unsigned int Dim>
class ChargedParticles : public ParticleBase< ParticleSpatialLayout<T,Dim> > {
  
private:
    
    T coupling_m;

    unsigned long N0_m;
    double current_m;    

    unsigned int numberOfClones;

public:
    
    typedef IntCIC IntrplCIC_t; 
    typedef IntNGP IntrplNGP_t;
    
    typedef typename ParticleSpatialLayout<T,Dim>::ParticlePos_t Ppos_t;
    typedef typename ParticleSpatialLayout<T,Dim>::ParticleIndex_t PID_t;
    
    typedef UniformCartesian<Dim,T> Mesh_t;
    
    typedef typename ParticleSpatialLayout<T,Dim>::SingleParticlePos_t Vector_t;
    
    typedef ParticleSpatialLayout<T, Dim, Mesh_t>      Layout_t;
  
    typedef Cell                                       Center_t;
    
    typedef CenteredFieldLayout<Dim, Mesh_t, Center_t> FieldLayout_t;
    typedef Field<T, Dim, Mesh_t, Center_t>            Field_t;
    typedef Field<Vector_t, Dim, Mesh_t, Center_t>     VField_t;
    
    ParticleAttrib<T>          q_m; // charge (=mass if gravity)
    ParticleAttrib<T>          mass_m; //  mass
    ParticleAttrib<T>          h; // stepsize
    ParticleAttrib<T>          Ti;
    ParticleAttrib<T>          tact; // actual time
    ParticleAttrib<Vector_t>   P;    // particle velocity 
    ParticleAttrib<Vector_t>   Ef;   // electric field at particle position
    ParticleAttrib<int>        bunchNo; // which bunch am  I

    Vector_t rmin_m, rmax_m;
    Vector_t CenterOfMass_m;
    
    inline T getq() { return q_m[0]; }
    inline T getmass() { return mass_m[0]; }
    
/** Clone business
    
*/

    Vector_t MeanR_m[MAXCLONE];
    Vector_t MeanR_loc[MAXCLONE];

    Vector_t MeanP_m[MAXCLONE];
    Vector_t MeanP_loc[MAXCLONE];

    double gamma_m[MAXCLONE];

    unsigned long int locNum_m[MAXCLONE];
    unsigned long int totNum_m[MAXCLONE];
    

    inline void setNumOfClones(unsigned int nc) {numberOfClones=nc;}
    inline unsigned int getNumOfClones() {return numberOfClones;}
    
    inline Vector_t getMeanR() { 
        /** assume all clones have the same number of particles
            and they are equal the number of particles of the mother bunch
        */
        double den = 1.0/(getTotalNum()/(numberOfClones+1));
  
        for (int j=0; j<MAXCLONE; j++) 
            MeanR_loc[j] = MeanR_m[j] = 0.0;
        
        for (unsigned long int i=0; i<getLocalNum(); i++) 
            MeanR_loc[bunchNo[i]] += R[i]*den;
        reduce(MeanR_loc, MeanR_loc + MAXCLONE, MeanR_m, OpAddAssign());
        return MeanR_m[0];
    }
    /** assume a previous call to  getMeanR() */
    inline Vector_t getMeanR(int k) { return MeanR_m[k]; }

    inline Vector_t getMeanP() { 
        /** assume all clones have the same number of particles
            and they are equal the number of particles of the mother bunch
        */
        double den = 1.0/(getTotalNum()/(numberOfClones+1));
  
        for (int j=0; j<MAXCLONE; j++) 
            MeanP_loc[j] = MeanP_m[j] = 0.0;
        
        for (unsigned long int i=0; i<getLocalNum(); i++) 
            MeanP_loc[bunchNo[i]] += P[i]*den;
        reduce(MeanP_loc, MeanP_loc + MAXCLONE, MeanP_m, OpAddAssign());
        return MeanP_m[0];
    }
    /** assume a previous call to  getMeanP() */
    inline Vector_t getMeanP(int k) { return MeanP_m[k]; }


    void compGamma() {
        for (int i=0; i< MAXCLONE; i++) {
            double b2 = dot(MeanP_m[i],MeanP_m[i])/(CLIGHT*CLIGHT);
            gamma_m[i] = 1.0/sqrt(1.0-b2);
        }
    }

    inline T getGamma(int k) { 
        return gamma_m[k];
    }

    inline T getGamma() { 
        return gamma_m[0];
    }
    
    inline T getBeta() { 
        return sqrt(1.0-(1.0/(gamma_m[0]*gamma_m[0])));
    }

    void updateBunches() {        
        /** 
            Must be done after each integration step
            in consequence all getXXX(k) are valid
        */
        
        Vector_t cm = getMeanP();
        Vector_t cv = getMeanR();
        compGamma();
        
        for (unsigned long int i=0; i<MAXCLONE ; i++) 
            locNum_m[i] = totNum_m[i] = 0;
        
        if (getNumOfClones()==0) {
            locNum_m[0]=getLocalNum();
            totNum_m[0]=getTotalNum();
        }
        else {
            for (unsigned long int i=0; i<getLocalNum(); i++) 
                locNum_m[bunchNo[i]]++;
            reduce(locNum_m,locNum_m+MAXCLONE,totNum_m,OpAddAssign());
        }
    }
    

    inline unsigned long int getLocalnum(unsigned int clone) {return locNum_m[clone];}
    inline unsigned long int getTotalnum(unsigned int clone) {return totNum_m[clone];}

    inline T getTime(int clone) {
        unsigned nPartLoc = 0;
        double ti = 0.0;
        
        for (unsigned k=0; k<getLocalNum(); k++) {
            if (bunchNo[k] == clone) {
                nPartLoc++;
                ti += Ti[k];
            }
        }
        return ti/nPartLoc;
    }
    
    inline void setTime(int clone, double t) {
        for (unsigned k=0; k<getLocalNum(); k++) 
            if (bunchNo[k] == clone)
                Ti[k] = t;
    }
    
    inline Vector_t getCenterOfMass() { return CenterOfMass_m; }
    
    inline Vector_t getGridSize() { return nr_m; }
    inline Vector_t getMeshSpacing() { return hr_m; }

    inline Vector_t getRmin() { return rmin_m; }
    inline Vector_t getRmax() { return rmax_m; }

    inline void setCurrent(double I) { current_m = I; }
    inline double getCurrent() { return  current_m; }    

    inline void setN0() {N0_m = getTotalNum();}
    inline void setN0(unsigned long n) {N0_m = n;}
    inline unsigned long getN0() {return N0_m;}
    
    inline bool isRoot() { return (Ippl::Comm->myNode() == 0); }


/**
   Convertion stuff from MCU to mks and vice versa
*/

    inline double mcu2mks(double p) {
        double gamma2 = pow(p*1.0e-3,2) + 1.0;     // MCU*1.0e-3 -> CU
        double beta = sqrt(1.0 - (1.0/gamma2));
        if (p<0.0)
            return  -1.0*CLIGHT*beta;
        else
            return  CLIGHT*beta;
    }
    
    inline double cu2mks(double p) {
        double gamma2 = pow(p,2) + 1.0;     
        double beta = sqrt(1.0 - (1.0/gamma2));
        if (p<0.0)
            return  -1.0*CLIGHT*beta;
        else
            return  CLIGHT*beta;
    }
    
    inline double mks2mcu(double p) {
        double beta = p/CLIGHT;
        double gamma2 = 1.0 / (1.0 - (beta*beta));
        if(p<0.0)
            return  -1.0e3*sqrt(gamma2-1.0);
        else
            return  1.0e3*sqrt(gamma2-1.0);
    }
    

    // Grid stuff
    Field_t  rho_m;                       // charge density/electric potential
    VField_t  eg_m;                       // electric field on grid
    
    BCT bcType_m;                        // type of BC
   
    Vector_t hr_m;
    Vektor<int,Dim> nr_m;
    
    T act_s;
    T act_phi;
    T frequ_m;
    
    BConds<T,Dim,Mesh_t,Center_t> bc_m;
    BConds<Vector_t,Dim,Mesh_t,Center_t> vbc_m;
  
    

    
    inline T   getFrequ() { return frequ_m; }
    inline void setFrequ(T f) { frequ_m = f; }
    
    string getTitle() { return string("pcyclint"); }
    
      
    //constructor
    ChargedParticles(Layout_t *playout, BCT bc) : 
            ParticleBase< ParticleSpatialLayout<T,Dim> >(playout),
            bcType_m(bc),
            numberOfClones(0)
        {
        
        // register the particle attributes
        addAttribute(q_m);
        addAttribute(P);
        addAttribute(Ef);
        addAttribute(mass_m);
        addAttribute(h);
        addAttribute(Ti);
        addAttribute(tact);
        addAttribute(bunchNo); 

        if (bcType_m == PPP) {
            setBCAllPeriodic();
        }
        else if (bcType_m == OOP) {
            setBCPeriodicInZOpenInXY();
        }
        else if (bcType_m == OOO) {
            setBCAllOpen();
        }
        else {
            INFOMSG("BC type not implemented use OOO instead" << endl);
            setBCAllOpen();
        }

        IpplRandom.SetSeed(static_cast<unsigned long>(234244113131719.0));

	q_m    = 0.0;
	mass_m = 0.0;
	h      = 0.0;
	Ti     = 0.0;
	tact   = 0.0;
	Ef = Vector_t(0.0);
	P  = Vector_t(0.0);
	R  = Vector_t(0.0);
    }

    inline const Mesh_t& getMesh() const { 
        return getLayout().getLayout().getMesh(); }
    inline Mesh_t& getMesh() { 
        return getLayout().getLayout().getMesh(); }
    inline const FieldLayout_t& getFieldLayout() const {
        return dynamic_cast<FieldLayout_t&>(getLayout().getLayout().getFieldLayout());
    }
    inline FieldLayout_t& getFieldLayout() {
        return dynamic_cast<FieldLayout_t&>(getLayout().getLayout().getFieldLayout());
    }
    
    inline void setCoupling(T c) {
        coupling_m = c;
    }
    
    inline T getCoupling() {
        return coupling_m;
    }

    inline void getBB(Vector_t &lowerCorner, Vector_t &upperCorner) {
        bounds(R, lowerCorner, upperCorner);
    }

    
    void clean() {
        scaleIn();
        double s=0.08;
        int count = 0;
        for (unsigned long int i=0; i<getLocalNum(); i++) {
            if ( (R[i](0)*R[i](0) > s*s) || (R[i](1)*R[i](1) > s*s) || (R[i](2)*R[i](2) > s*s)) {
                destroy(1,i);
                count ++;
            }
        }
        reduce(count,count,OpAddAssign());
        if (count>0) {
            update();
            INFOMSG("clean deleted " << count << "particles"<<endl);
        }
        scaleOut();
    }
    
    void boundp() {
      Inform msg2all("boundp ", INFORM_ALL_NODES);
      bounds(R, rmin_m, rmax_m);
      
      NDIndex<Dim> domain = getFieldLayout().getDomain(); 
      for(int i=0; i<Dim; i++)
	nr_m[i] = domain[i].length();
      
      // enlarge domaine 
      Vector_t dr = 1.05 * (rmax_m - rmin_m);
      
      for(int i=0; i<Dim; i++)
	hr_m[i] = (rmax_m[i] - rmin_m[i]) / (nr_m[i] - 1.0); 
      
      // rescale mesh 
      getMesh().set_meshSpacing( &(hr_m[0])  );
      getMesh().set_origin( rmin_m ); 
      // (re) initialize the fields
      rho_m.initialize(getMesh(), 
		       getFieldLayout(), 
		       GuardCellSizes<Dim>(1), 
		       bc_m);
      eg_m.initialize(getMesh(), 
		      getFieldLayout(), 
		      GuardCellSizes<Dim>(1), 
		      vbc_m);
      rho_m = 0.0;
      eg_m = Vector_t(0.0);
      update();

      if (getNumOfClones()==0) {
	locNum_m[0]=getLocalNum();
	totNum_m[0]=getTotalNum();
      }
      else {
	for (unsigned long int i=0; i<getLocalNum(); i++) 
	  locNum_m[bunchNo[i]]++;
	reduce(locNum_m,locNum_m+MAXCLONE,totNum_m,OpAddAssign());
      }
    }

    inline void scaleIn () { 
        for (unsigned long int i=0; i<getLocalNum(); i++)
            R[i] -= MeanR_m[bunchNo[i]];
        boundp();
    }
    

    inline void scaleOut (Vector_t center)  {
        /** special scaleOut for the first bunch in the machine
         */
        R += center;
    }
    
    inline void scaleOut ()  { 
        for (unsigned long int i=0; i<getLocalNum(); i++)
            R[i] += MeanR_m[bunchNo[i]];
    }
    

    inline void do_binaryRepart(int maxDev) {
        Inform msg("do_binaryRepart ");
        if(maxDev>0) {            
#ifdef USE_PBE
            pbe_start(2,"binaryRepart");
#endif  
            double idealNp     = getTotalNum()/Ippl::getNodes();
            double localDev    = 100.0*std::abs((idealNp-getLocalNum())/idealNp);
            double actMaxDev = 0.0;
            
            reduce(localDev, actMaxDev, OpMaxAssign());
            
            if (actMaxDev > maxDev) {
                msg << "do binaryRepart threshold is " << maxDev << " max inbalance reached is " << actMaxDev << endl;
                do_binaryRepart();
                localDev = 100.0*std::abs((idealNp-getLocalNum())/idealNp);
                actMaxDev = 0.0;
                reduce(localDev, actMaxDev, OpMaxAssign());
                msg << "after binaryRepart  max inbalance is now " << actMaxDev << endl;
            }
#ifdef USE_PBE
            pbe_stop(2);
#endif
        }
    }
    
    inline void do_binaryRepart() {
        scaleIn(); 
        BinaryRepartition(*this);
        scaleOut();
    }

    inline long getMinPartPerNode() {
        long  localmin = getLocalNum();
        long  globalMin = 0;
        reduce(localmin, globalMin, OpMinAssign());
        return globalMin;
    }
    



    inline void setBCAllOpen() {
        for (int i=0; i < 2*Dim; ++i) {
            bc_m[i] = new ZeroFace<T,Dim,Mesh_t,Center_t>(i);
            vbc_m[i] = new ZeroFace<Vector_t,Dim,Mesh_t,Center_t>(i);
            getBConds()[i] = ParticleNoBCond;
        }
    }
  
    inline void setBCAllPeriodic() {
        for (int i=0; i < 2*Dim; ++i) {
            bc_m[i] = new PeriodicFace<T,Dim,Mesh_t,Center_t>(i);
            vbc_m[i] = new PeriodicFace<Vector_t,Dim,Mesh_t,Center_t>(i);
            getBConds()[i] = ParticlePeriodicBCond;
        }
    }
  
    inline void setBCPeriodicInZOpenInXY() {
        bc_m[0] = new ZeroFace<T,Dim,Mesh_t,Center_t>(0);
        bc_m[1] = new ZeroFace<T,Dim,Mesh_t,Center_t>(1);
        bc_m[2] = new ZeroFace<T,Dim,Mesh_t,Center_t>(2);
        bc_m[3] = new ZeroFace<T,Dim,Mesh_t,Center_t>(3);
    
        vbc_m[0] = new ZeroFace<Vector_t,Dim,Mesh_t,Center_t>(0);
        vbc_m[1] = new ZeroFace<Vector_t,Dim,Mesh_t,Center_t>(1);
        vbc_m[2] = new ZeroFace<Vector_t,Dim,Mesh_t,Center_t>(2);
        vbc_m[3] = new ZeroFace<Vector_t,Dim,Mesh_t,Center_t>(3);
    
        getBConds()[0] = ParticleNoBCond;
        getBConds()[1] = ParticleNoBCond;
        getBConds()[2] = ParticleNoBCond;
        getBConds()[3] = ParticleNoBCond;
        getBConds()[4] = ParticlePeriodicBCond;
        getBConds()[5] = ParticlePeriodicBCond;

        if (Ippl::getNodes() == 1) {
            bc_m[4] = new PeriodicFace<T,Dim,Mesh_t,Center_t>(4);
            bc_m[5] = new PeriodicFace<T,Dim,Mesh_t,Center_t>(5);
            vbc_m[4] = new PeriodicFace<Vector_t,Dim,Mesh_t,Center_t>(4);
            vbc_m[5] = new PeriodicFace<Vector_t,Dim,Mesh_t,Center_t>(5);
    
        } else {
            bc_m[4] = new PeriodicFace<T,Dim,Mesh_t,Center_t>(4);
            bc_m[5] = new PeriodicFace<T,Dim,Mesh_t,Center_t>(5);
            vbc_m[4] = new PeriodicFace<Vector_t,Dim,Mesh_t,Center_t>(4);
            vbc_m[5] = new PeriodicFace<Vector_t,Dim,Mesh_t,Center_t>(5);
        }    
    }
     
    inline void writeHeader(ofstream &os, unsigned int idx, double ga, 
        unsigned long int N, Vector_t pmin, Vector_t pmax,  Vector_t cm, double t,
        double I, unsigned long N0,  Vector_t rcm, double r, double phi, unsigned int turn) {
        
        os << "# " << getTitle() << " format: x/m,y/m,z/m,px,y,pz,id (p == v/(beta*c)" << endl;
        os << "# s: "<< 0.0  << " rcm= " << rcm << " r= " << r << " phi= " << phi << " turn= " << turn << endl;
        os << "# struture lenght: 0 ... not used t= " << t << endl;
        os << "# Mesh spacing  " << getMeshSpacing() << endl;
        os << "# Origin not implemented " << endl;
        os << "# Max/Min position " <<  getRmax() <<  getRmin() << endl;
        os << "# Max/Min momenta " <<  pmax <<  pmin << endl;
        os << "# GridSize " << getGridSize() << " Centroid drive beam" << cm << endl;
        os << "# Protons  0 ... " <<  N << endl;
        os << "# Electrons 0" << endl;
        os << "# Initial particles in simulation = " << N0 << " I0= " << I << " /A" << endl;
        os << "# Energy drive beam (gamma) " << ga << endl;
        os << "# Data set  " << idx << endl;
    }
        
    inline void writePhaseSpace(string fn,unsigned int fnum, double rcm, double phi, unsigned int turn, unsigned int clone) {
        
        stringstream ff;
        ofstream of;
	Inform msg("writePhaseSpace: ", INFORM_ALL_NODES);
	
#ifdef USE_PBE
        pbe_start(10,"writePhaseSpace");
#endif                    
        if (clone==0)
            ff << fn  << setw(4) << setfill('0') << fnum << ".dat";
        else
            ff << fn + string("-clone-") << setw(4) << setfill('0') << fnum << ".dat";
        
        int tag = Ippl::Comm->next_tag(IPPL_APP_TAG4, IPPL_APP_CYCLE);
        

        double ga = getGamma(clone);
        double be = sqrt(1.0-(1.0/(ga*ga)));  
        
        double v0 = be*CLIGHT;
        double t = getTime(clone);
        
        Vector_t pmin, pmax;
        bounds(P,pmin, pmax);
        
        Vector_t cm = getCenterOfMass();
        scaleIn();
        
        double current = getCurrent();
        unsigned long N0 = getN0();
        
        if(isRoot()) {
            of.open(ff.str().c_str(),ios::out);
            of.precision(15);
            of.setf(ios::scientific, ios::floatfield);
            writeHeader(of,fnum,ga,getTotalnum(clone),pmax,pmin,cm,t,current,N0,MeanR_m[clone],rcm,phi,turn);
            
            unsigned int dataBlocks=0;
            unsigned int id=0;
            
            double x,y,z,px,py,pz;
            
            for (unsigned k = 0; k < getLocalNum(); k++) {
                if (bunchNo[k] == clone)
                    of <<  R[k](0) << "  " << R[k](1)   << "  " << R[k](2) << "  "
                       <<  P[k](0) << "  " << P[k](1)   << "  " << P[k](2) << "  " << ID[k] << endl;
            }
            
            int notReceived =  Ippl::getNodes() - 1; 
            while (notReceived > 0) {
                int node = COMM_ANY_NODE;
                Message* rmsg =  Ippl::Comm->receive_block(node, tag);
                if (rmsg == 0)
                    ERRORMSG("Could not receive from client nodes in main." << endl);
                notReceived--;
                rmsg->get(&dataBlocks);
                for (unsigned i=0; i < dataBlocks; i+=7) {
                    rmsg->get(&x);
                    rmsg->get(&px);
                    rmsg->get(&y);
                    rmsg->get(&py);
                    rmsg->get(&z);
                    rmsg->get(&pz);
                    rmsg->get(&id);
                    of << x  << " " << y << " " << z << " "   \
                       << px << " " << py << " " << pz << " " << id << endl;
                }
                delete rmsg;
            }
            of.close();
        }
        else {
            
            Message* smsg = new Message();
            unsigned dataBlocks = 7*getLocalnum(clone);
            smsg->put(dataBlocks);
            for (unsigned k=0; k<getLocalNum(); k++) {
                if (bunchNo[k] == clone) {
                    for (unsigned i=0; i < Dim; i++) {
                        smsg->put(R[k](i));
                        smsg->put(P[k](i));
                    }
                    smsg->put(ID[k]);
                }
            }
            bool res = Ippl::Comm->send(smsg, 0, tag);	
            if (! res) 
                ERRORMSG("Ippl::Comm->send(smsg, 0, tag) failed " << endl;);
        }
        IpplInfo::Comm->barrier();
        scaleOut();
#ifdef USE_PBE
        pbe_stop(10);
#endif                    
    }
    
    
    inline void writeStatistics(string fn, int turn, bool calcTune, double r, double phi, double dT) {
        ofstream ofs1,ofs2;
        Vector_t MeanR;
        Vector_t MeanP;
        Vector_t MeanR2;
        Vector_t MeanP2;
        Vector_t MeanRP;
        Vector_t Mean2RP;
        Vector_t rmsEmit;
        Vector_t TotalP;
        Vector_t TotalL;
        
        Vector_t rmsR;
        Vector_t rmsP;

#ifdef USE_PBE
            pbe_start(9,"writeStatistics");
#endif                    
            scaleIn();  
            
            T denN = 1.0 / getTotalNum();
            int Np  = getTotalNum();
        
            MeanR = sum(R) * denN;
            MeanP = sum(P) * denN;
            
            MeanR2 = sum(R*R) * denN;
            MeanP2 = sum(P*P) * denN;
            
            MeanRP = sum(R*P) * denN;
            Mean2RP = MeanRP * MeanRP;
            
            double ga = getGamma(0);
            

            // Calculate rms emittance and standard deviation of R and P: rmsR and rmsP
            for (int i=0;i<Dim;i++) {
                rmsEmit(i) = sqrt( (MeanR2*MeanP2)(i) - Mean2RP(i)  ) / ( getBeta()*CLIGHT );
                rmsEmit(2) /= ga*ga;
                rmsR(i) = sqrt( MeanR2(i) - (MeanR*MeanR)(i)  ); 
                rmsP(i) = sqrt( MeanP2(i) - (MeanP*MeanP)(i)  ); 
            }
            
            CenterOfMass_m = sum(mass_m*R) / sum(mass_m);
            
            //Calculate total momentum and angular momentum
            
            TotalP = sum(P);
            TotalL = sum(cross(R,P));
            
            // maximal radius
            Vector_t maxR = max(R);
            Vector_t maxE = max(Ef);
            
            Vector_t meanPrad = sum(MeanP - P)/MeanP;
            Vector_t rmsPrad = sum(rmsP - P)/rmsP;
            T actTime = getTime(0);
            
            
            if (isRoot()) {
                ofs1.precision(14);  
                ofs1.open(fn.c_str(),ios::app);
                ofs1.setf(ios_base::scientific); 
                ofs1 << MeanR[0]  << "  "  << MeanR[1]  << "  "   << MeanR[2]  << "  " // 1 2 3      [m]
                     << meanPrad[0] << "  "  << meanPrad[1] << "  "   << meanPrad[2] << "  " // 4 5 6     [rad]
                     << actTime  << "  "                                                     // 7          [s]
                     << rmsR[0]   << "  "  << rmsR[1]    << "  "   << rmsR[2]   << "  " // 8 9 10    [m]
                     << rmsPrad[0]  << "  "  << rmsPrad[1]   << "  "   << rmsPrad[2]  << "  " // 11 12 13    [rad]    
                     << maxR[0]     << "  "  << maxR[1]      << "  "   << maxR[2]     << "  " // 14 15 16    [m]
                     << maxE[0]     << "  "  << maxE[1]      << "  "   << maxE[2]     << "  " // 17 18 19     [V/m]
                     << rmsEmit[0] << "  " << rmsEmit[1] << "  "   << rmsEmit[2] << "  " // 20 21 22    
                     << TotalP[0] << "  "  << TotalP[1]  << "  "   << TotalP[2] << "  "  // 23 24 25 
                     << TotalL[0] << "  "  << TotalL[1]  << "  "   << TotalL[2] << "  "  // 26 27 28
                     << MeanR_m[0](0) << "  " << MeanR_m[0](1) << "  " << MeanR_m[0](2) << "  "                        // 29 30 31    [m]
                     << MeanP_m[0](0) << "  " << MeanP_m[0](1) << "  " << MeanP_m[0](2) << "  "                        // 32 33 34    [m]
                     << ga << "   "                 // 35
                     << r  << "   "                 // 36
                     << phi  << "   "               // 37
                     << endl;
                ofs1.close();
            }

            long locNum = getLocalNum();
            long minLocNum = 0;
            long maxLocNum = 0;
            
            reduce(locNum, minLocNum, OpMinAssign());
            reduce(locNum, maxLocNum, OpMaxAssign());

            if (isRoot()) {
                ofs1.precision(14);  
                ofs1.open((fn+string("-comp")).c_str(),ios::app);
                ofs1.setf(ios_base::scientific); 
                
                ofs1 << dT << "  " <<  minLocNum << "  " << maxLocNum << endl;
                
                ofs1.close();
            }
            
            IpplInfo::Comm->barrier();
            scaleOut();
#ifdef USE_PBE
            pbe_stop(9);
#endif                    
    }
  
};
#endif
