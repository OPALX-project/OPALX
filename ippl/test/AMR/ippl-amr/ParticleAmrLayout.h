// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 *
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/
#ifndef PARTICLE_AMR_LAYOUT_H
#define PARTICLE_AMR_LAYOUT_H

/*
 * ParticleAmrLayout - particle layout based on AMReX AMR framework.
 *
 * This is a specialized version of ParticleLayout, which places particles
 * on processors based on their spatial location relative to multilevel grid.
 * The grids are defined by AMR framework (AMReX) and can contain multiple levels.
 * Layout uses specialized AmrParticleBase class and AMReXs ParGDB to determine to 
 * which grid and level the particle belongs to.
 */

#include "Particle/ParticleLayout.h"
#include "AmrParticleBase.h"

#include "Region/RegionLayout.h"
#include "Message/Message.h"
#include "FieldLayout/FieldLayoutUser.h"
#include <cstddef>

#include <vector>
#include <iostream>
#include <map>

#include "Message/Formatter.h"
#include <mpi.h>

#include <AMReX_ParGDB.H>
#include <AMReX_REAL.H>
#include <AMReX_IntVect.H>
#include <AMReX_Array.H>
#include <AMReX_Utility.H>
#include <AMReX_Geometry.H>
#include <AMReX_VisMF.H>
#include <AMReX_Particles.H>
#include <AMReX_RealBox.H>

using namespace amrex;

template <class T, unsigned Dim>
class ParticleAmrLayout : public ParticleLayout<T, Dim>
{

public:
    // pair iterator definition ... this layout does not allow for pairlists
    typedef int pair_t;
    typedef pair_t* pair_iterator;
    typedef typename ParticleLayout<T, Dim>::SingleParticlePos_t
    SingleParticlePos_t;
    typedef typename ParticleLayout<T, Dim>::Index_t Index_t;
  
    // type of attributes this layout should use for position and ID
    typedef ParticleAttrib<SingleParticlePos_t> ParticlePos_t;
    typedef ParticleAttrib<Index_t>             ParticleIndex_t;
    
    static bool do_tiling;
    static IntVect tile_size;
    
private:
  
    //ParGDBBase class from AMReX is used to determine grid, level 
    //and node where each particle belongs
    ParGDBBase* m_gdb;
    ParGDB m_gdb_object;

    // Function from AMReX adjusted to work with Ippl AmrParticleBase class
    // Checks/sets a particles location on levels lev_min and higher.
    // Returns false if the particle does not exist on that level.
    bool Where (AmrParticleBase< ParticleAmrLayout<T,Dim> >& p,
                const unsigned int ip, 
                int lev_min = 0, int lev_max = -1, int nGrow=0) const;

    // Function from AMReX adjusted to work with Ippl AmrParticleBase class
    // Checks/sets whether the particle has crossed a periodic boundary in such a way
    // that it is on levels lev_min and higher.
    bool EnforcePeriodicWhere (AmrParticleBase< ParticleAmrLayout<T,Dim> >& prt,
                               const unsigned int ip,
                               int lev_min = 0, int lev_max = -1) const;

    // Function from AMReX adjusted to work with Ippl AmrParticleBase class
    // Returns true if the particle was shifted.
    bool PeriodicShift (SingleParticlePos_t R) const;

public:

    //empty constructor: Define needs to be called in order to initialize ParGDB object 
    //for redistribute to function properly
    ParticleAmrLayout() : 
        m_gdb(nullptr) {}

    //constructor: takes ParGDBBase as argument
    ParticleAmrLayout(ParGDBBase* gdb) : 
        m_gdb(gdb) {}

    //constructor: takes AMReX Geometry, DistributionMapping and BoxArray objects as arguments and
    //defines a new ParGDBBase object
    ParticleAmrLayout(const Geometry &geom, 
                      const DistributionMapping &dmap, 
                      const BoxArray &ba) :
        m_gdb_object(geom, dmap, ba)
    {
        m_gdb = &m_gdb_object;
    }

    //constructor: takes AMReX Geometry, DistributionMapping, BoxArray and an array of refinements
    //at each level and constructs a new ParGDBBase object
    ParticleAmrLayout(const Array<Geometry>            & geom, 
                      const Array<DistributionMapping> & dmap,
                      const Array<BoxArray>            & ba,
                      const Array<int>                 & rr):
        m_gdb_object(geom,dmap,ba,rr)
    {
        m_gdb = & m_gdb_object;
    }

    //define the ParGDBBase object
    void Define (ParGDBBase* gdb) 
    {
        m_gdb = gdb;
    }
  
    //create new ParGDBBase using Geometry, DistributionMapping and BoxArray
    void Define (const Geometry &geom, 
                 const DistributionMapping &dmap, 
                 const BoxArray &ba) 
    {
        m_gdb_object = ParGDB(geom, dmap, ba);
        m_gdb = &m_gdb_object;
    }

    //create new ParGDBBase using Geometry, DistributionMapping, BoxArray and 
    //array with refinements at each level
    void Define (const Array<Geometry>            & geom, 
                 const Array<DistributionMapping> & dmap,
                 const Array<BoxArray>            & ba,
                 const Array<int>                 & rr)
    {
        m_gdb_object = ParGDB(geom, dmap, ba, rr);
        m_gdb = &m_gdb_object;
    }

    //set a new BoxArray for the ParGDBBase at level lev
    void SetParticleBoxArray(int lev, const BoxArray& new_ba)
    {
        m_gdb->SetParticleBoxArray(lev, new_ba);
    }

    //set a new DistributionMapping for ParGDBBase at level lev
    void SetParticleDistributionMap(int lev, const DistributionMapping& new_dmap)
    {
        m_gdb->SetParticleDistributionMap(lev, new_dmap);
    }

    //get the BoxArray at level lev
    const BoxArray& ParticleBoxArray (int lev) const 
    { 
        return m_gdb->ParticleBoxArray(lev); 
    }

    //get the Distribution mapping at level lev
    const DistributionMapping& ParticleDistributionMap (int lev) const 
    { 
        return m_gdb->ParticleDistributionMap(lev); 
    }
    
    const Geometry& Geom (int lev) const { return m_gdb->Geom(lev); }
    
    int finestLevel () const { return m_gdb->finestLevel(); }
    int maxLevel ()    const { return m_gdb->maxLevel(); }

    //get the PartGDBBase object
    const ParGDBBase* GetParGDB () const
    { 
        return m_gdb; 
    }

    // Function from AMReX adjusted to work with Ippl AmrParticleBase class
    //get the cell of the particle
    IntVect Index (AmrParticleBase< ParticleAmrLayout<T,Dim> >& p,
                   const unsigned int ip, int leve) const;

    // Function from AMReX adjusted to work with Ippl AmrParticleBase class
    //get the cell of the particle
    IntVect Index (SingleParticlePos_t &R, int lev) const;

    // Function from AMReX adjusted to work with Ippl AmrParticleBase class
    //redistribute the particles using AMReXs ParGDB class to determine where particle should go
    void Redistribute(AmrParticleBase< ParticleAmrLayout<T,Dim> >& PData,
                      int lev_min = 0, int lev_max = -1, int nGrow = 0)
    {
        unsigned N = Ippl::getNodes();
        unsigned myN = Ippl::myNode();

        int theEffectiveFinestLevel = m_gdb->finestLevel();
        while (!m_gdb->LevelDefined(theEffectiveFinestLevel))
            theEffectiveFinestLevel--;
        
        if (lev_max == -1)
            lev_max = theEffectiveFinestLevel;

        //loop trough the particles and assigne the grid and level where each particle belongs
        size_t LocalNum = PData.getLocalNum();
        
        typedef typename AmrParticleBase< ParticleAmrLayout<T,Dim> >::LevelNumCounter_t LevelNumCounter_t;
        LevelNumCounter_t& LocalNumPerLevel = PData.getLocalNumPerLevel();
        
        std::multimap<unsigned, unsigned> p2n; //node ID, particle 

        int *msgsend = new int[N];
        std::fill(msgsend, msgsend+N, 0);
        int *msgrecv = new int[N];
        std::fill(msgrecv, msgrecv+N, 0);

        unsigned sent = 0;
        unsigned particlesLeft = LocalNum;
        
        int lBegin = LocalNumPerLevel.begin(lev_min);
        int lEnd   = LocalNumPerLevel.end(lev_max);
        
        //loop trough particles and assign grid and level to each particle
        //if particle doesn't belong to this process save the index of the particle to be sent
        for (int ip=lBegin; ip < lEnd; ++ip) {
            
            bool particleLeftDomain = false;
            
            size_t lold = PData.m_lev[ip];
            
            //check to which level and grid the particle belongs to
            locateParticle(PData, ip, lev_min, lev_max, nGrow, particleLeftDomain);
            
            if ( !particleLeftDomain ) {
                // The owner of the particle is the CPU owning the finest grid
                // in state data that contains the particle.
                const unsigned int who = m_gdb->ParticleDistributionMap(PData.m_lev[ip])[PData.m_grid[ip]];
                
                if (who != myN) {
                    msgsend[who] = 1;
                    p2n.insert(std::pair<unsigned, unsigned>(who, ip));
                    sent++;
                    particlesLeft--;
                    
                    --LocalNumPerLevel[lold];
                } else {
                    // if we own it it may have moved to another level
                    --LocalNumPerLevel[lold];
                    ++LocalNumPerLevel[PData.m_lev[ip]];
                }
            } else {
                // a particle left the domain
                std::cerr << "\033[01;31mParticle left domain.\033[0m" << std::endl;
                --LocalNumPerLevel[ lold ];
            }
        }

        //reduce message count so every node knows how many messages to receive
        MPI_Allreduce(msgsend, msgrecv, N, MPI_INT, MPI_SUM, Ippl::getComm());

        int tag = Ippl::Comm->next_tag(P_SPATIAL_TRANSFER_TAG,P_LAYOUT_CYCLE);

        typename std::multimap<unsigned, unsigned>::iterator i = p2n.begin();

        Format *format = PData.getFormat();

        std::vector<MPI_Request> requests;
        std::vector<MsgBuffer*> buffers;

        //create a message and send particles to nodes they belong to
        while (i!=p2n.end())
        {
            unsigned cur_destination = i->first;
      
            MsgBuffer *msgbuf = new MsgBuffer(format, p2n.count(i->first));

            for (; i!=p2n.end() && i->first == cur_destination; ++i)
            {
                Message msg;
                PData.putMessage(msg, i->second);
                PData.destroy(1, i->second);
                msgbuf->add(&msg);
            }

      
            MPI_Request request = Ippl::Comm->raw_isend(msgbuf->getBuffer(), 
                                                        msgbuf->getSize(), 
                                                        cur_destination, tag);

            //remember request and buffer so we can delete them later
            requests.push_back(request);
            buffers.push_back(msgbuf);
      
        }

        //destroy the particles that are sent to other domains
        LocalNum -= PData.getDestroyNum();  // update local num
        PData.performDestroy();
        
        for (int lev = lev_min; lev < lev_max; ++lev) {
            if ( LocalNumPerLevel[lev] < 0 )
                amrex::Abort("ParticleAmrLayout::Redistribute(): Negative particle level count.");
        }
        
        //receive new particles
        for (int k = 0; k<msgrecv[myN]; ++k)
        {
            int node = Communicate::COMM_ANY_NODE;
            char *buffer = 0;
            int bufsize = Ippl::Comm->raw_probe_receive(buffer, node, tag);
            MsgBuffer recvbuf(format, buffer, bufsize);
      
            Message *msg = recvbuf.get();
            while (msg != 0)
                {
                    /* pBeginIdx is the start index of the new particle data
                     * pEndIdx is the end index of the new particle data
                     */
                    size_t pBeginIdx = LocalNum;
                    
                    LocalNum += PData.getSingleMessage(*msg);
                    
                    size_t pEndIdx = LocalNum;
                    
                    for (size_t i = pBeginIdx; i < pEndIdx; ++i)
                        ++LocalNumPerLevel[ PData.m_lev[i] ];
                    
                    delete msg;
                    msg = recvbuf.get();
                }  
        }
        
        //wait for communication to finish and clean up buffers
        MPI_Waitall(requests.size(), &(requests[0]), MPI_STATUSES_IGNORE);
        for (unsigned int j = 0; j<buffers.size(); ++j)
        {
            delete buffers[j];
        }
    
        delete[] msgsend;
        delete[] msgrecv;
        delete format;

        // there is extra work to do if there are multipple nodes, to distribute
        // the particle layout data to all nodes
        //TODO: do we need info on how many particles are on each node?

        //save how many total particles we have
        size_t TotalNum = 0;
        MPI_Allreduce(&LocalNum, &TotalNum, 1, MPI_INT, MPI_SUM, Ippl::getComm());

        // update our particle number counts
        PData.setTotalNum(TotalNum);	// set the total atom count
        PData.setLocalNum(LocalNum);	// set the number of local atoms
    
    }

    //update the location and indices of all atoms in the given AmrParticleBase object.
    //uses the Redistribute function to swap particles among processes if needed
    //handles create and destroy requests. When complete all nodes have correct layout information
    void update(AmrParticleBase< ParticleAmrLayout<T,Dim> >& PData,
                const ParticleAttrib<char>* canSwap=0,
                int lbase = 0, int lfine = -1)
    {
        unsigned N = Ippl::getNodes();
        unsigned myN = Ippl::myNode();
    
        size_t LocalNum = PData.getLocalNum();
        size_t DestroyNum = PData.getDestroyNum();
        size_t TotalNum;
        int node;

        Redistribute(PData, lbase, lfine);

    }

  
    void update(IpplParticleBase< ParticleAmrLayout<T,Dim> >& PData, 
                const ParticleAttrib<char>* canSwap=0)
    {
        std::cout << "IpplBase update" << std::endl;
        //TODO: exit since we need AmrParticleBase with grids and levels for particles for this layout
        //if IpplParticleBase is used something went wrong
    }
    
    
private:
    int getTileIndex(const IntVect& iv, const Box& box, Box& tbx);
    
    void locateParticle(AmrParticleBase< ParticleAmrLayout<T,Dim> >& p, 
                        const unsigned int ip,
                        int lev_min, int lev_max, int nGrow,
                        bool &particleLeftDomain) const
    {
        bool outside = D_TERM( p.R[ip](0) <  Geometry::ProbLo(0)
                            || p.R[ip](0) >= Geometry::ProbHi(0),
                            || p.R[ip](1) <  Geometry::ProbLo(1)
                            || p.R[ip](1) >= Geometry::ProbHi(1),
                            || p.R[ip](2) <  Geometry::ProbLo(2)
                            || p.R[ip](2) >= Geometry::ProbHi(2));
        
        bool success;
        if (outside)
        {
            // Note that EnforcePeriodicWhere may shift the particle if it is successful.
            success = EnforcePeriodicWhere(p, ip, lev_min, lev_max);
            if (!success && lev_min == 0)
            {
                // The particle has left the domain; invalidate it.
                particleLeftDomain = true;
                p.destroy(1, ip);
                success = true;
            }
        }
        else
        {
            success = Where(p, ip, lev_min, lev_max);
        }
    
        if (!success)
        {
            success = (nGrow > 0) && Where(p, ip, lev_min, lev_min, nGrow);
        }
    
        if (!success)
        {
            amrex::Abort("ParticleContainer<NR, NI, NA>::locateParticle(): invalid particle.");
        }
    }

};

#include "ParticleAmrLayout.hpp"

#endif
