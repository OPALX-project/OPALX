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
 * ParticleAmrLayout - particle layout based on BoxLib AMR framework.
 *
 * This is a specialized version of ParticleLayout, which places particles
 * on processors based on their spatial location relative to multilevel grid.
 * The grids are defined by AMR framework (BoxLib) and can contain multiple levels.
 * Layout uses specialized AmrParticleBase class and BoxLibs ParGDB to determine to 
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

#include <ParGDB.H>
#include <REAL.H>
#include <IntVect.H>
#include <Array.H>
#include <Utility.H>
#include <Geometry.H>
#include <VisMF.H>
#include <Particles.H>
#include <RealBox.H>

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

private:

  ParGDBBase* m_gdb;
  ParGDB m_gdb_object;

  //
  // Checks/sets a particles location on levels lev_min and higher.
  // Returns false if the particle does not exist on that level.
  //
  static bool Where (AmrParticleBase< ParticleAmrLayout<T,Dim> >& prt, 
		     const unsigned int ip, const ParGDBBase* gdb, 
		     int lev_min = 0, int finest_level = -1);
  //
  // Checks/sets whether the particle has crossed a periodic boundary in such a way
  // that it is on levels lev_min and higher.
  //
  static bool PeriodicWhere (AmrParticleBase< ParticleAmrLayout<T,Dim> >& prt, 
			     const unsigned int ip, const ParGDBBase* gdb, 
			     int lev_min = 0, int finest_level = -1);
  //
  // Checks/sets whether a particle is within its grid (including grow cells).
  //
  static bool RestrictedWhere (AmrParticleBase< ParticleAmrLayout<T,Dim> >& p, 
			       const unsigned int ip, const ParGDBBase* gdb, int ngrow);

  static bool PeriodicShift (SingleParticlePos_t R, const ParGDBBase* gdb);

public:

  ParticleAmrLayout() : 
    m_gdb(nullptr) {}

  ParticleAmrLayout(ParGDBBase* gdb) : 
    m_gdb(gdb) {}

  ParticleAmrLayout(const Geometry &geom, 
		    const DistributionMapping &dmap, 
		    const BoxArray &ba) :
    m_gdb_object(geom, dmap, ba)
  {
      m_gdb = &m_gdb_object;
  }

  ParticleAmrLayout(const Array<Geometry>            & geom, 
		    const Array<DistributionMapping> & dmap,
		    const Array<BoxArray>            & ba,
		    const Array<int>                 & rr):
    m_gdb_object(geom,dmap,ba,rr)
  {
    m_gdb = & m_gdb_object;
  }

  void Define (ParGDBBase* gdb) 
  {
    m_gdb = gdb;
  }

  void Define (const Geometry &geom, 
	       const DistributionMapping &dmap, 
	       const BoxArray &ba) 
  {
    m_gdb_object = ParGDB(geom, dmap, ba);
    m_gdb = &m_gdb_object;
  }

  void Define (const Array<Geometry>            & geom, 
	       const Array<DistributionMapping> & dmap,
	       const Array<BoxArray>            & ba,
	       const Array<int>                 & rr)
  {
    m_gdb_object = ParGDB(geom, dmap, ba, rr);
    m_gdb = &m_gdb_object;
  }

  void SetParticleBoxArray(int lev, const BoxArray& new_ba)
  {
    m_gdb->SetParticleBoxArray(lev, new_ba);
  }

  void SetParticleDistributionMap(int lev, const DistributionMapping& new_dmap)
  {
    m_gdb->SetParticleDistributionMap(lev, new_dmap);
  }

  const BoxArray& ParticleBoxArray (int lev) const 
  { 
    return m_gdb->ParticleBoxArray(lev); 
  }

  const DistributionMapping& ParticleDistributionMap (int lev) const 
  { 
    return m_gdb->ParticleDistributionMap(lev); 
  }

  const ParGDBBase* GetParGDB () const
  { 
    return m_gdb; 
  }

  //get the cell of the particle
  static IntVect Index (AmrParticleBase< ParticleAmrLayout<T,Dim> >& p, 
			const unsigned int ip,
			const Geometry& geom);

  //get the cell of the particle
  static IntVect Index (SingleParticlePos_t &R, const Geometry& geom);

  void Redistribute(AmrParticleBase< ParticleAmrLayout<T,Dim> >& PData,
		    bool where_already_called = false,
		    bool full_where = false,
		    int lev_min = 0,
		    int nGrow = 0)
  {

    unsigned N = Ippl::getNodes();
    unsigned myN = Ippl::myNode();

    int theEffectiveFinestLevel = m_gdb->finestLevel();
    while (!m_gdb->LevelDefined(theEffectiveFinestLevel))
      theEffectiveFinestLevel--;

    //loop trough the particles and and assigne the grid and level where each particle belongs
    size_t LocalNum = PData.getLocalNum();

    std::multimap<unsigned, unsigned> p2n; //node ID, particle 

    int *msgsend = new int[N];
    std::fill(msgsend, msgsend+N, 0);
    int *msgrecv = new int[N];
    std::fill(msgrecv, msgrecv+N, 0);

    unsigned sent = 0;
    unsigned particlesLeft = LocalNum;
  
    for (unsigned int ip=0; ip < LocalNum; ++ip) {
      bool particleLeftDomain = false;
       //check to which level and grid the particle belongs to
      if (!where_already_called) {
	if (!Where(PData, ip, m_gdb, lev_min, theEffectiveFinestLevel)) {
	  if (full_where) {
	    if (!PeriodicWhere(PData, ip, m_gdb, lev_min, theEffectiveFinestLevel)) {
	      if (lev_min != 0) {
		if (!RestrictedWhere(PData, ip, m_gdb, nGrow))
		  BoxLib::Abort("RestrictedWhere failed in redistribute");
	      } else {
		//the particle has left the domain - invalidate it.
		particleLeftDomain = true;
		PData.destroy(1, ip);
	      }
	    }
	  } else {
	    std::cout << "Bad Particle: " << PData.R[ip] << std::endl;
	  }
	}
      }

      if (!particleLeftDomain) {
	//get the node the particle belongs to
	unsigned int who = m_gdb->ParticleDistributionMap(PData.m_lev[ip])[PData.m_grid[ip]];
	if (who == myN) {
	  //TODO: if grid and level of the particle has changed move it to the right place
	  //to keep the local particles sorted
	} else {
	  msgsend[who] = 1;
	  p2n.insert(std::pair<unsigned, unsigned>(who, ip));
	  sent++;
	  particlesLeft--;
	}
      }
    }	      

    //reduce message count so every node knows how many messages to receive
    MPI_Allreduce(msgsend, msgrecv, N, MPI_INT, MPI_SUM, Ippl::getComm());

    int tag = Ippl::Comm->next_tag(P_SPATIAL_TRANSFER_TAG,P_LAYOUT_CYCLE);

    typename std::multimap<unsigned, unsigned>::iterator i = p2n.begin();

    Format *format = PData.getFormat();

    std::vector<MPI_Request> requests;
    std::vector<MsgBuffer*> buffers;

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

      
      MPI_Request request = Ippl::Comm->raw_isend( msgbuf->getBuffer(), 
						   msgbuf->getSize(), 
						   cur_destination, tag);

      //remember request and buffer so we can delete them later
      requests.push_back(request);
      buffers.push_back(msgbuf);
      
    }

    //destroy the particles that are sent to other domains
    LocalNum -= PData.getDestroyNum();  // update local num
    PData.performDestroy();

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
	LocalNum += PData.getSingleMessage(*msg);
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

  void update(AmrParticleBase< ParticleAmrLayout<T,Dim> >& PData, 
	      const ParticleAttrib<char>* canSwap=0)
  {
    unsigned N = Ippl::getNodes();
    unsigned myN = Ippl::myNode();
    
    size_t LocalNum = PData.getLocalNum();
    size_t DestroyNum = PData.getDestroyNum();
    size_t TotalNum;
    int node;

    //delete particles in destroy list, update local num
    //PData.performDestroy();
    //LocalNum -= DestroyNum;

    Redistribute(PData, false, false, 0, 0);

    std::cout << "AmrBase update" << std::endl;

  }

  
  void update(IpplParticleBase< ParticleAmrLayout<T,Dim> >& PData, 
	      const ParticleAttrib<char>* canSwap=0)
  {
    std::cout << "IpplBase update" << std::endl;
    //TODO: exit since we need AmrParticleBase with grids and levels for particles for this layout
    //if IpplParticleBase is used something went wrong
  }
  

};

#include "ParticleAmrLayout.hpp"

#endif
