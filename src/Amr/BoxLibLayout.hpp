#ifndef BoxLibLayout_HPP
#define BoxLibLayout_HPP

#include "BoxLibLayout.h"

#include "Message/Formatter.h"

#include <map>
#include <vector>

// #include <mpi.h>

// // Function from BoxLib adjusted to work with Ippl AmrParticleBase class
// // redistribute the particles using BoxLibs ParGDB class to determine where particle should go
template<class T, unsigned Dim>
void BoxLibLayout<T, Dim::update(ParticleBase_t& PData,
                                 int lev_min,
                                 const ParticleAttrib<char> canSwap)
{
    // Input parameters of ParticleContainer::Redistribute of BoxLib
    bool where_already_called = false;
    bool full_where = false;
    int nGrow = 0;
    //
    

    unsigned N = Ippl::getNodes();
    unsigned myN = Ippl::myNode();

    int theEffectiveFinestLevel = this->finestLevel();
    while (!this->LevelDefined(theEffectiveFinestLevel))
        theEffectiveFinestLevel--;
    
    //loop trough the particles and assigne the grid and level where each particle belongs
    size_t LocalNum = PData.getLocalNum();

    std::multimap<unsigned, unsigned> p2n; //node ID, particle 

    int *msgsend = new int[N];
    std::fill(msgsend, msgsend+N, 0);
    int *msgrecv = new int[N];
    std::fill(msgrecv, msgrecv+N, 0);

    unsigned sent = 0;
    unsigned particlesLeft = LocalNum;
  
    //loop trough particles and assign grid and level to each particle
    //if particle doesn't belong to this process save the index of the particle to be sent
    for (unsigned int ip=0; ip < LocalNum; ++ip) {
        bool particleLeftDomain = false;
        //check to which level and grid the particle belongs to
        if (!where_already_called) {
            if (!Where(PData, ip, lev_min, theEffectiveFinestLevel)) {
                if (full_where) {
                    if (!PeriodicWhere(PData, ip, lev_min, theEffectiveFinestLevel)) {
                        if (lev_min != 0) {
                            if (!RestrictedWhere(PData, ip, nGrow))
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
            unsigned int who = this->ParticleDistributionMap(PData.level[ip])[PData.grid[ip]];
            if (who != myN) {
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


template <class T, unsigned Dim, class PT, class CacheData>
void BoxLibLayout::scatter(AmrField_t& f,
                           const ParticleAttrib<Vektor<PT, Dim>& pp,
                           ParticleAttrib<CacheData>& cache,
                           int lbase, int lfine) const
{
    
    
    
}


template <class T, unsigned Dim, class PT, class CacheData>
void BoxLibLayout::gather(AmrField_t& f,
                          const ParticleAttrib<Vektor<PT, Dim>& pp,
                          ParticleAttrib<CacheData>& cache,
                          int lbase, int lfine) const
{
    
    
    
}


// Function from BoxLib adjusted to work with Ippl AmrParticleBase class
//get the cell where particle is located - uses AmrParticleBase object and particle id
template <class T, unsigned Dim>
IntVect BoxLibLayout<T, Dim>::Index (ParticleBase_t>& p,
                                     const unsigned int ip,
                                     const Geometry&     geom)
{
    IntVect iv;

    D_TERM(iv[0]=floor((p.R[ip][0]-geom.ProbLo(0))/geom.CellSize(0));,
           iv[1]=floor((p.R[ip][1]-geom.ProbLo(1))/geom.CellSize(1));,
           iv[2]=floor((p.R[ip][2]-geom.ProbLo(2))/geom.CellSize(2)););

    iv += geom.Domain().smallEnd();

    return iv;
}


//get the cell where particle is located - uses the particle position vector R
template <class T, unsigned Dim>
IntVect BoxLibLayout<T, Dim>::Index (SingleParticlePos_t &R,
                                     const Geometry&     geom)
{
    IntVect iv;

    D_TERM(iv[0]=floor((R[0]-geom.ProbLo(0))/geom.CellSize(0));,
           iv[1]=floor((R[1]-geom.ProbLo(1))/geom.CellSize(1));,
           iv[2]=floor((R[2]-geom.ProbLo(2))/geom.CellSize(2)););

    iv += geom.Domain().smallEnd();

    return iv;
}


//sets the grid and level where particle belongs - returns flase if prticle is outside the domain
template <class T, unsigned Dim>
bool BoxLibLayout<T, Dim>::Where_m(ParticleBase_t& p,
                                   const unsigned int ip,
                                   int lev_min,
                                   int finest_level)
{
    BL_ASSERT(this != 0);

    if (finest_level == -1)
        finest_level = this->finestLevel();

    BL_ASSERT(finest_level <= this->finestLevel());

    std::vector< std::pair<int,Box> > isects;

    for (unsigned int lev = (unsigned)finest_level; lev >= (unsigned)lev_min; lev--)
    {
      const IntVect& iv = Index(p, ip, this->Geom(lev));

        if (lev == p.level[ip]) { 
            // We may take a shortcut because the fact that we are here means 
            // this particle does not belong to any finer grids.
            const BoxArray& ba = this->ParticleBoxArray(p.level[ip]);
            if (0 <= p.grid[ip] && p.grid[ip] < ba.size() && 
                ba[p.grid[ip]].contains(iv)) 
            {
                return true;
            }
        }

        this->ParticleBoxArray(lev).intersections(Box(iv,iv),isects,true,0);

        if (!isects.empty())
        {
            p.level[ip]  = lev;
            p.grid[ip] = isects[0].first;

            return true;
        }
    }
    return false;
}


//Function from BoxLib adjusted to work with Ippl AmrParticleBase class
//Checks/sets whether the particle has crossed a periodic boundary in such a way
//that it is on levels lev_min and higher.
template <class T, unsigned Dim>
bool BoxLibLayout<T, Dim>::PeriodicWhere_m(ParticleBase_t& p,
                                           const unsigned int ip,
                                           int lev_min,
                                           int finest_level)
{
    BL_ASSERT(this != 0);

    if (!this->Geom(0).isAnyPeriodic()) return false;

    if (finest_level == -1)
        finest_level = this->finestLevel();

    BL_ASSERT(finest_level <= this->finestLevel());
    //
    // Create a copy "dummy" particle to check for periodic outs.
    //
    SingleParticlePos_t R = p.R[ip];

    if (PeriodicShift(R))
    {
        std::vector< std::pair<int,Box> > isects;

        for (int lev = finest_level; lev >= lev_min; lev--)
        {
            const IntVect& iv = Index(R, this->Geom(lev));

            this->ParticleBoxArray(lev).intersections(Box(iv,iv),isects,true,0);

            if (!isects.empty())
            {
                D_TERM(p.R[ip][0] = R[0];,
                       p.R[ip][1] = R[1];,
                       p.R[ip][2] = R[2];);

                p.level[ip]  = lev;
                p.grid[ip] = isects[0].first;

                return true;
            }
        }
    }

    return false;
}


// Function from BoxLib adjusted to work with Ippl AmrParticleBase class
// Checks/sets whether a particle is within its grid (including grow cells).
template <class T, unsigned Dim>
bool BoxLibLayout<T, Dim>::RestrictedWhere_m(ParticleBase_t& p,
                                             const unsigned int ip,
                                             int ngrow)
{
    BL_ASSERT(this != 0);

    const IntVect& iv = Index(p,ip,this->Geom(p.level[ip]));

    if (BoxLib::grow(this->ParticleBoxArray(p.level[ip])[p.grid[ip]], ngrow).contains(iv))
    {
        return true;
    }

    return false;
}


// Function from BoxLib adjusted to work with Ippl AmrParticleBase class
// Returns true if the particle was shifted.
template <class T, unsigned Dim>
bool BoxLibLayout<T, Dim>::PeriodicShift_m(SingleParticlePos_t R)
{
    //
    // This routine should only be called when Where() returns false.
    //
    BL_ASSERT(this != 0);
    //
    // We'll use level 0 stuff since ProbLo/ProbHi are the same for every level.
    //
    const Geometry& geom    = this->Geom(0);
    const Box&      dmn     = geom.Domain();
    const IntVect&  iv      = Index(R,this->Geom(0));
    bool            shifted = false;  

    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        if (!geom.isPeriodic(i)) continue;

        if (iv[i] > dmn.bigEnd(i))
        {
            if (R[i] == geom.ProbHi(i))
                //
                // Don't let particles lie exactly on the domain face.
                // Force the particle to be outside the domain so the
                // periodic shift will bring it back inside.
                //
                R[i] += .125*geom.CellSize(i);

            R[i] -= geom.ProbLength(i);

            if (R[i] <= geom.ProbLo(i))
                //
                // This can happen due to precision issues.
                //
                R[i] += .125*geom.CellSize(i);

            BL_ASSERT(R[i] >= geom.ProbLo(i));

            shifted = true;
        }
        else if (iv[i] < dmn.smallEnd(i))
        {
            if (R[i] == geom.ProbLo(i))
                //
                // Don't let particles lie exactly on the domain face.
                // Force the particle to be outside the domain so the
                // periodic shift will bring it back inside.
                //
                R[i] -= .125*geom.CellSize(i);

            R[i] += geom.ProbLength(i);

            if (R[i] >= geom.ProbHi(i))
                //
                // This can happen due to precision issues.
                //
                R[i] -= .125*geom.CellSize(i);

            BL_ASSERT(R[i] <= geom.ProbHi(i));

            shifted = true;
        }
    }
    //
    // The particle may still be outside the domain in the case
    // where we aren't periodic on the face out which it travelled.
    //
    return shifted;
}

#endif
