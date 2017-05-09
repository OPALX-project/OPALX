#ifndef BoxLibLayout_HPP
#define BoxLibLayout_HPP

#include "BoxLibLayout.h"

#include "Message/Formatter.h"

template<class T, unsigned Dim>
BoxLibLayout<T, Dim>::BoxLibLayout()
{
    /* figure out how big the maximum grid can
     * be in order that all processes have some data at the beginning
     */
    int nProcs = Ippl::getNodes();
    int nGridPoints = 16 * nProcs;      //FIXME
    int maxGridSize = 16;
    double lower = -0.1; // m
    double upper =  0.1; // m
    
    this->initDefaultBox(nGridPoints, maxGridSize, lower, upper);
}


template<class T, unsigned Dim>
BoxLibLayout<T, Dim>::BoxLibLayout(const Geometry &geom,
                                   const DistributionMapping &dmap,
                                   const BoxArray &ba)
    : ParGDB(geom, dmap, ba)
{ }


template<class T, unsigned Dim>
BoxLibLayout<T, Dim>::BoxLibLayout(const Array<Geometry> &geom,
                                   const Array<DistributionMapping> &dmap,
                                   const Array<BoxArray> &ba,
                                   const Array<int> &rr)
    : ParGDB(geom, dmap, ba, rr)
{ }


template<class T, unsigned Dim>
void BoxLibLayout<T, Dim>::update(IpplParticleBase< BoxLibLayout<T,Dim> >& PData,
                                  const ParticleAttrib<char>* canSwap)
{
    std::cout << "IpplBase update" << std::endl;
    //TODO: exit since we need AmrParticleBase with grids and levels for particles for this layout
    //if IpplParticleBase is used something went wrong
}

// // Function from BoxLib adjusted to work with Ippl AmrParticleBase class
// // redistribute the particles using BoxLibs ParGDB class to determine where particle should go
template<class T, unsigned Dim>
void BoxLibLayout<T, Dim>::update(AmrParticleBase< BoxLibLayout<T,Dim> >& PData,
                                  int lev_min,
                                  const ParticleAttrib<char>* canSwap)
{
    std::cout << "BoxLibLayout::update()" << std::endl;
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
            if (!Where_m(PData, ip, lev_min, theEffectiveFinestLevel)) {
                if (full_where) {
                    if (!PeriodicWhere_m(PData, ip, lev_min, theEffectiveFinestLevel)) {
                        if (lev_min != 0) {
                            if (!RestrictedWhere_m(PData, ip, nGrow))
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
            unsigned int who = this->ParticleDistributionMap(PData.Level[ip])[PData.Grid[ip]];
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
    
    std::cout << "Sent: " << sent << std::endl;
    if ( LocalNum < PData.getDestroyNum() )
        std::cout << "Can't destroy more particles" << std::endl;
    else {
        LocalNum -= PData.getDestroyNum();  // update local num
        PData.performDestroy();
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


// Function from BoxLib adjusted to work with Ippl AmrParticleBase class
//get the cell where particle is located - uses AmrParticleBase object and particle id
template <class T, unsigned Dim>
IntVect BoxLibLayout<T, Dim>::Index (AmrParticleBase< BoxLibLayout<T,Dim> >& p,
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
bool BoxLibLayout<T, Dim>::Where_m(AmrParticleBase< BoxLibLayout<T,Dim> >& p,
                                   const unsigned int ip,
                                   int lev_min,
                                   int finest_level)
{
//     BL_ASSERT(this != 0);

    if (finest_level == -1)
        finest_level = this->finestLevel();

    BL_ASSERT(finest_level <= this->finestLevel());

    std::vector< std::pair<int,Box> > isects;

    for (unsigned int lev = (unsigned)finest_level; lev >= (unsigned)lev_min; lev--)
    {
      const IntVect& iv = Index(p, ip, this->Geom(lev));

        if (lev == p.Level[ip]) { 
            // We may take a shortcut because the fact that we are here means 
            // this particle does not belong to any finer grids.
            const BoxArray& ba = this->ParticleBoxArray(p.Level[ip]);
            if (0 <= p.Grid[ip] && p.Grid[ip] < ba.size() && 
                ba[p.Grid[ip]].contains(iv)) 
            {
                return true;
            }
        }

        this->ParticleBoxArray(lev).intersections(Box(iv,iv),isects,true,0);

        if (!isects.empty())
        {
            p.Level[ip]  = lev;
            p.Grid[ip] = isects[0].first;

            return true;
        }
    }
    return false;
}


//Function from BoxLib adjusted to work with Ippl AmrParticleBase class
//Checks/sets whether the particle has crossed a periodic boundary in such a way
//that it is on levels lev_min and higher.
template <class T, unsigned Dim>
bool BoxLibLayout<T, Dim>::PeriodicWhere_m(AmrParticleBase< BoxLibLayout<T,Dim> >& p,
                                           const unsigned int ip,
                                           int lev_min,
                                           int finest_level)
{
//     BL_ASSERT(this != 0);

    if (!this->Geom(0).isAnyPeriodic()) return false;

    if (finest_level == -1)
        finest_level = this->finestLevel();

    BL_ASSERT(finest_level <= this->finestLevel());
    //
    // Create a copy "dummy" particle to check for periodic outs.
    //
    SingleParticlePos_t R = p.R[ip];

    if (PeriodicShift_m(R))
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

                p.Level[ip]  = lev;
                p.Grid[ip] = isects[0].first;

                return true;
            }
        }
    }

    return false;
}


// Function from BoxLib adjusted to work with Ippl AmrParticleBase class
// Checks/sets whether a particle is within its grid (including grow cells).
template <class T, unsigned Dim>
bool BoxLibLayout<T, Dim>::RestrictedWhere_m(AmrParticleBase< BoxLibLayout<T,Dim> >& p,
                                             const unsigned int ip,
                                             int ngrow)
{
//     BL_ASSERT(this != 0);

    const IntVect& iv = Index(p,ip,this->Geom(p.Level[ip]));

    if (BoxLib::grow(this->ParticleBoxArray(p.Level[ip])[p.Grid[ip]], ngrow).contains(iv))
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
//     BL_ASSERT(this != 0);
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

// template<class T, unsigned Dim>
// void BoxLibLayout::update(AmrParticleBase< BoxLibLayout<T,Dim> >& PData, const ParticleAttrib<char> canSwap = 0)
// {
//     
//     
// }


template <class T, unsigned Dim>
void BoxLibLayout<T, Dim>::initDefaultBox(int nGridPoints, int maxGridSize,
                                          double lower, double upper)
{
    // physical box [-0.1, 0.1]^3 (in meters)
    RealBox real_box;
    for (int d = 0; d < BL_SPACEDIM; ++d) {
        real_box.setLo(d, -lower);
        real_box.setHi(d,  upper);
    }
    
    // define underlying box for physical domain
    IntVect domain_lo(0 , 0, 0); 
    IntVect domain_hi(nGridPoints - 1, nGridPoints - 1, nGridPoints - 1); 
    const Box domain(domain_lo, domain_hi);

    // use Cartesian coordinates
    int coord = 0;

    // Dirichelt boundary conditions
    int is_per[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++) 
        is_per[i] = 0;
    
    Geometry geom;
    geom.define(domain, &real_box, coord, is_per);

    BoxArray ba;
    ba.define(domain);
    // break the BoxArrays at both levels into max_grid_size^3 boxes
    ba.maxSize(maxGridSize);
    
    DistributionMapping dmap;
    dmap.define(ba, Ippl::getNodes());
    
    // set protected ParGDB member variables
    this->m_geom.resize(1);
    this->m_geom[0] = geom;
    
    this->m_dmap.resize(1);
    this->m_dmap[0] = dmap;
    
    this->m_ba.resize(1);
    this->m_ba[0] = ba;
    
    this->m_nlevels = 1;
}

#endif
