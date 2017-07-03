#ifndef BoxLibLayout_HPP
#define BoxLibLayout_HPP

#include "BoxLibLayout.h"

#include "Message/Formatter.h"
#include "Utilities/OpalException.h"

#include <cmath>

template <class T, unsigned Dim>
const Vector_t BoxLibLayout<T, Dim>::lowerBound = - Vector_t(1.025, 1.025, 1.025);


template <class T, unsigned Dim>
const Vector_t BoxLibLayout<T, Dim>::upperBound = Vector_t(1.025, 1.025, 1.025);


template<class T, unsigned Dim>
BoxLibLayout<T, Dim>::BoxLibLayout()
    : ParGDB(),
      finestLevel_m(0),
      maxLevel_m(0),
      forbidTransform_m(false),
      refRatio_m(0),
      scale_m(1.0)
{
    /* FIXME There might be a better solution
     * 
     * 
     * Figure out the number of grid points in each direction
     * such that all processes have some data at the beginning
     * 
     * ( nGridPoints / maxGridSize ) ^3 = max. #procs
     * 
     */
    int nProcs = Ippl::getNodes();
    int maxGridSize = 16;
    
    int nGridPoints = std::ceil( std::cbrt( nProcs ) ) * maxGridSize;
    
    this->initDefaultBox(nGridPoints, maxGridSize);
}


template<class T, unsigned Dim>
BoxLibLayout<T, Dim>::BoxLibLayout(int nGridPoints, int maxGridSize)
    : ParGDB(),
      finestLevel_m(0),
      maxLevel_m(0),
      forbidTransform_m(false),
      refRatio_m(0),
      scale_m(1.0)
{
    this->initDefaultBox(nGridPoints, maxGridSize);
}


template<class T, unsigned Dim>
BoxLibLayout<T, Dim>::BoxLibLayout(const AmrGeometry_t &geom,
                                   const AmrProcMap_t &dmap,
                                   const AmrGrid_t &ba)
    : ParGDB(geom, dmap, ba),
      finestLevel_m(0),
      maxLevel_m(0),
      forbidTransform_m(false),
      refRatio_m(0),
      scale_m(1.0)
{ }


template<class T, unsigned Dim>
BoxLibLayout<T, Dim>::BoxLibLayout(const AmrGeomContainer_t &geom,
                                   const AmrProcMapContainer_t &dmap,
                                   const AmrGridContainer_t &ba,
                                   const AmrIntArray_t &rr)
    : ParGDB(geom, dmap, ba, rr),
      finestLevel_m(0),
      maxLevel_m(0),
      forbidTransform_m(false),
      refRatio_m(0),
      scale_m(1.0)
{ }


template<class T, unsigned Dim>
void BoxLibLayout<T, Dim>::update(IpplParticleBase< BoxLibLayout<T,Dim> >& PData,
                                  const ParticleAttrib<char>* canSwap)
{
    /* Exit since we need AmrParticleBase with grids and levels for particles for this layout
     * if IpplParticleBase is used something went wrong
     */
    throw OpalException("BoxLibLayout::update(IpplParticleBase, ParticleAttrib) ",
                        "Wrong update method called.");
}

// // Function from AMReX adjusted to work with Ippl AmrParticleBase class
// // redistribute the particles using BoxLibs ParGDB class to determine where particle should go
template<class T, unsigned Dim>
void BoxLibLayout<T, Dim>::update(AmrParticleBase< BoxLibLayout<T,Dim> >& PData,
                                  int lev_min, int lev_max,
                                  const ParticleAttrib<char>* canSwap)
{
    // in order to avoid transforms when already done
    if ( !forbidTransform_m ) {
        // we need to update on Amr domain, has to be undone at end of function
        this->domainMapping(PData);
    }
    
    int nGrow = 0;

    unsigned N = Ippl::getNodes();
    unsigned myN = Ippl::myNode();
    
    int theEffectiveFinestLevel = this->finestLevel();
    while (!this->LevelDefined(theEffectiveFinestLevel)) {
        theEffectiveFinestLevel--;
    }
    
    if (lev_max == -1)
        lev_max = theEffectiveFinestLevel;
    
    //loop trough the particles and assigne the grid and level where each particle belongs
    size_t LocalNum = PData.getLocalNum();
    
    auto& LocalNumPerLevel = PData.getLocalNumPerLevel();
    
    if ( LocalNum != LocalNumPerLevel.getLocalNumAllLevel() )
        throw OpalException("BoxLibLayout::update()",
                            "Local #particles disagrees with sum over levels");
    
    std::multimap<unsigned, unsigned> p2n; //node ID, particle 

    int *msgsend = new int[N];
    std::fill(msgsend, msgsend+N, 0);
    int *msgrecv = new int[N];
    std::fill(msgrecv, msgrecv+N, 0);

    unsigned sent = 0;
    size_t lBegin = LocalNumPerLevel.begin(lev_min);
    size_t lEnd   = LocalNumPerLevel.end(lev_max);
  
    //loop trough particles and assign grid and level to each particle
    //if particle doesn't belong to this process save the index of the particle to be sent
    for (unsigned int ip = lBegin; ip < lEnd; ++ip) {
        bool particleLeftDomain = false;
        
        // old level
        const size_t& lold = PData.Level[ip];
        
//         /*
//          * AMReX sets m_grid = -1 and m_lev = -1
//          */
//         PData.Level[ip] = -1;
//         PData.Grid[ip] = -1;
        
        //check to which level and grid the particle belongs to
        locateParticle(PData, ip, lev_min, lev_max, nGrow, particleLeftDomain);
        
        if ( !particleLeftDomain ) {
            // The owner of the particle is the CPU owning the finest grid
            // in state data that contains the particle.
            const size_t& lnew = PData.Level[ip];
            
            const unsigned int who = ParticleDistributionMap(lnew)[PData.Grid[ip]];
            
            --LocalNumPerLevel[lold];
            
            if (who != myN) {
                // we lost the particle to another process
                msgsend[who] = 1;
                p2n.insert(std::pair<unsigned, unsigned>(who, ip));
                sent++;
            } else {
                /* if we still own the particle it may have moved to
                 * another level
                 */
                ++LocalNumPerLevel[lnew];
            }
        } else {
            /* a particle left the domain
             * This should NEVER happen! An exception is thrown in
             * BoxLibLayout::locateParticle().
             */
            --LocalNumPerLevel[lold];
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
    if ( LocalNum < PData.getDestroyNum() )
        throw OpalException("BoxLibLayout::update()",
                            "Rank " + std::to_string(myN) +
                            " can't destroy more particles than possessed.");
    else {
        LocalNum -= PData.getDestroyNum();  // update local num
        PData.performDestroy();
    }
    
    for (int lev = lev_min; lev <= lev_max; ++lev) {
        if ( LocalNumPerLevel[lev] < 0 )
            throw OpalException("BoxLibLayout::update()",
                                "Negative particle level count.");
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
                 * pEndIdx is the last index of the new particle data
                 */
                size_t pBeginIdx = LocalNum;
                
                LocalNum += PData.getSingleMessage(*msg);
                
                size_t pEndIdx = LocalNum;
                
                for (size_t idx = pBeginIdx; idx < pEndIdx; ++idx)
                    ++LocalNumPerLevel[ PData.Level[idx] ];
                
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
    PData.setTotalNum(TotalNum);    // set the total atom count
    PData.setLocalNum(LocalNum);    // set the number of local atoms
    
    if ( !forbidTransform_m ) {
        // undo domain transformation
        this->domainMapping(PData, true);
    }
}


// Function from AMReX adjusted to work with Ippl AmrParticleBase class
//get the cell where particle is located - uses AmrParticleBase object and particle id
template <class T, unsigned Dim>
typename BoxLibLayout<T, Dim>::AmrIntVect_t
    BoxLibLayout<T, Dim>::Index(AmrParticleBase< BoxLibLayout<T,Dim> >& p,
                                const unsigned int ip,
                                int lev) const
{
    return Index(p.R[ip], lev);
}

//get the cell where particle is located - uses the particle position vector R
template <class T, unsigned Dim>
typename BoxLibLayout<T, Dim>::AmrIntVect_t
    BoxLibLayout<T, Dim>::Index(SingleParticlePos_t &R,
                                int lev) const
{
    AmrIntVect_t iv;
    const AmrGeometry_t& geom = Geom(lev);

    D_TERM(iv[0]=floor((R[0]-geom.ProbLo(0))/geom.CellSize(0));,
           iv[1]=floor((R[1]-geom.ProbLo(1))/geom.CellSize(1));,
           iv[2]=floor((R[2]-geom.ProbLo(2))/geom.CellSize(2)););

    iv += geom.Domain().smallEnd();

    return iv;
}


// template <class T, unsigned Dim>
// int BoxLibLayout<T, Dim>::getTileIndex(const AmrIntVect_t& iv, const Box& box, Box& tbx) {
//     if (do_tiling == false) {
//         tbx = box;
//         return 0;
//     } else {
//         //
//         // This function must be consistent with FabArrayBase::buildTileArray function!!!
//         //
//         auto tiling_1d = [](int i, int lo, int hi, int tilesize,
//                             int& ntile, int& tileidx, int& tlo, int& thi) {
//             int ncells = hi-lo+1;
//             ntile = std::max(ncells/tilesize, 1);
//             int ts_right = ncells/ntile;
//             int ts_left  = ts_right+1;
//             int nleft = ncells - ntile*ts_right;
// 	    int ii = i - lo;
//             int nbndry = nleft*ts_left;
//             if (ii < nbndry) {
//                 tileidx = ii / ts_left; // tiles on the left of nbndry have size of ts_left
//                 tlo = lo + tileidx * ts_left;
//                 thi = tlo + ts_left - 1;
//             } else {
//                 tileidx = nleft + (ii-nbndry) / ts_right;  // tiles on the right: ts_right
//                 tlo = lo + tileidx * ts_right + nleft;
//                 thi = tlo + ts_right - 1;
//             }
//         };
//         const AmrIntVect_t& small = box.smallEnd();
//         const AmrIntVect_t& big   = box.bigEnd();
//         AmrIntVect_t ntiles, ivIndex, tilelo, tilehi;
// 
//         D_TERM(int iv0 = std::min(std::max(iv[0], small[0]), big[0]);,
//                int iv1 = std::min(std::max(iv[1], small[1]), big[1]);,
//                int iv2 = std::min(std::max(iv[2], small[2]), big[2]););
// 
//         D_TERM(tiling_1d(iv0, small[0], big[0], tile_size[0], ntiles[0], ivIndex[0], tilelo[0], tilehi[0]);,
//                tiling_1d(iv1, small[1], big[1], tile_size[1], ntiles[1], ivIndex[1], tilelo[1], tilehi[1]);,
//                tiling_1d(iv2, small[2], big[2], tile_size[2], ntiles[2], ivIndex[2], tilelo[2], tilehi[2]););
// 
//         tbx = Box(tilelo, tilehi);
// 
//         return D_TERM(ivIndex[0], + ntiles[0]*ivIndex[1], + ntiles[0]*ntiles[1]*ivIndex[2]);
//     }
// }


//sets the grid and level where particle belongs - returns false if particle is outside the domain
template <class T, unsigned Dim>
bool BoxLibLayout<T, Dim>::Where (AmrParticleBase< BoxLibLayout<T,Dim> >& p,
				       const unsigned int ip,
				       int lev_min,
                                       int lev_max,
                                       int nGrow) const
{

    if (lev_max == -1)
        lev_max = finestLevel();
  
    BL_ASSERT(lev_max <= finestLevel());
    
    BL_ASSERT(nGrow == 0 || (nGrow >= 0 && lev_min == lev_max));

    std::vector< std::pair<int, AmrBox_t> > isects;

    for (unsigned int lev = (unsigned)lev_max; lev >= (unsigned)lev_min; lev--)
    {
        const AmrIntVect_t& iv = Index(p, ip, lev);
        const AmrGrid_t& ba = ParticleBoxArray(lev);
        BL_ASSERT(ba.ixType().cellCentered());

	if (lev == p.Level[ip]) { 
            // The fact that we are here means this particle does not belong to any finer grids.
	    if (0 <= p.Grid[ip] && p.Grid[ip] < ba.size())
	    {
		const AmrBox_t& bx = ba.getCellCenteredBox(p.Grid[ip]);
                const AmrBox_t& gbx = amrex::grow(bx,nGrow);
                if (gbx.contains(iv))
                {
//                     if (bx != pld.m_gridbox || !pld.m_tilebox.contains(iv)) {
//                         pld.m_tile = getTileIndex(iv, bx, pld.m_tilebox);
//                         pld.m_gridbox = bx;
//                     }
                    return true;
                }
            }
        }
        
        ba.intersections(AmrBox_t(iv, iv), isects, true, nGrow);
        
        if (!isects.empty())
        {
            p.Level[ip]  = lev;
            p.Grid[ip] = isects[0].first;

            return true;
        }
    }
    return false;
}


//Function from AMReX adjusted to work with Ippl AmrParticleBase class
//Checks/sets whether the particle has crossed a periodic boundary in such a way
//that it is on levels lev_min and higher.
template <class T, unsigned Dim>
bool BoxLibLayout<T, Dim>::EnforcePeriodicWhere (AmrParticleBase< BoxLibLayout<T,Dim> >& p,
                                                      const unsigned int ip,
                                                      int lev_min,
                                                      int lev_max) const
{
    if (!Geom(0).isAnyPeriodic()) return false;

    if (lev_max == -1)
        lev_max = finestLevel();

    BL_ASSERT(lev_max <= finestLevel());
    //
    // Create a copy "dummy" particle to check for periodic outs.
    //
    SingleParticlePos_t R = p.R[ip];

    if (PeriodicShift(R))
    {
        std::vector< std::pair<int, AmrBox_t> > isects;

        for (int lev = lev_max; lev >= lev_min; lev--)
        {
            const AmrIntVect_t& iv = Index(R, lev);
            const AmrGrid_t& ba = ParticleBoxArray(lev);
            
            ba.intersections(AmrBox_t(iv,iv),isects,true,0);

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


// Function from AMReX adjusted to work with Ippl AmrParticleBase class
// Returns true if the particle was shifted.
template <class T, unsigned Dim>
bool BoxLibLayout<T, Dim>::PeriodicShift (SingleParticlePos_t R) const
{
  //
    // This routine should only be called when Where() returns false.
    //
    //
    // We'll use level 0 stuff since ProbLo/ProbHi are the same for every level.
    //
    const AmrGeometry_t& geom    = Geom(0);
    const AmrBox_t&      dmn     = geom.Domain();
    const AmrIntVect_t&  iv      = Index(R, 0);
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


template <class T, unsigned Dim>
void BoxLibLayout<T, Dim>::locateParticle(
    AmrParticleBase< BoxLibLayout<T,Dim> >& p,
    const unsigned int ip,
    int lev_min, int lev_max, int nGrow,
    bool &particleLeftDomain) const
{
    bool outside = D_TERM( p.R[ip](0) <  AmrGeometry_t::ProbLo(0)
                        || p.R[ip](0) >= AmrGeometry_t::ProbHi(0),
                        || p.R[ip](1) <  AmrGeometry_t::ProbLo(1)
                        || p.R[ip](1) >= AmrGeometry_t::ProbHi(1),
                        || p.R[ip](2) <  AmrGeometry_t::ProbLo(2)
                        || p.R[ip](2) >= AmrGeometry_t::ProbHi(2));
        
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
            
            /* We shouldn't lose particles since they are mapped to be within
             * [-1, 1]^3.
             */
            throw OpalException("BoxLibLayout::locateParticle()",
                                "We're losing particles although we shouldn't");
            
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
        throw OpalException("BoxLibLayout::locateParticle()",
                            "Invalide particle.");
    }
}


template <class T, unsigned Dim>
void BoxLibLayout<T, Dim>::initDefaultBox(int nGridPoints, int maxGridSize)
{
    // physical box (in meters)
    AmrDomain_t real_box;
    for (int d = 0; d < BL_SPACEDIM; ++d) {
        real_box.setLo(d, lowerBound[d]);
        real_box.setHi(d, upperBound[d]);
    }
    
    // define underlying box for physical domain
    AmrIntVect_t domain_lo(0 , 0, 0); 
    AmrIntVect_t domain_hi(nGridPoints - 1, nGridPoints - 1, nGridPoints - 1); 
    const AmrBox_t domain(domain_lo, domain_hi);

    // use Cartesian coordinates
    int coord = 0;

    // Dirichlet boundary conditions
    int is_per[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++) 
        is_per[i] = 0;
    
    AmrGeometry_t geom;
    geom.define(domain, &real_box, coord, is_per);

    AmrGrid_t ba;
    ba.define(domain);
    // break the BoxArrays at both levels into max_grid_size^3 boxes
    ba.maxSize(maxGridSize);
    
    AmrProcMap_t dmap;
    dmap.define(ba, Ippl::getNodes());
    
    // set protected ParGDB member variables
    this->m_geom.resize(1);
    this->m_geom[0] = geom;
    
    this->m_dmap.resize(1);
    this->m_dmap[0] = dmap;
    
    this->m_ba.resize(1);
    this->m_ba[0] = ba;
    
    this->m_nlevels = ba.size();
    
//     std::cout << "----------------------------" << std::endl
//               << "BoxLibLayout::initDefaultBox" << std::endl
//               << this->m_geom[0] << std::endl
//               << this->m_ba[0] << std::endl
//               << "----------------------------" << std::endl;
}


template <class T, unsigned Dim>
void BoxLibLayout<T, Dim>::setFinestLevel(int finestLevel) {
    finestLevel_m = finestLevel;
}


template <class T, unsigned Dim>
void BoxLibLayout<T, Dim>::setMaxLevel(int maxLevel) {
    maxLevel_m = maxLevel;
}


template <class T, unsigned Dim>
void BoxLibLayout<T, Dim>::setForbidTransform(bool forbidTransform) {
    forbidTransform_m = forbidTransform;
}


template <class T, unsigned Dim>
bool BoxLibLayout<T, Dim>::isForbidTransform() const {
    return forbidTransform_m;
}


// overwritten functions
template <class T, unsigned Dim>
bool BoxLibLayout<T, Dim>::LevelDefined (int level) const {
    return level <= maxLevel_m && !m_ba[level].empty() && !m_dmap[level].empty();
}


template <class T, unsigned Dim>
int BoxLibLayout<T, Dim>::finestLevel () const {
    return finestLevel_m;
}


template <class T, unsigned Dim>
int BoxLibLayout<T, Dim>::maxLevel () const {
    return maxLevel_m;
}


template <class T, unsigned Dim>
typename BoxLibLayout<T, Dim>::AmrIntVect_t
    BoxLibLayout<T, Dim>::refRatio (int level) const
{
    return refRatio_m[level];
}


template <class T, unsigned Dim>
int BoxLibLayout<T, Dim>::MaxRefRatio (int level) const {
    int maxval = 0;
    for (int n = 0; n<BL_SPACEDIM; n++) 
        maxval = std::max(maxval, refRatio_m[level][n]);
    return maxval;
}


template <class T, unsigned Dim>
const double& BoxLibLayout<T, Dim>::domainMapping(
    AmrParticleBase< BoxLibLayout<T,Dim> >& PData,
    bool inverse)
{
    Vector_t rmin, rmax;
    bounds(PData.R, rmin, rmax);
    
    double scale = scale_m;
    
    if ( !inverse ) {
        Vector_t tmp = Vector_t(std::max( std::abs(rmin[0]), std::abs(rmax[0]) ),
                                std::max( std::abs(rmin[1]), std::abs(rmax[1]) ),
                                std::max( std::abs(rmin[2]), std::abs(rmax[2]) )
                               );
        
        scale = std::max( tmp[0], tmp[1] );
        scale = std::max( scale, tmp[2] );
    }
    
    Vector_t vscale = Vector_t(scale, scale, scale);
    
    for (unsigned int i = 0; i < PData.getLocalNum(); ++i)
        PData.R[i] /= vscale;
    
    
    scale_m = 1.0 / scale;
    
    return scale_m;
}

#endif
