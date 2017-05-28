
template <class T, unsigned Dim>
bool ParticleAmrLayout<T, Dim>::do_tiling = false;

template <class T, unsigned Dim>
IntVect ParticleAmrLayout<T, Dim>::tile_size   { D_DECL(1024000,8,8) };

// Function from AMReX adjusted to work with Ippl AmrParticleBase class
//get the cell where particle is located - uses AmrParticleBase object and particle id
template <class T, unsigned Dim>
IntVect ParticleAmrLayout<T, Dim>::Index (AmrParticleBase< ParticleAmrLayout<T,Dim> >& p,
					  const unsigned int ip,
					  int lev) const
{
    return Index(p.R[ip], lev);
}

//get the cell where particle is located - uses the particle position vector R
template <class T, unsigned Dim>
IntVect ParticleAmrLayout<T, Dim>::Index (SingleParticlePos_t &R,
					  int lev) const
{
    IntVect iv;
    const Geometry& geom = Geom(lev);

    D_TERM(iv[0]=floor((R[0]-geom.ProbLo(0))/geom.CellSize(0));,
           iv[1]=floor((R[1]-geom.ProbLo(1))/geom.CellSize(1));,
           iv[2]=floor((R[2]-geom.ProbLo(2))/geom.CellSize(2)););

    iv += geom.Domain().smallEnd();

    return iv;
}

template <class T, unsigned Dim>
int ParticleAmrLayout<T, Dim>::getTileIndex(const IntVect& iv, const Box& box, Box& tbx) {
    if (do_tiling == false) {
        tbx = box;
        return 0;
    } else {
        //
        // This function must be consistent with FabArrayBase::buildTileArray function!!!
        //
        auto tiling_1d = [](int i, int lo, int hi, int tilesize,
                            int& ntile, int& tileidx, int& tlo, int& thi) {
            int ncells = hi-lo+1;
            ntile = std::max(ncells/tilesize, 1);
            int ts_right = ncells/ntile;
            int ts_left  = ts_right+1;
            int nleft = ncells - ntile*ts_right;
	    int ii = i - lo;
            int nbndry = nleft*ts_left;
            if (ii < nbndry) {
                tileidx = ii / ts_left; // tiles on the left of nbndry have size of ts_left
                tlo = lo + tileidx * ts_left;
                thi = tlo + ts_left - 1;
            } else {
                tileidx = nleft + (ii-nbndry) / ts_right;  // tiles on the right: ts_right
                tlo = lo + tileidx * ts_right + nleft;
                thi = tlo + ts_right - 1;
            }
        };
        const IntVect& small = box.smallEnd();
        const IntVect& big   = box.bigEnd();
        IntVect ntiles, ivIndex, tilelo, tilehi;

        D_TERM(int iv0 = std::min(std::max(iv[0], small[0]), big[0]);,
               int iv1 = std::min(std::max(iv[1], small[1]), big[1]);,
               int iv2 = std::min(std::max(iv[2], small[2]), big[2]););

        D_TERM(tiling_1d(iv0, small[0], big[0], tile_size[0], ntiles[0], ivIndex[0], tilelo[0], tilehi[0]);,
               tiling_1d(iv1, small[1], big[1], tile_size[1], ntiles[1], ivIndex[1], tilelo[1], tilehi[1]);,
               tiling_1d(iv2, small[2], big[2], tile_size[2], ntiles[2], ivIndex[2], tilelo[2], tilehi[2]););

        tbx = Box(tilelo, tilehi);

        return D_TERM(ivIndex[0], + ntiles[0]*ivIndex[1], + ntiles[0]*ntiles[1]*ivIndex[2]);
    }
}

//sets the grid and level where particle belongs - returns false if particle is outside the domain
template <class T, unsigned Dim>
bool ParticleAmrLayout<T, Dim>::Where (AmrParticleBase< ParticleAmrLayout<T,Dim> >& p,
				       const unsigned int ip,
				       int lev_min,
                                       int lev_max,
                                       int nGrow) const
{
    BL_ASSERT(m_gdb != 0);

    if (lev_max == -1)
    lev_max = finestLevel();
  
    BL_ASSERT(lev_max <= finestLevel());
    
    BL_ASSERT(nGrow == 0 || (nGrow >= 0 && lev_min == lev_max));

    std::vector< std::pair<int,Box> > isects;

    for (unsigned int lev = (unsigned)lev_max; lev >= (unsigned)lev_min; lev--)
    {
        const IntVect& iv = Index(p, ip, lev);
        const BoxArray& ba = ParticleBoxArray(lev);
        BL_ASSERT(ba.ixType().cellCentered());

	if (lev == p.m_lev[ip]) { 
            // The fact that we are here means this particle does not belong to any finer grids.
	    if (0 <= p.m_grid[ip] && p.m_grid[ip] < ba.size())
	    {
		const Box& bx = ba.getCellCenteredBox(p.m_grid[ip]);
                const Box& gbx = amrex::grow(bx,nGrow);
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
        
        ba.intersections(Box(iv, iv), isects, true, nGrow);
        
        if (!isects.empty())
        {
            p.m_lev[ip]  = lev;
            p.m_grid[ip] = isects[0].first;

            return true;
        }
    }
    return false;
}

//Function from AMReX adjusted to work with Ippl AmrParticleBase class
//Checks/sets whether the particle has crossed a periodic boundary in such a way
//that it is on levels lev_min and higher.
template <class T, unsigned Dim>
bool ParticleAmrLayout<T, Dim>::EnforcePeriodicWhere (AmrParticleBase<ParticleAmrLayout<T,Dim> >& p,
                                                      const unsigned int ip,
                                                      int lev_min,
                                                      int lev_max) const
{
    BL_ASSERT(m_gdb != 0);

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
        std::vector< std::pair<int,Box> > isects;

        for (int lev = lev_max; lev >= lev_min; lev--)
        {
            const IntVect& iv = Index(R, lev);
            const BoxArray& ba = ParticleBoxArray(lev);
            
            ba.intersections(Box(iv,iv),isects,true,0);

            if (!isects.empty())
            {
                D_TERM(p.R[ip][0] = R[0];,
                       p.R[ip][1] = R[1];,
                       p.R[ip][2] = R[2];);

                p.m_lev[ip]  = lev;
                p.m_grid[ip] = isects[0].first;

                return true;
            }
        }
    }

    return false;
}

// Function from AMReX adjusted to work with Ippl AmrParticleBase class
// Returns true if the particle was shifted.
template <class T, unsigned Dim>
bool ParticleAmrLayout<T, Dim>::PeriodicShift (SingleParticlePos_t R) const
{
  //
    // This routine should only be called when Where() returns false.
    //
    BL_ASSERT(m_gdb != 0);
    //
    // We'll use level 0 stuff since ProbLo/ProbHi are the same for every level.
    //
    const Geometry& geom    = Geom(0);
    const Box&      dmn     = geom.Domain();
    const IntVect&  iv      = Index(R, 0);
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
