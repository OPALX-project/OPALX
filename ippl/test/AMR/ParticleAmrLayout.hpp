template <class T, unsigned Dim>
IntVect ParticleAmrLayout<T, Dim>::Index (AmrParticleBase< ParticleAmrLayout<T,Dim> >& p,
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

template <class T, unsigned Dim>
IntVect ParticleAmrLayout<T, Dim>::Index (SingleParticlePos_t &R,
					  const Geometry&     geom)
{
    IntVect iv;

    D_TERM(iv[0]=floor((R[0]-geom.ProbLo(0))/geom.CellSize(0));,
           iv[1]=floor((R[1]-geom.ProbLo(1))/geom.CellSize(1));,
           iv[2]=floor((R[2]-geom.ProbLo(2))/geom.CellSize(2)););

    iv += geom.Domain().smallEnd();

    return iv;
}

template <class T, unsigned Dim>
bool ParticleAmrLayout<T, Dim>::Where (AmrParticleBase< ParticleAmrLayout<T,Dim> >& p,
				       const unsigned int ip,
				       const ParGDBBase* gdb,
				       int lev_min,
				       int finest_level)
{
    BL_ASSERT(gdb != 0);

    if (finest_level == -1)
        finest_level = gdb->finestLevel();

    BL_ASSERT(finest_level <= gdb->finestLevel());

    std::vector< std::pair<int,Box> > isects;

    for (unsigned int lev = (unsigned)finest_level; lev >= (unsigned)lev_min; lev--)
    {
      const IntVect& iv = Index(p, ip, gdb->Geom(lev));

	if (lev == p.m_lev[ip]) { 
            // We may take a shortcut because the fact that we are here means 
            // this particle does not belong to any finer grids.
	    const BoxArray& ba = gdb->ParticleBoxArray(p.m_lev[ip]);
	    if (0 <= p.m_grid[ip] && p.m_grid[ip] < ba.size() && 
		ba[p.m_grid[ip]].contains(iv)) 
	    {
		return true;
	    }
	}

        gdb->ParticleBoxArray(lev).intersections(Box(iv,iv),isects,true,0);

        if (!isects.empty())
        {
            p.m_lev[ip]  = lev;
            p.m_grid[ip] = isects[0].first;

            return true;
        }
    }
    return false;
}

template <class T, unsigned Dim>
bool ParticleAmrLayout<T, Dim>::PeriodicWhere (AmrParticleBase<ParticleAmrLayout<T,Dim> >& p,
					       const unsigned int ip,
					       const ParGDBBase* gdb,
					       int lev_min,
					       int finest_level)
{
    BL_ASSERT(gdb != 0);

    if (!gdb->Geom(0).isAnyPeriodic()) return false;

    if (finest_level == -1)
        finest_level = gdb->finestLevel();

    BL_ASSERT(finest_level <= gdb->finestLevel());
    //
    // Create a copy "dummy" particle to check for periodic outs.
    //
    SingleParticlePos_t R = p.R[ip];

    if (PeriodicShift(R, gdb))
    {
        std::vector< std::pair<int,Box> > isects;

        for (int lev = finest_level; lev >= lev_min; lev--)
        {
            const IntVect& iv = Index(R, gdb->Geom(lev));

            gdb->ParticleBoxArray(lev).intersections(Box(iv,iv),isects,true,0);

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

template <class T, unsigned Dim>
bool ParticleAmrLayout<T, Dim>::RestrictedWhere (AmrParticleBase<ParticleAmrLayout<T,Dim> >& p,
						 const unsigned int ip,
						 const ParGDBBase* gdb,
						 int ngrow)
{
    BL_ASSERT(gdb != 0);

    const IntVect& iv = Index(p,ip,gdb->Geom(p.m_lev[ip]));

    if (BoxLib::grow(gdb->ParticleBoxArray(p.m_lev[ip])[p.m_grid[ip]], ngrow).contains(iv))
    {
        return true;
    }

    return false;
}

template <class T, unsigned Dim>
bool ParticleAmrLayout<T, Dim>::PeriodicShift (SingleParticlePos_t R,
					       const ParGDBBase* gdb)
{
  //
    // This routine should only be called when Where() returns false.
    //
    BL_ASSERT(gdb != 0);
    //
    // We'll use level 0 stuff since ProbLo/ProbHi are the same for every level.
    //
    const Geometry& geom    = gdb->Geom(0);
    const Box&      dmn     = geom.Domain();
    const IntVect&  iv      = Index(R,gdb->Geom(0));
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
