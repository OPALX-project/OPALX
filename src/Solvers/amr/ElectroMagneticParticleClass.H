#include "Particles.H"
#include "RealBox.H"

template <int N>
class EMPC
: public ParticleContainer<N>
{

public:

    EMPC (Amr* amr)
        : ParticleContainer<N>((ParGDBBase*)amr->GetParGDB()), m_amr(amr) {}

    void FillEvec (MultiFab& field_vector, int level, Real betaC,             int e_comp, int b_comp);
    void MKD      (int level, Real timestep, int field_ngrow    , int v_comp, int e_comp, int b_comp);
    void MK       (int level, Real timestep                     , int v_comp, int e_comp, int b_comp);

    Real sumParticleMomenta (int bz_comp);

    //
    // The type of Particles we hold.
    //
    typedef Particle<N> ParticleType;

    //
    // We want to store the particles on a level by level and grid by grid basis.  This will
    // make accessing them and doing operations on them more memory efficient since most of our
    // operations on particles are done on a level by level basis or grid by grid basis.
    //
    typedef typename ParticleContainer<N>::PBox PBox;

    //
    // A level of particles is stored in a map indexed by the grid number.
    //
    typedef typename  ParticleContainer<N>::PMap PMap;

protected:
    
private:
    Amr* m_amr;

};

template <int N>
void
EMPC<N>::FillEvec (MultiFab& field_vector, int lev, Real betaC, int e_comp, int b_comp)
{
    // If there are no particles at this level
    if (lev >= int(this->m_particles.size()))
        return;

    // Make sure we have enough room to store the charge, velocity and field.
    BL_ASSERT(N >= 2*BL_SPACEDIM+1);

    PMap&      pmap      = this->m_particles[lev];

    MultiFab* gv_pointer;
    if (((ParGDB*)this->m_amr->GetParGDB())->OnSameGrids(lev, field_vector))
    {
        gv_pointer = &field_vector;
    }
    else 
    {
//         gv_pointer = new MultiFab(this->m_amr->ParticleBoxArray(lev),
//                                   field_vector.nComp(),field_vector.nGrow(),Fab_allocate);
        for (MFIter mfi(*gv_pointer); mfi.isValid(); ++mfi)
            gv_pointer->setVal(0.);
        gv_pointer->copy(field_vector,0,0,field_vector.nComp());
        gv_pointer->FillBoundary(); // DO WE NEED GHOST CELLS FILLED ???
    }

    for (typename PMap::iterator pmap_it = pmap.begin(), 
         pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it)
    {
        const int        grid = pmap_it->first;
        PBox&            pbox = pmap_it->second;
        const int        n    = pbox.size();
        const FArrayBox& efab = (*gv_pointer)[grid];

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < n; i++)
        {
            ParticleType& p = pbox[i];

            if (p.m_id > 0)
            {
                BL_ASSERT(p.m_grid == grid);
                Real field[BL_SPACEDIM];
                ParticleBase::GetGravity(efab, this->m_amr->Geom(lev), p, field);

                // These are the components of the electric field
                p.m_data[e_comp  ] = field[0];
                p.m_data[e_comp+1] = field[1];
                p.m_data[e_comp+2] = field[2];

                // These are the components of the magnetic field
                p.m_data[b_comp  ] = -betaC * p.m_data[e_comp+1];
                p.m_data[b_comp+1] =  betaC * p.m_data[e_comp  ];
                p.m_data[b_comp+2] = 0.;
            }
        }
    }
}

//
// This takes as input the field vector stored in the particle data
//
template <int N>
void
EMPC<N>::MKD (int lev, Real dt, int field_ngrow, int v_comp, int e_comp, int b_comp)
{
    BL_ASSERT(N >= BL_SPACEDIM+1);
    BL_ASSERT(lev >= 0);

    //If there are no particles at this level
    if (lev >= int(this->m_particles.size()))
        return;

    Real v_cross_b[3];
    const Real half_dt       = 0.5 * dt;
    PMap&      pmap          = this->m_particles[lev];

    for (typename PMap::iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it)
    {
        PBox&            pbox = pmap_it->second;
        const int        n    = pbox.size();

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < n; i++)
        {
            ParticleType& p = pbox[i];

            if (p.m_id <= 0) continue;

            //
            // Compute the cross-product V x B before we update V
            //
            v_cross_b[0] = p.m_data[v_comp+1] * p.m_data[b_comp+2] - p.m_data[v_comp+2] * p.m_data[b_comp+1];
            v_cross_b[1] = p.m_data[v_comp+2] * p.m_data[b_comp+0] - p.m_data[v_comp  ] * p.m_data[b_comp+2];
            v_cross_b[2] = p.m_data[v_comp  ] * p.m_data[b_comp+1] - p.m_data[v_comp+1] * p.m_data[b_comp  ];

            //
            // Update the velocity: V^{n+1/2} = V^n + (dt/2) E^n
            //
            p.m_data[v_comp  ] += half_dt * (p.m_data[e_comp  ] + v_cross_b[0]);
            p.m_data[v_comp+1] += half_dt * (p.m_data[e_comp+1] + v_cross_b[1]);
            p.m_data[v_comp+2] += half_dt * (p.m_data[e_comp+2] + v_cross_b[2]);

            //
            // Now update x^new = x^old + dt * V^{n+1/2}
            //
            D_TERM(p.m_pos[0] += dt * p.m_data[1];,
                   p.m_pos[1] += dt * p.m_data[2];,
                   p.m_pos[2] += dt * p.m_data[3];);
        } 
    }

    if (lev > 0 && this->m_amr->subCycle())
    {
        for (typename PMap::iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it)
        {
            PBox&            pbox = pmap_it->second;
            const int        n    = pbox.size();
 
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (int i = 0; i < n; i++)
            {
                ParticleType& p = pbox[i];
                if (p.m_id <= 0) continue;

                // Move the particle to the proper ghost cell. 
                //      and remove any *ghost* particles that have gone too far
                // Note that this should only negate ghost particles, not real particles.
                if (!ParticleBase::RestrictedWhere(p,(ParGDBBase*)this->m_amr->GetParGDB(), field_ngrow-2))
                {
                    // Assert that the particle being removed is a ghost particle;
                    // the ghost particle is no longer in relevant ghost cells for this grid.
                    if (p.m_id == GhostParticleID)
                    {
                        p.m_id = -1;
                    }
                    else
                    {
                        std::cout << "Oops -- removing particle " << p.m_id << std::endl;
                        BoxLib::Error("Trying to get rid of a non-ghost particle in moveKickDrift");
                    }
                }
            }
        }
    }
}

//
// This takes as input the field vector stored in the particle data
//
template <int N>
void
EMPC<N>::MK (int lev, Real dt, int v_comp, int e_comp, int b_comp)
{
    BL_ASSERT(N >= BL_SPACEDIM+1);
    BL_ASSERT(lev >= 0 && lev < this->m_particles.size());

    Real v_cross_b[3];
    const Real half_dt   = 0.5 * dt;
    PMap&      pmap      = this->m_particles[lev];

    for (typename PMap::iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it)
    {
        PBox&            pbox = pmap_it->second;
        const int        n    = pbox.size();

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < n; i++)
        {
            ParticleType& p = pbox[i];

            if (p.m_id > 0)
            {
                //
                // Compute the cross-product V x B before we update V
                //
                v_cross_b[0] = p.m_data[v_comp+1] * p.m_data[b_comp+2] - p.m_data[v_comp+2] * p.m_data[b_comp+1];
                v_cross_b[1] = p.m_data[v_comp+2] * p.m_data[b_comp+0] - p.m_data[v_comp  ] * p.m_data[b_comp+2];
                v_cross_b[2] = p.m_data[v_comp  ] * p.m_data[b_comp+1] - p.m_data[v_comp+1] * p.m_data[b_comp  ];

                //
                // Update the velocity: v^{n+1} = v^{n+1/2} + (dt/2) E^{n+1}
                //
                p.m_data[v_comp  ] += half_dt * (p.m_data[e_comp  ] + v_cross_b[0]);
                p.m_data[v_comp+1] += half_dt * (p.m_data[e_comp+1] + v_cross_b[1]);
                p.m_data[v_comp+2] += half_dt * (p.m_data[e_comp+2] + v_cross_b[2]);

            }
        }
    }
}

//
// This sums the z-momenta of all the particles and divides by the number of particles
//
template <int N>
Real
EMPC<N>::sumParticleMomenta (int bz_comp)
{

    // This is the total number of particles on *all* processors
    long npart = this->TotalNumberOfParticles(true,false);
    std::cout << "Summing the z-momentum over " << npart << " particles " << std::endl;

    Real sum = 0.;

    for (int lev = 0; lev <= this->m_amr->finestLevel(); lev++)
    {
        PMap& pmap = this->m_particles[lev];

        for (typename PMap::iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it)
        {
            PBox&     pbx = pmap_it->second;
            const int n    = pbx.size();

            for (int i = 0; i < n; i++)
            {
                ParticleType& p = pbx[i];
                if (p.m_id > 0)
                    sum += p.m_data[bz_comp];
            }
        }
    }

    ParallelDescriptor::ReduceRealSum(sum);
    sum /= npart;
    return sum;
}
