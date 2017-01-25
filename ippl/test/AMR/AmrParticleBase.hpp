#include "Ippl.h"

#include <Particles.H>
#include <ParmParse.H>
#include <limits>

template<class PLayout>
void AmrParticleBase<PLayout>::CIC_Cells_Fracs_Basic(const SingleParticlePos_t &R, 
						     const Real* plo, 
						     const Real* dx, 
						     Real* fracs,  
						     IntVect* cells)
{

  const Real len[BL_SPACEDIM] = { D_DECL((R[0]-plo[0])/dx[0] + Real(0.5),
					 (R[1]-plo[1])/dx[1] + Real(0.5),
					 (R[2]-plo[2])/dx[2] + Real(0.5)) };

  const IntVect cell(D_DECL(floor(len[0]), floor(len[1]), floor(len[2])));
  
  const Real frac[BL_SPACEDIM] = { D_DECL(len[0]-cell[0], len[1]-cell[1], len[2]-cell[2]) };

  ParticleBase::CIC_Fracs(frac, fracs);
  ParticleBase::CIC_Cells(cell, cells);

}

template<class PLayout>
int AmrParticleBase<PLayout>::CIC_Cells_Fracs (const SingleParticlePos_t &R,
					       const Real*         plo,
					       const Real*         dx_geom,
					       const Real*         dx_part,
					       Array<Real>&        fracs,
					       Array<IntVect>&     cells)
{
    BL_PROFILE("AmrParticleBase::CIC_Cells_Fracs()");
    if (dx_geom == dx_part)
    {
        const int M = D_TERM(2,+2,+4);
        fracs.resize(M);
        cells.resize(M);
        AmrParticleBase::CIC_Cells_Fracs_Basic(R,plo,dx_geom,fracs.dataPtr(),cells.dataPtr());
        return M;
    }
    //
    // The first element in fracs and cells is the lowest corner, the last is the highest.
    //
    const Real hilen[BL_SPACEDIM] = { D_DECL((R[0]-plo[0]+dx_part[0]/2)/dx_geom[0],
                                             (R[1]-plo[1]+dx_part[1]/2)/dx_geom[1],
                                             (R[2]-plo[2]+dx_part[2]/2)/dx_geom[2]) };

    const Real lolen[BL_SPACEDIM] = { D_DECL((R[0]-plo[0]-dx_part[0]/2)/dx_geom[0],
                                             (R[1]-plo[1]-dx_part[1]/2)/dx_geom[1],
                                             (R[2]-plo[2]-dx_part[2]/2)/dx_geom[2]) };

    const IntVect hicell(D_DECL(floor(hilen[0]), floor(hilen[1]), floor(hilen[2])));
    
    const IntVect locell(D_DECL(floor(lolen[0]), floor(lolen[1]), floor(lolen[2])));
    
    const Real cell_density = D_TERM(dx_geom[0]/dx_part[0],*dx_geom[1]/dx_part[1],*dx_geom[2]/dx_part[2]);
    
    const int M = D_TERM((hicell[0]-locell[0]+1),*(hicell[1]-locell[1]+1),*(hicell[2]-locell[2]+1));

    fracs.resize(M);
    cells.resize(M);
    //
    // This portion might be slightly inefficient. Feel free to redo it if need be.
    //
    int i = 0;
    for (int zi = locell[2]; zi <= hicell[2]; zi++)
    {
        const Real zf = std::min(hilen[2]-zi,Real(1))-std::max(lolen[2]-zi,Real(0));
        for (int yi = locell[1]; yi <= hicell[1]; yi++)
        {
            const Real yf = std::min(hilen[1]-yi,Real(1))-std::max(lolen[1]-yi,Real(0));
            for (int xi = locell[0]; xi <= hicell[0]; xi++)
            {
                cells[i][0] = xi;
                cells[i][1] = yi;
                cells[i][2] = zi;
                fracs[i] = zf * yf * (std::min(hilen[0]-xi,Real(1))
				      -std::max(lolen[0]-xi,Real(0))) * cell_density;
                i++;
            }
        }
    }

    return M;
}



