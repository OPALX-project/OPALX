#ifndef AMROPAL_H
#define AMROPAL_H

#include <AmrCore.H>
#include "AmrPartBunch.h"

#include <MultiFabUtil.H>

#include <memory>

/*!
 * @file AmrOpal.h
 * @authors Matthias Frey
 *          Ann Almgren
 *          Weiqun Zhang
 * @date October 2016, LBNL
 * @details This class implements the abstract base class
 * AmrCore of the BoxLib library. It defines functionalities
 * to do Amr calculations.
 * @brief Concrete AMR implementation
 */

/// Concrete AMR implementation
class AmrOpal : public AmrCore {
    
private:
//     typedef Array<std::unique_ptr<MultiFab> > mfs_mt; ///< instead of using PArray<MultiFab>
    typedef PArray<MultiFab > mfs_mt;
    
    
    
//     typedef PArray<MultiFab> mp_mt;
    
public:
    /*!
     * Create an AMR object.
     * @param rb is the physical domain
     * @param max_level_in is the max. number of allowed AMR levels
     * @param n_cell_in is the number of grid cells at the coarsest level
     * @param coord is the coordinate system (0: cartesian)
     * @param bunch is the particle bunch
     */
    AmrOpal(const RealBox* rb, int max_level_in, const Array<int>& n_cell_in, int coord, PartBunchBase* bunch);
    
    /*!
     * Create an AMR object.
     * @param rb is the physical domain
     * @param max_level_in is the max. number of allowed AMR levels
     * @param n_cell_in is the number of grid cells at the coarsest level
     * @param coord is the coordinate system (0: cartesian)
     */
    AmrOpal(const RealBox* rb, int max_level_in, const Array<int>& n_cell_in, int coord);
    
    virtual ~AmrOpal();     ///< does nothing
    
    void initBaseLevel();   ///< defines BoxArray, DistributionMapping of level 0
    
    /*!
     * @param lbase is the level on top of which a new level is created
     * @param time not used
     */
    void regrid (int lbase, Real time);
    
    /*!
     * Update the grids and the distributionmapping for a specific level
     * @param lev is the current level
     * @param time not used
     * @param new_grids are the new created grids for this level
     * @param new_dmap is the new distribution to processors
     */
    void RemakeLevel (int lev, Real time,
                      const BoxArray& new_grids, const DistributionMapping& new_dmap);
    
    /*!
     * Create completeley new grids for a level
     * @param lev is the current level
     * @param time not used
     * @param new_grids are the new created grids for this level
     * @param new_dmap is the new distribution to processors
     */
    void MakeNewLevel (int lev, Real time,
                       const BoxArray& new_grids, const DistributionMapping& new_dmap);
    
    
    void ClearLevel(int lev);
    
    void setBunch(AmrPartBunch* bunch) {
        bunch_m = bunch;
    }
    
    /*!
     * Print the number of particles per cell (minimum and maximum)
     */
    void info() {
        for (int i = 0; i < finest_level; ++i)
            std::cout << "density level " << i << ": "
#ifdef UNIQUE_PTR
                      << nPartPerCell_m[i]->min(0) << " "
                      << nPartPerCell_m[i]->max(0) << std::endl;
#else
                      
                      << nPartPerCell_m[i].min(0) << " "
                      << nPartPerCell_m[i].max(0) << std::endl;
#endif
                      
    }
    
    /*!
     * Write a timestamp file for displaying with yt.
     */
    void writePlotFileYt(std::string filename, int step);
    
    /*!
     * Write a timestamp file for displaying with AmrVis.
     */
    void writePlotFile(std::string filename, int step);
    
    mfs_mt* getPartPerCell() {
        return &nPartPerCell_m;
    }
    
    void assignDensity() {
        
        for (int i = 0; i < finest_level; ++i)
#ifdef UNIQUE_PTR
            chargeOnGrid_m[i]->setVal(0.0);
#else
            chargeOnGrid_m[i].setVal(0.0);
#endif
        
        bunch_m->AssignDensity(0, false, chargeOnGrid_m, 0, 1, finest_level);
        
//         double assign_sum = 0.0;
//         double charge_sum = 0.0;
//         std::cout << "---------------------------------------------" << std::endl;
//         std::cout << "          CHARGE CONSERVATION TEST           " << std::endl;
//         for (int i = 0; i <= finest_level; ++i) {
//             Real charge = bunch_m->sumParticleMass(0 /*attribute*/, i /*level*/);
//             Real invVol = (*(Geom(i).CellSize()) * *(Geom(i).CellSize()) * *(Geom(i).CellSize()) );
//             std::cout << "dx * dy * dz = " << invVol << std::endl;
//             assign_sum += chargeOnGrid_m[i]->sum() * invVol;
//             std::cout << "Level " << i << " MultiFab sum * dx * dy * dz: " << chargeOnGrid_m[i]->sum() * invVol
//                       << " Charge sum: " << charge
//                       << " Spacing: " << *(Geom(i).CellSize()) << std::endl;
//             charge_sum += charge;
//         }
//         std::cout << "Total charge: " << assign_sum << " " << charge_sum << std::endl;
//         std::cout << "---------------------------------------------" << std::endl;
    }

protected:
    /*!
     * Is called in the AmrCore function for performing tagging.
     */
    virtual void ErrorEst(int lev, TagBoxArray& tags, Real time, int /*ngrow*/) override;
    
private:
    AmrPartBunch* bunch_m;      ///< Particle bunch
    mfs_mt/*mp_mt*/ nPartPerCell_m;      ///< used in tagging.
    mfs_mt chargeOnGrid_m;
    
};

#endif