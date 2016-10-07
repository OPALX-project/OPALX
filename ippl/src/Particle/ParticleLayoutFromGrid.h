//
// C++ Interface: ParticleLayoutFromGrid
//
// Description: 
//
//
// Author: Roman Geus <geus@maxwell>, (C) 2005
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef ParticleLayoutFromGrid_H
#define ParticleLayoutFromGrid_H

#include "Ippl.h"
#include "extpde.h"

class ParticleLayoutFromGrid : public ParticleLayout<double,3> {
public:
    typedef int pair_t;
    typedef pair_t* pair_iterator;
    typedef ParticleAttrib<SingleParticlePos_t> ParticlePos_t;
    typedef ParticleAttrib<Index_t>             ParticleIndex_t;
    typedef Vektor<double,3>                    Vector_t;
    /**
     * Constructor.
     * @param grid      Grid.
     * @param domain    Domain on which the grid is based.
     * @param local_min Corner of axiparallel box defining local domain (with minimal coordinates)
     * @param local_max Corner of axiparallel box defining local domain (with maximal coordinates)
     */
    ParticleLayoutFromGrid(Grid* grid, Domain* domain, D3vector& local_min, D3vector& local_max) :
        grid_(grid),
        geom_domain_(domain),
        my_corner_min_(local_min),
        my_corner_max_(local_max)
    {}    
    /**
     * Does "book--keeping" after particles have been integrated one time step.
     *
     * "Book--keeping" includes 
     *    - Detecting particles which have moved across a boundary and 
     *      taking appropriate measures. (not yet implemented)
     *    - Redistributing particles which have moved out of the domain
     *      assigned to the local processor.
     *    - Updating the total particle count.
     * @param particles Particle container.
     */
    void update(IpplParticleBase< ParticleLayoutFromGrid >& particles);    
    /**
     * Test if x is inside the domain assigned to my processor.
     * @param x Position to test.
     * @return  True if inside.
     */
    inline bool is_local_pos(const Vector_t& x) const {
      const D3vector y(x[0],x[1],x[2]);      
      return is_inside_box(my_corner_min_, my_corner_max_, y);
    }    
    /**
     * Test if x is inside the domain assigned to my processor.
     * @param x Position to test.
     * @return  True if inside.
     */
    inline bool is_local_pos(const D3vector& x) const {
      return is_inside_box(my_corner_min_, my_corner_max_, x);
    }    

    /**
     * Test if x is inside the geometric domain
     * @param x Position to test.
     * @return  True if inside.
     */
    bool is_in_domain(const D3vector& x) const {
        return geom_domain_->point_in_domain(x);
    }
    /**
     * Test if x is inside the geometric domain
     * @param x Position to test.
     * @return  True if inside.
     */
    bool is_in_domain(const Vector_t& x) const {
      const D3vector y(x[0],x[1],x[2]);
      return geom_domain_->point_in_domain(y);
    }

    /**
     * Test whether x is inside axiparalllel box defined by corner_min and corner_max.
     * @param corner_min Coordinates of box corner with minimal values.
     * @param corner_max Coordinates of box corner with maximal values.
     * @param x          Position to test.
     * @return           True if x is inside.
     */
    inline static bool is_inside_box(const D3vector& corner_min, const D3vector& corner_max, const D3vector& x) {
        return corner_min[0] <= x[0] && x[0] < corner_max[0] &&
               corner_min[1] <= x[1] && x[1] < corner_max[1] &&
               corner_min[2] <= x[2] && x[2] < corner_max[2];
    }
   
private:
    /**
     * Dummy function for boundary treatment (currently not used)
     * @param R Particle positions.
     */
    void apply_bconds(ParticlePos_t& R);
    /**
     * Go through all local particles, and send particles which are no longer in 
     * the local bounding box to the corresponding processors.
     *
     * @param particles Particle container
     * @return Number of locally stored particles after redistribution
     */
    size_t redistribute_particles(IpplParticleBase< ParticleLayoutFromGrid >& particles);
    /**
     * Pointer to the grid data structure.
     */
    Grid* grid_;
    /**
     * Pointer to the geometric domain data structure.
     */
    Domain* geom_domain_;
    /**
     * Corner of axiparallel box defining local domain (with minimal coordinates).
     */
    D3vector my_corner_min_;
    /**
     * Corner of axiparallel box defining local domain (with maximal coordinates).
     */
    D3vector my_corner_max_;
};

#endif
