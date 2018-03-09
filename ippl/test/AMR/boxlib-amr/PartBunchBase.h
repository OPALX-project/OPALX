#ifndef PARTBUNCHBASE__H
#define PARTBUNCHBASE__H

#include <Ippl.h>
#include <iomanip>

// dimension of our positions
const unsigned Dim = 3;

// some typedefs
typedef ParticleSpatialLayout<double,Dim>::SingleParticlePos_t Vector_t;
typedef ParticleSpatialLayout<double,Dim> playout_t;
typedef UniformCartesian<Dim,double> Mesh_t;
typedef Cell                                       Center_t;
typedef CenteredFieldLayout<Dim, Mesh_t, Center_t> FieldLayout_t;
typedef Field<double, Dim, Mesh_t, Center_t>       Field_t;
typedef Field<Vector_t, Dim, Mesh_t, Center_t>     VField_t;
typedef IntCIC IntrplCIC_t;

const double pi = acos(-1.0);


/*!
 * @file PartBunchBase.h
 * @details In order to work properly with OPAL it has
 * to provide the functionalities of the
 * ippl/src/Particle/IpplParticleBase
 * and src/Classic/Algorithms/PartBunch class.
 * @authors Matthias Frey \n
 *          Andreas Adelmann \n
 *          Ann Almgren \n
 *          Weiqun Zhang
 * @date LBNL, October 2016
 * @brief Abstract base class for particle bunches
 */


/// Abstract base class for particle bunches.
class PartBunchBase {

public:

    /// Does nothing.
    virtual ~PartBunchBase() {}


    /// After tracking the particle need to be redistributed to guarantee load balancing
    virtual void myUpdate() = 0;

    /// Each processor creates m particles (see IpplParticleBase)
    /*! @param m is the number of particles to be created.
     */
    virtual void create(size_t m) = 0;

    /// Local number of particles (see IpplParticleBase)
    /*! @returns the local number of particles
     */
    virtual size_t getLocalNum() const = 0;

    /// Total number of particles (see IpplParticleBase)
    /*! @returns the total number of particles
     */
    virtual size_t getTotalNum() const = 0;

    // Load balancing statistics
    /*!
     * Obtain load balancing information and print
     * it to the shell.
     */
    virtual void gatherStatistics() = 0;

    /// Lower corner of domain.
    /*! @returns the lower corner of the box domain.
     */
    virtual Vector_t getRMin() = 0;

    /// Upper corner of domain.
    /*! @returns the upper corner of the box domain.
     */
    virtual Vector_t getRMax() = 0;

    /// Mesh spacing
    /*! @returns the mesh spacing.
     */
    virtual Vector_t getHr() = 0;

    /// Access the position of a particle
    /*! @returns the i-th particle position (x, y, z)
     */
    virtual Vector_t getR(int i) = 0;

    /// Access the charge-to-mass ratio of a particle
    /*! @returns the charge-to-mass ratio of the i-th particle
     */
    virtual double getQM(int i) = 0;
    
    /// Access the mass of a particle
    /*! @returns the mass of the i-th particle
     */
    virtual double getMass(int i) = 0;
    
    /// Access the velocity of a particle
    /*! @returns the i-th particle velocity (px, py, pz)
     */
    virtual Vector_t getP(int i) = 0;

    /// Access the electric field at the particle location
    /*! @returns (Ex, Ey, Ez) at i-th particle location
     */
    virtual Vector_t getE(int i) = 0;

    /// Access the magnetic field a the particle location
    /*! @returns (Bx, By, Bz) at i-th particle location
     */
    virtual Vector_t getB(int i) = 0;

    /// Set the particle position
    /*!
     * @param pos is the position (x, y, z)
     * @param i specifies the i-th particle
     */
    virtual void setR(Vector_t pos, int i) = 0;

    /// Set the particle charge-to-mass ratio
    /*!
     * @param q is the new charge-to-mass ratio
     * @param i specifies the i-th particle
     */
    virtual void setQM(double q, int i) = 0;
    
    /// Set the particle charge-to-mass ratio
    /*!
     * @param m is the new mass
     * @param i specifies the i-th particle
     */
    virtual void setMass(double q, int i) = 0;
    
    /// Set the particle velocity
    /*!
     * @param v is the velocity (vx, vy, vz)
     * @param i specifies the i-th particle
     */
    virtual void setP(Vector_t v, int i) = 0;

    /// Set the E-field
    /*!
     * @param Ef is the electric field (Ex, Ey, Ez)
     * @param i specifies the i-th particle
     */
    virtual void setE(Vector_t Ef, int i) = 0;

    /// Set the particle position
    /*!
     * @param Bf is the magnetic field (Bx, By, Bz)
     * @param i specifies the i-th particle
     */
    virtual void setB(Vector_t Bf, int i) = 0;

    virtual double scatter() = 0;
    virtual void initFields() = 0;
    virtual void gatherCIC() = 0;

    /// Writes the particles (x, y, z) to the shell.
    void print();

    /// Delete all particles and free allocated memory
    virtual void destroyAll() = 0;
    
    /*!
     * As gatherStatistics() but it dumps the data to file.
     * @param filename is the pathname
     */
    virtual void dumpStatistics(const std::string& filename) = 0;
};


inline
void PartBunchBase::print() {
    for (int p = 0; p < Ippl::getNodes(); ++p) {

        if ( p == Ippl::myNode() )
            for (size_t i = 0; i < getLocalNum(); ++i) {
                std::cout << std::setprecision(8) << getR(i)(0) << " "
                          << getR(i)(1) << " " << getR(i)(2) << std::endl;
            }

        Ippl::Comm->barrier();
    }
}


#endif