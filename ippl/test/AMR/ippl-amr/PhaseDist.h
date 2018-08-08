#ifndef PHASE_DISTRIBUTION_H
#define PHASE_DISTRIBUTION_H

#include "../AmrOpal.h"

#include <boost/serialization/vector.hpp>
#include <boost/mpi.hpp>

class Particle {
    
public:
    typedef Vektor<double, AMREX_SPACEDIM> Vector_t;
    typedef std::vector<double> vector_t;
    
public:
    
    Particle() {}
    
    Particle(const Vector_t& x,
             const Vector_t& v,
             const double& q)
        : q_m(q)
    {
        for (std::size_t i = 0; i < AMREX_SPACEDIM; ++i) {
            x_m.push_back(x(i));
            v_m.push_back(v(i));
        }
    }
    
    Particle(const vector_t& x,
             const vector_t& v,
             const double& q)
        : x_m(x)
        , v_m(v)
        , q_m(q)
    { }
    
    Particle(const Particle& p)
        : x_m(p.x_m)
        , v_m(p.v_m)
        , q_m(p.q_m)
    { }
    
private:
    friend class boost::serialization::access;
    
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & x_m;
        ar & v_m;
        ar & q_m;
    }
    
public:
    vector_t x_m;
    vector_t v_m;
    double   q_m;
};

// BOOST_IS_MPI_DATATYPE(Particle)

// namespace boost { namespace mpi {
//   template <>
//   struct is_mpi_datatype<Particle> : mpl::true_ { };
// } }


class PhaseDist {
    
public:
    typedef AmrOpal::amrplayout_t amrplayout_t;
    typedef Vektor<double, AMREX_SPACEDIM> Vector_t;
    
public:
    
    PhaseDist();
    
    PhaseDist(const Vector_t& left,
              const Vector_t& right,
              const Vector_t& nx,
              const Vector_t& vmin,
              const Vector_t& vmax,
              const Vector_t& nv);
    
    void define(const Vector_t& left,
                const Vector_t& right,
                const Vector_t& nx,
                const Vector_t& vmin,
                const Vector_t& vmax,
                const Vector_t& nv);
    
    void deposit(const ParticleAttrib<double>& q,
                 const amrplayout_t::ParticlePos_t& x,
                 const amrplayout_t::ParticlePos_t& v,
                 std::size_t localnum);
private:
    
    void fill_m(const ParticleAttrib<double>& q,
                const amrplayout_t::ParticlePos_t& x,
                const amrplayout_t::ParticlePos_t& v,
                std::size_t localnum);
    
    amrex::IntVect index_m (const Particle::vector_t& x,
                            const Particle::vector_t& v) const;
    
    int where_m(const Particle::vector_t& x,
                const Particle::vector_t& v);
    
    void redistribute_m();
    
    
    
private:
    Vector_t left_m;
    Vector_t nx_m;
    Vector_t dx_m;
    
    Vector_t vmin_m;
    Vector_t nv_m;
    Vector_t dv_m;
    
    amrex::BoxArray fba_m;
    amrex::DistributionMapping fdmap_m;
    amrex::MultiFab fmf_m;
    
    std::vector<Particle> particles_m;
};

#endif
