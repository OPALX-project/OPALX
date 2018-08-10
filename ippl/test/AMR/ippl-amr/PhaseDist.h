#ifndef PHASE_DISTRIBUTION_H
#define PHASE_DISTRIBUTION_H

#include "../AmrOpal.h"

#include <boost/serialization/vector.hpp>
#include <boost/mpi.hpp>

#include "Particle/ParticleLayout.h"


template <class T, unsigned Dim>
class DummyLayout : public ParticleLayout<T, Dim>
{

public:
    // pair iterator definition ... this layout does not allow for pairlists
    typedef int pair_t;
    typedef pair_t* pair_iterator;
    typedef typename ParticleLayout<T, Dim>::SingleParticlePos_t
    SingleParticlePos_t;
    typedef typename ParticleLayout<T, Dim>::Index_t Index_t;
  
    // type of attributes this layout should use for position and ID
    typedef ParticleAttrib<SingleParticlePos_t> ParticlePos_t;
    typedef ParticleAttrib<Index_t>             ParticleIndex_t;
    
    
    void update(IpplParticleBase< DummyLayout<T,Dim> >& PData, 
                const ParticleAttrib<char>* canSwap=0)
    {
        throw std::runtime_error("Not provided");
    }
};



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


class PhaseDist : public IpplParticleBase<DummyLayout<double, AMREX_SPACEDIM> > {
    
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
              const Vector_t& nv,
              const int& maxgrid);
    
    void define(const Vector_t& left,
                const Vector_t& right,
                const Vector_t& nx,
                const Vector_t& vmin,
                const Vector_t& vmax,
                const Vector_t& nv,
                const int& maxgrid);
    
    void deposit(const ParticleAttrib<double>& q,
                 const amrplayout_t::ParticlePos_t& x,
                 const amrplayout_t::ParticlePos_t& v,
                 std::size_t localnum);
    
    void write(const std::string& fname);
    
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
    
    
    void amrex_deposit_m();
    
#ifdef USE_IPPL
    void ippl_deposit_m();
    
    void ippl_init_m();
#endif
    
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
    
    
    typedef Cell                                        Center_t;
    typedef IntCIC                                      IntrplCIC_t;
    typedef UniformCartesian<2, double>                 Mesh2d_t;
    typedef CenteredFieldLayout<2, Mesh2d_t, Center_t>  FieldLayout2d_t;
    typedef Field<double, 2, Mesh2d_t, Center_t>        Field2d_t;
    
    Mesh2d_t  mesh2d_m;
    std::shared_ptr<FieldLayout2d_t> layout2d_m;
    Field2d_t field2d_m;
    
    
    ParticleAttrib<Vektor<double,2> > xphase_m;
    ParticleAttrib<double> q_m;
};

#endif
