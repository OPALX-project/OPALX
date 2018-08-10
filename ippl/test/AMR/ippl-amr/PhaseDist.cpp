#include "PhaseDist.h"

#include <boost/mpi/collectives.hpp>
#include <boost/mpi/status.hpp>

#include <algorithm>

PhaseDist::PhaseDist()
    : left_m(Vector_t())
    , nx_m(Vector_t())
    , dx_m(Vector_t())
    , vmin_m(Vector_t())
    , nv_m(Vector_t())
    , dv_m(Vector_t())
{
#ifdef USE_IPPL
    this->addAttribute(xphase_m);
    this->addAttribute(q_m);
#endif
}

PhaseDist::PhaseDist(const Vector_t& left,
              const Vector_t& right,
              const Vector_t& nx,
              const Vector_t& vmin,
              const Vector_t& vmax,
              const Vector_t& nv,
              const int& maxgrid)
{
    this->define(left, right, nx,
                 vmin, vmax, nv,
                 maxgrid);
}


void PhaseDist::define(const Vector_t& left,
                       const Vector_t& right,
                       const Vector_t& nx,
                       const Vector_t& vmin,
                       const Vector_t& vmax,
                       const Vector_t& nv,
                       const int& maxgrid)
{
    left_m = left;
    nx_m = nx;
    dx_m = (right - left) / nx;
    
    vmin_m = vmin;
    nv_m = nv;
    dv_m = (vmax - vmin) / nv;
    
    amrex::IntVect low(D_DECL(0, 0, 0));
    amrex::IntVect high(D_DECL(nx_m[0] - 1, nv_m[0] - 1, 0));    
    amrex::Box bx(low, high);
    
    fba_m.define(bx);
    fba_m.maxSize( maxgrid );
    
    fdmap_m.define(fba_m, amrex::ParallelDescriptor::NProcs());
    
#ifdef USE_IPPL
    ippl_init_m();
#endif
    
    fmf_m.define(fba_m, fdmap_m, 1, 1);
}


void PhaseDist::deposit(const ParticleAttrib<double>& q,
                        const amrplayout_t::ParticlePos_t& x,
                        const amrplayout_t::ParticlePos_t& v,
                        std::size_t localnum)
{
    this->fill_m(q, x, v, localnum);
    
    this->redistribute_m();
    
#ifdef USE_IPPL
    this->ippl_deposit_m();
#else
    this->amrex_deposit_m();
#endif
    
    particles_m.clear();
}


void PhaseDist::write(const std::string& fname)
{    
#ifdef USE_IPPL
    this->ippl_write_m(fname);
#else
    this->amrex_write_m(fname);
#endif
}


void PhaseDist::fill_m(const ParticleAttrib<double>& q,
                       const amrplayout_t::ParticlePos_t& x,
                       const amrplayout_t::ParticlePos_t& v,
                       std::size_t localnum)
{
    particles_m.clear();
    
    for (std::size_t i = 0; i < localnum; ++i)
        particles_m.push_back( Particle(x[i], v[i], q[i]) );
}

amrex::IntVect PhaseDist::index_m (const Particle::vector_t& x,
                                   const Particle::vector_t& v) const
{
    amrex::IntVect iv;
        
    D_TERM(iv[0]=std::floor((x[0]-left_m[0])/dx_m[0]);,
           iv[1]=std::floor((v[0]-vmin_m[0])/dv_m[0]);,
           iv[2]=0;);
//         iv += geom.Domain().smallEnd();
    return iv;
}

    
int PhaseDist::where_m(const Particle::vector_t& x,
                       const Particle::vector_t& v)
{
    std::vector< std::pair<int,amrex::Box> > isects;
    
    int nGrow = 0;
    int grid = 0;
    
    const amrex::IntVect iv = index_m(x, v);
    const amrex::Box& bx = fba_m.getCellCenteredBox(grid);
    const amrex::Box& gbx = amrex::grow(bx,nGrow);
    if (gbx.contains(iv))
    {
        return grid;
    }
    
    fba_m.intersections(amrex::Box(iv, iv), isects, true, nGrow);
    
    if (!isects.empty())
    {
        grid = isects[0].first;
    }
    
    return grid;
}


void PhaseDist::redistribute_m()
{
    if ( Ippl::getNodes() < 2 )
        return;
    
    boost::mpi::environment env;
    boost::mpi::communicator world(Ippl::getComm(), boost::mpi::comm_duplicate);
    
    typedef std::multimap<std::size_t, std::size_t> multimap_t;
    
    multimap_t p2n; //node ID, particle
    std::size_t N = Ippl::getNodes();
    
    int *msgsend = new int[N];
    std::fill(msgsend, msgsend+N, 0);
    int *msgrecv = new int[N];
    std::fill(msgrecv, msgrecv+N, 0);
    
    std::vector<Particle> tmp;
    
    for (std::size_t i = 0; i < particles_m.size(); ++i)
    {
        int grid = this->where_m(particles_m[i].x_m, particles_m[i].v_m);
        
        const int cpu = fdmap_m[grid];
        
        if ( cpu != Ippl::myNode() ) {
            msgsend[cpu] = 1;
            p2n.insert(std::pair<std::size_t, std::size_t>(cpu, i));
        } else {
            tmp.push_back( particles_m[i] );
        }
    }
    
    
    boost::mpi::all_reduce(world, msgsend, N, msgrecv, std::plus<int>());
    
    std::vector<boost::mpi::request> requests;
    
    typename multimap_t::iterator it = p2n.begin();
    
    int tag = 11;
    
    while ( it != p2n.end() ) {
        std::size_t cur_destination = it->first;
        
        std::vector<Particle> psend;
        
        for (; it != p2n.end() && it->first == cur_destination; ++it) {
            psend.push_back( particles_m[it->second] );
        }
        
        requests.push_back( world.isend(cur_destination, tag, psend) );
    }
    
    
    for (int k = 0; k < msgrecv[Ippl::myNode()]; ++k)
    {
        boost::mpi::status stat = world.probe(boost::mpi::any_source, tag);
        
        std::vector<Particle> precv;
        
        world.recv(stat.source(), stat.tag(), precv);
        
        std::copy(precv.begin(), precv.end(), std::back_inserter(tmp));
    }
    
//     boost::mpi::wait_all(requests, requests.size());
    particles_m.swap(tmp);
}


void PhaseDist::amrex_deposit_m()
{
    amrex::Periodicity periodicity(amrex::IntVect(D_DECL(0, 0, 0)));
    
    int grid = 0;
    double lxv[2] = { 0.0, 0.0 };
    double wxv_hi[2] = { 0.0, 0.0 };
    double wxv_lo[2] = { 0.0, 0.0 };
    int ij[2] = { 0, 0 };
    
    fmf_m.setVal(0.0, 1);
    
    amrex::MultiFab fmf(fba_m, fdmap_m, 1, 2);
    fmf.setVal(0.0, 2);
    
    double inv_dx[2] = { 1.0 / dx_m[0], 1.0 / dv_m[0] };
    
    for (std::vector<Particle>::iterator it = particles_m.begin();
         it != particles_m.end(); ++it)
    {
        lxv[0] = ( it->x_m[0] - left_m[0] ) * inv_dx[0] + 0.5;
        ij[0] = lxv[0];
        wxv_hi[0] = lxv[0] - ij[0];
        wxv_lo[0] = 1.0 - wxv_hi[0];
        
        lxv[1] = ( it->v_m[0] - vmin_m[0] ) * inv_dx[1] + 0.5;
        ij[1] = lxv[1];
        wxv_hi[1] = lxv[1] - ij[1];
        wxv_lo[1] = 1.0 - wxv_hi[1];
        
        int& i = ij[0];
        int& j = ij[1];
        
        amrex::IntVect i1(D_DECL(i-1, j-1, 0));
        amrex::IntVect i2(D_DECL(i-1, j,   0));
        amrex::IntVect i3(D_DECL(i,   j-1, 0));
        amrex::IntVect i4(D_DECL(i,   j,   0));
        
        grid = this->where_m(it->x_m, it->v_m);
        
        amrex::FArrayBox& fab = fmf[grid];
        
        fab(i1, 0) += wxv_lo[0] * wxv_lo[1] * it->q_m;
        fab(i2, 0) += wxv_lo[0] * wxv_hi[1] * it->q_m;
        fab(i3, 0) += wxv_hi[0] * wxv_lo[1] * it->q_m;
        fab(i4, 0) += wxv_hi[0] * wxv_hi[1] * it->q_m;
    }
    
    fmf.SumBoundary(periodicity);
    
    const amrex::Real vol = dx_m[0] * dv_m[0];
    fmf.mult(-1.0/vol, 0, 1, fmf.nGrow()); // minus --> to make positive (charges < 0)
    
    amrex::MultiFab::Copy(fmf_m, fmf, 0, 0, 1, 0);
    
    if ( fmf_m.contains_nan() || fmf_m.contains_inf() )
        throw std::runtime_error("\033[1;31mError: NANs or INFs on charge grid.\033[0m");
}


void PhaseDist::amrex_write_m(const std::string& fname) {
    
    std::ofstream out;
    out.precision(10);
    out.setf(std::ios::scientific, std::ios::floatfield);
    
    bool first = true;
    
    for (amrex::MFIter mfi(fmf_m); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        const amrex::FArrayBox& fab = fmf_m[mfi];
        
        for (int p = 0; p < amrex::ParallelDescriptor::NProcs(); ++p) {
            
            if ( p == amrex::ParallelDescriptor::MyProc() ) {
            
                if ( p == 0 && first) {
                    first = false;
                    out.open(fname.c_str(), std::ios::out);
                    out << "x, vx, f" << std::endl;
                } else {
                    out.open(fname.c_str(), std::ios::app);
                }
                
                for (int i = bx.loVect()[0]; i <= bx.hiVect()[0]; ++i) {
                    for (int j = bx.loVect()[1]; j <= bx.hiVect()[1]; ++j) {
                        
                        amrex::IntVect iv(D_DECL(i, j, 0));
                        
                        out << (i + 0.5) * dx_m[0]
                            << ", "
                            << (j + 0.5) * dv_m[0] + vmin_m[0]
                            << ", "
                            << fab(iv, 0)
                            << std::endl;
                    }
                }
                out.close();
            }
            amrex::ParallelDescriptor::Barrier();
        }
    }
}


#ifdef USE_IPPL
void PhaseDist::ippl_deposit_m()
{
    std::size_t new_localnum = particles_m.size();
    
    xphase_m.create(new_localnum);
    q_m.create(new_localnum);
    
    for (std::size_t i = 0; i < new_localnum; ++i)
    {
        xphase_m[i](0) = particles_m[i].x_m[0];
        xphase_m[i](1) = particles_m[i].v_m[0];
        q_m[i]         = particles_m[i].q_m;
    }
    
    
    field2d_m = 0;
    
    q_m.scatter(field2d_m, xphase_m, IntrplCIC_t());
    
    xphase_m.destroy(new_localnum, 0);
    q_m.destroy(new_localnum, 0);
}


void PhaseDist::ippl_init_m() {
    amrex::IntVect low(D_DECL(0, 0, 0));
    amrex::IntVect high(D_DECL(nx_m[0] - 1, nv_m[0] - 1, 0));    
    amrex::Box bx(low, high);
    amrex::IntVect iv = bx.size();
    
    NDIndex<2> domain;
    domain[0] = Index(0, iv[0]);
    domain[1] = Index(0, iv[1]);
    
    
    double spacings[2] = { dx_m[0], dv_m[0] };    
    Vektor<double,2> origin = { left_m[0], vmin_m[0] };
    
    mesh2d_m = Mesh2d_t(domain, spacings, origin);
    
    BConds<double, 2, Mesh2d_t, Center_t> bc;
    
    if (Ippl::getNodes() > 1) {
        bc[0] = new ParallelInterpolationFace<double,2,Mesh2d_t, Center_t>(0);
        bc[1] = new ParallelInterpolationFace<double,2,Mesh2d_t, Center_t>(1);
        bc[2] = new ParallelInterpolationFace<double,2,Mesh2d_t, Center_t>(2);
        bc[3] = new ParallelInterpolationFace<double,2,Mesh2d_t, Center_t>(3);
    }
    else {
        bc[0] = new InterpolationFace<double,2,Mesh2d_t, Center_t>(0);
        bc[1] = new InterpolationFace<double,2,Mesh2d_t, Center_t>(1);
        bc[2] = new InterpolationFace<double,2,Mesh2d_t, Center_t>(2);
        bc[3] = new InterpolationFace<double,2,Mesh2d_t, Center_t>(3);
    }
    
    auto pmap = fdmap_m.ProcessorMap();
    
    std::vector< NDIndex<2> > regions;
    std::vector< int > nodes;
    for (uint i = 0; i < pmap.size(); ++i) {
        amrex::Box bx = fba_m[i];
        
        NDIndex<2> range;
        for (int j = 0; j < 2; ++j)
            range[j] = Index(bx.smallEnd(j), bx.bigEnd(j));
        
        regions.push_back( range );
        nodes.push_back( pmap[i] );
    }
    
    layout2d_m = std::shared_ptr<FieldLayout2d_t>(new FieldLayout2d_t(mesh2d_m,
                                                                      &regions[0],
                                                                      &regions[0] + regions.size(),
                                                                      &nodes[0],
                                                                      &nodes[0] + nodes.size()));
    
    field2d_m.initialize(mesh2d_m, *(layout2d_m.get()), GuardCellSizes<2>(1), bc);
}


void PhaseDist::ippl_write_m(const std::string& fname) {
    std::ofstream out;
    out.precision(10);
    out.setf(std::ios::scientific, std::ios::floatfield);
    
    if ( Ippl::myNode() == 0 ) {
        out.open(fname.c_str(), std::ios::out);
        out << "x, vx, f" << std::endl;
    }
    
    NDIndex<2> lDom = layout2d_m->getDomain();    
    
    for (int i=lDom[0].first(); i<=lDom[0].last(); i++) {
        for (int j=lDom[1].first(); j<=lDom[1].last(); j++) {
             // multiply by -1, since charge is < 0
            double fval = -field2d_m[i][j].get();
            
            if ( Ippl::myNode() == 0 ) {
                out << (i + 0.5) * dx_m[0] << ", "
                    << (j + 0.5) * dv_m[0] + vmin_m[0] << ", "
                    << fval<< std::endl;
            }
        }
    }
    
    if ( Ippl::myNode() == 0 )
        out.close();
}
#endif
