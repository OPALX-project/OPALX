#include "AmrPartBunch.h"

#include <memory>

// ----------------------------------------------------------------------------
// STATIC MEMBER VARIABLES
// ----------------------------------------------------------------------------
size_t AmrPartBunch::nAttributes = 5;

// ----------------------------------------------------------------------------
// PUBLIC MEMBER FUNCTIONS
// ----------------------------------------------------------------------------

AmrPartBunch::AmrPartBunch(const Geometry& geom,
                           const DistributionMapping & dmap,
                           const BoxArray & ba)
    : ParticleContainer(geom, dmap, ba),
      nLocalParticles_m(0)
{ }

AmrPartBunch::AmrPartBunch(const Array<Geometry>& geom,
                           const Array<DistributionMapping>& dmap,
                           const Array<BoxArray>& ba,
                           const Array<int>& rr)
    : ParticleContainer(geom, dmap, ba, rr),
      nLocalParticles_m(0)
{ }


AmrPartBunch::~AmrPartBunch()
{ }

// ----------------------------------------------------------------------------
// PRIVATE MEMBER FUNCTIONS
// ----------------------------------------------------------------------------

void AmrPartBunch::buildIndexMapping_m() {
    
    idxMap_m.clear();
    
    int i = 0;
    nLocalParticles_m = 0;
    for (unsigned int l = 0; l < m_particles.size(); ++l) {
        for (unsigned int g = 0; g <  m_particles[l].size(); ++g) {
            
            nLocalParticles_m += m_particles[l][g].size();
            
            for (unsigned int d = 0; d < m_particles[l][g].size(); ++d) {
                idxMap_m[i] = std::tuple<int, int, int>(l,g,d);
                ++i;
            }
        }
    }
}


// ----------------------------------------------------------------------------
// INHERITED MEMBER FUNCTIONS
// ----------------------------------------------------------------------------


void AmrPartBunch::myUpdate() {
    
    // 2nd argument true due to periodic shift
    Redistribute(false, true);
    
    buildIndexMapping_m();
}


void AmrPartBunch::create(size_t m) {
    /* Each processor constructs m
     * particles. These have to be
     * redistributed when giving
     * the "real" values since they
     * belong to different grids
     * and thus different processors.
     */
    if ( m_particles.size() == 0) {
        m_particles.reserve(15);
        m_particles.resize(m_gdb->finestLevel()+1);
    }
    
    for (size_t i = 0; i < m; ++i) {
        ParticleType p;
        
        p.m_id  = ParticleBase::NextID();
        p.m_cpu = ParallelDescriptor::MyProc();
        
        // fill position with zeros
        for (int j = 0; j < BL_SPACEDIM; ++j)
            p.m_pos[j] = 0.0;
        
        // fill zeros for all floating point attributes
        for (size_t att = 0; att < nAttributes; ++att)
            p.m_data[att] = 0.0;
        
        // assign a particle to a grid --> sets p.m_grid
        ParticleBase::Where(p,m_gdb);
        
        m_particles[p.m_lev][p.m_grid].push_back(p);
    }
    
    buildIndexMapping_m();
}


void AmrPartBunch::gatherStatistics() {
    Inform m("gatherStatistics ");
    Inform m2a("gatherStatistics ", INFORM_ALL_NODES);
    
    double *partPerNode = new double[ParallelDescriptor::NProcs()];
    double *globalPartPerNode = new double[ParallelDescriptor::NProcs()];
    
    for (int i = 0; i < ParallelDescriptor::NProcs(); ++i)
        partPerNode[i] = globalPartPerNode[i] = 0.0;
    
    partPerNode[ParallelDescriptor::MyProc()] = this->getLocalNum();
    
    reduce(partPerNode,
           partPerNode + ParallelDescriptor::NProcs(),
           globalPartPerNode,
           OpAddAssign());
    
    for (int i = 0; i < ParallelDescriptor::NProcs(); ++i)
        m << "Node " << i << " has "
          <<   globalPartPerNode[i]/this->getTotalNum()*100.0 << " \% of the total particles " << endl;
    
    
    for (unsigned int lev = 0; lev < m_particles.size(); ++lev)
        m << "#particles at level " << lev << ": "
          << NumberOfParticlesAtLevel(lev) << endl;
}

void AmrPartBunch::dumpStatistics(const std::string& filename) {
    
    
    std::unique_ptr<double[]> partPerNode(new double[ParallelDescriptor::NProcs()]);
    std::unique_ptr<double[]> globalPartPerNode(new double[ParallelDescriptor::NProcs()]);
    
    for (int i = 0; i < ParallelDescriptor::NProcs(); ++i)
        partPerNode[i] = globalPartPerNode[i] = 0.0;
    
    partPerNode[ParallelDescriptor::MyProc()] = this->getLocalNum();
    
    reduce(partPerNode.get(),
           partPerNode.get() + ParallelDescriptor::NProcs(),
           globalPartPerNode.get(),
           OpAddAssign());
    
    std::size_t total = this->getTotalNum();
    
    if ( Ippl::myNode() == 0 ) {
        std::ofstream out(filename, std::ios::app);
        for (int i = 0; i < ParallelDescriptor::NProcs(); ++i)
            out << globalPartPerNode[i]/double(total)*100.0 << " ";
        out << std::endl;
        out.close();
    }
}

size_t AmrPartBunch::getLocalNum() const {
    return nLocalParticles_m;
}

size_t AmrPartBunch::getTotalNum() const {
    ///@bug In ParticleContainer of BoxLib (check second boolean)
    return TotalNumberOfParticles(true, true);
}


Vector_t AmrPartBunch::getRMin() {
    return Vector_t(0.0, 0.0, 0.0);
}

// 
Vector_t AmrPartBunch::getRMax() {
    return Vector_t(0.0, 0.0, 0.0);
}


Vector_t AmrPartBunch::getHr() {
    return Vector_t(0.0, 0.0, 0.0);
}

double AmrPartBunch::scatter() {
    return 0;
}


void AmrPartBunch::initFields() {}


void AmrPartBunch::gatherCIC() {}


Vector_t AmrPartBunch::getR(int i) {
    int l, g, dq;
    std::tie(l,g,dq) = idxMap_m[i];
    return Vector_t(m_particles[l][g][dq].m_pos[0],
                    m_particles[l][g][dq].m_pos[1],
                    m_particles[l][g][dq].m_pos[2]);
}


double AmrPartBunch::getQM(int i) {
    int l, g, dq;
    std::tie(l,g,dq) = idxMap_m[i];
    return m_particles[l][g][dq].m_data[0];
}

double AmrPartBunch::getMass(int i) {
    int l, g, dq;
    std::tie(l,g,dq) = idxMap_m[i];
    return m_particles[l][g][dq].m_data[4];
}



Vector_t AmrPartBunch::getP(int i) {
    int l, g, dq;
    std::tie(l,g,dq) = idxMap_m[i];
    return Vector_t(m_particles[l][g][dq].m_data[1],
                    m_particles[l][g][dq].m_data[2],
                    m_particles[l][g][dq].m_data[3]);
}


Vector_t AmrPartBunch::getE(int i) {
//     int l, g, dq;
//     std::tie(l,g,dq) = idxMap_m[i];
//     return Vector_t(m_particles[l][g][dq].m_data[4],
//                     m_particles[l][g][dq].m_data[5],
//                     m_particles[l][g][dq].m_data[6]);
    return Vector_t(0.0, 0.0, 0.0);
}
// 
// 
Vector_t AmrPartBunch::getB(int i) {
//     int l, g, dq;
//     std::tie(l,g,dq) = idxMap_m[i];
//     return Vector_t(m_particles[l][g][dq].m_data[7],
//                     m_particles[l][g][dq].m_data[8],
//                     m_particles[l][g][dq].m_data[9]);
    return Vector_t(0.0, 0.0, 0.0);
}


void AmrPartBunch::setR(Vector_t pos, int i) {
    int l, g, dq;
    std::tie(l,g,dq) = idxMap_m[i];
    
    for (int d = 0; d < 3; ++d)
        m_particles[l][g][dq].m_pos[d] = pos(d);
    
    // if going over boundary
    ParticleBase::PeriodicShift(m_particles[l][g][dq], m_gdb);
}


void AmrPartBunch::setQM(double q, int i) {
    int l, g, dq;
    std::tie(l,g,dq) = idxMap_m[i];
    
    m_particles[l][g][dq].m_data[0] = q;
}

void AmrPartBunch::setMass(double m, int i) {
    int l, g, dq;
    std::tie(l,g,dq) = idxMap_m[i];
    
    m_particles[l][g][dq].m_data[4] = m;
}

void AmrPartBunch::setP(Vector_t v, int i) {
    int l, g, dq;
    std::tie(l,g,dq) = idxMap_m[i];
    
    for (int d = 0; d < 3; ++d)
        m_particles[l][g][dq].m_data[d + 1] = v(d);
}


void AmrPartBunch::setE(Vector_t Ef, int i) {
//     int l, g, dq;
//     std::tie(l,g,dq) = idxMap_m[i];
//     
//     for (int d = 0; d < 3; ++d)
//         m_particles[l][g][dq].m_data[d + 4] = Ef(d);
}
// 
// 
void AmrPartBunch::setB(Vector_t Bf, int i) {
//     int l, g, dq;
//     std::tie(l,g,dq) = idxMap_m[i];
//     
//     for (int d = 0; d < 3; ++d)
//         m_particles[l][g][dq].m_data[d + 7] = Bf(d);
}
