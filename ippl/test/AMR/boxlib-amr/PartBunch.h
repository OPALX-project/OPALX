#ifndef PARTBUNCH_H
#define PARTBUNCH_H

/*!
 * @file PartBunch.h
 * @details Particle bunch class
 * for IPPL
 * @authors Matthias Frey \n
 *          Andreas Adelmann \n
 *          Ann Almgren \n
 *          Weiqun Zhang
 * @date LBNL, October 2016
 * @brief Concrete implementation of PartBunchBase. It
 * is the OPAL particle bunch class.
 */


#include "PartBunchBase.h"

/// Particle bunch class using IPPL
template<class PL>
class PartBunch : public IpplParticleBase<PL>,
                  public PartBunchBase
{
private:

    Field<Vektor<double,Dim>,Dim> EFD_m;
    Field<double,Dim> EFDMag_m;

    // Boundary conditions (all periodic)
    BConds<double,Dim,Mesh_t,Center_t> bc_m;
    BConds<Vector_t,Dim,Mesh_t,Center_t> vbc_m;

    Vektor<int,Dim> nr_m;           ///< Number of grid points in each direction (nx, ny, nz)

    bool fieldNotInitialized_m;     ///< Check in case of myUpdate() if field is initialized

    e_dim_tag decomp_m[Dim];

    Vector_t hr_m;                  ///< Mesh spacing
    Vector_t rmin_m;                ///< Lower corner of domain
    Vector_t rmax_m;                ///< Upper corner of domain


    ParticleAttrib<double>     qm; ///< charge-to-mass ratio
    typename PL::ParticlePos_t P;  ///< particle velocity
    typename PL::ParticlePos_t E;  ///< electric field at particle position
    typename PL::ParticlePos_t B;  ///< magnetic field at particle position

public:

    /*!
     * @param pl is the particle layout (defines the distribution mapping)
     * @param hr is the mesh spacing
     * @param rmin is the lower corner of the domain
     * @param rmax is the upper corner of the domain
     */
    PartBunch(PL* pl, Vector_t hr, Vector_t rmin,
              Vector_t rmax, e_dim_tag decomp[Dim]);


    inline const Mesh_t& getMesh() const;

    /*!
     * @returns the grid specifying the domain
     * decomposition
     */
    inline Mesh_t& getMesh();

    /*!
     * @return the layout of the field.
     */
    inline FieldLayout_t& getFieldLayout();

    void savePhaseSpace(std::string fn, int idx);

    inline void gatherCIC();

    // ------------------------------------------------------------------------
    // INHERITED MEMBER FUNCTIONS
    // ------------------------------------------------------------------------
    inline void create(size_t m);

    inline size_t getLocalNum() const;

    inline size_t getTotalNum() const;

    inline Vector_t getR(int i);

    inline double getQM(int i);
    
    inline double getMass(int i);
    
    inline Vector_t getP(int i);

    inline Vector_t getE(int i);

    inline Vector_t getB(int i);

    inline void setR(Vector_t pos, int i);

    inline void setQM(double q, int i);
    
    inline void setMass(double m, int i);
    
    inline void setP(Vector_t v, int i);

    inline void setE(Vector_t Ef, int i);

    inline void setB(Vector_t Bf, int i);

    double scatter();

    void myUpdate();

    void gatherStatistics();
    
    void dumpStatistics(const std::string& filename) {}

    void initFields();

    inline Vector_t getRMin();
    inline Vector_t getRMax();
    inline Vector_t getHr();


    void destroyAll() {

    }

private:

    inline void setBCAllPeriodic_m();

    inline void scatterCIC_m();
};




// ----------------------------------------------------------------------------
// PUBLIC MEMBER FUNCTIONS
// ----------------------------------------------------------------------------


template <class PL>
PartBunch<PL>::PartBunch(PL* pl, Vector_t hr, Vector_t rmin,
                         Vector_t rmax, e_dim_tag decomp[Dim])
    : IpplParticleBase<PL>(pl),
      fieldNotInitialized_m(true),
      hr_m(hr),
      rmin_m(rmin),
      rmax_m(rmax)
{
    // register the particle attributes
    this->addAttribute(qm);
    this->addAttribute(P);
    this->addAttribute(E);
    this->addAttribute(B);

    setBCAllPeriodic_m();

    for (unsigned int i = 0; i < Dim; ++i)
        decomp_m[i] = decomp[i];
}


template <class PL>
const Mesh_t& PartBunch<PL>::getMesh() const {
    return this->getLayout().getLayout().getMesh();
}


template <class PL>
Mesh_t& PartBunch<PL>::getMesh() {
    return this->getLayout().getLayout().getMesh();
}


template <class PL>
FieldLayout_t& PartBunch<PL>::getFieldLayout() {
    return dynamic_cast<FieldLayout_t&>(
        this->getLayout().getLayout().getFieldLayout());
}


template <class PL>
void PartBunch<PL>::savePhaseSpace(std::string fn, int idx) {
    int tag = Ippl::Comm->next_tag(IPPL_APP_TAG4, IPPL_APP_CYCLE);
    std::vector<double> tmp;
    tmp.clear();

    // every node ckecks if he has to dump particles
    for (unsigned i=0; i<this->getLocalNum(); i++) {
        tmp.push_back(this->ID[i]);
        tmp.push_back(this->R[i](0));
        tmp.push_back(this->R[i](1));
        tmp.push_back(this->R[i](2));
        tmp.push_back(this->P[i](0));
        tmp.push_back(this->P[i](1));
        tmp.push_back(this->P[i](2));
    }

    if(Ippl::myNode() == 0) {
        std::ofstream ostr;
        std::string Fn;
        char numbuf[6];
        sprintf(numbuf, "%05d", idx);
        Fn = fn  + std::string(numbuf) + ".dat";
        ostr.open(Fn.c_str(), std::ios::out );
        ostr.precision(15);
        ostr.setf(std::ios::scientific, std::ios::floatfield);

        ostr << " x, px, y, py t, pt, id, node" << std::endl;

        unsigned int dataBlocks=0;
        double x,y,z,px,py,pz,id;
        unsigned  vn;

        for (unsigned i=0; i < tmp.size(); i+=7)
            ostr << tmp[i+1] << " " << tmp[i+4] << " " << tmp[i+2]  << " " \
                    << tmp[i+5] << " " << tmp[i+3] << " " << tmp[i+6]  << " " \
                    << tmp[i]   << " " << Ippl::myNode() << std::endl;

        int notReceived =  Ippl::getNodes() - 1;
        while (notReceived > 0) {
            int node = COMM_ANY_NODE;
            Message* rmsg =  Ippl::Comm->receive_block(node, tag);
            if (rmsg == 0)
                ERRORMSG("Could not receive from client nodes in main." << endl);
            notReceived--;
            rmsg->get(&dataBlocks);
            rmsg->get(&vn);
            for (unsigned i=0; i < dataBlocks; i+=7) {
                rmsg->get(&id);
                rmsg->get(&x);
                rmsg->get(&y);
                rmsg->get(&z);
                rmsg->get(&px);
                rmsg->get(&py);
                rmsg->get(&pz);
                ostr << x  << "\t " << px  << "\t " << y  << "\t " \
                        << py << "\t " << z << "\t " << pz << "\t "   \
                        << id   << "\t " << vn << std::endl;
            }
            delete rmsg;
        }
        ostr.close();
    }
    else {
        unsigned dataBlocks = 0;
        Message* smsg = new Message();
        dataBlocks = tmp.size();
        smsg->put(dataBlocks);
        smsg->put(Ippl::myNode());
        for (unsigned i=0; i < tmp.size(); i++)
            smsg->put(tmp[i]);
        bool res = Ippl::Comm->send(smsg, 0, tag);
        if (! res)
            ERRORMSG("Ippl::Comm->send(smsg, 0, tag) failed " << endl;);
    }
}


template <class PL>
void PartBunch<PL>::gatherCIC() {
    // create interpolater object (cloud in cell method)
    IntCIC myinterp;
    E.gather(EFD_m, this->R, myinterp);
}


// ----------------------------------------------------------------------------
// PRIVATE MEMBER FUNCTIONS
// ----------------------------------------------------------------------------


template <class PL>
inline void PartBunch<PL>::setBCAllPeriodic_m() {
    for (unsigned i = 0; i < 2 * Dim; ++i) {
        this->getBConds()[i] = ParticlePeriodicBCond;
        bc_m[i]  = new PeriodicFace<double  ,Dim,Mesh_t,Center_t>(i);
        vbc_m[i] = new PeriodicFace<Vector_t,Dim,Mesh_t,Center_t>(i);
    }
}


template <class PL>
inline void PartBunch<PL>::scatterCIC_m() {
    // create interpolater object (cloud in cell method)
    IntCIC myinterp;
    qm.scatter(EFDMag_m, this->R, myinterp);
}


// ----------------------------------------------------------------------------
// INHERITED MEMBER FUNCTIONS
// ----------------------------------------------------------------------------


template <class PL>
void PartBunch<PL>::create(size_t m) {
    IpplParticleBase<PL>::create(m);
}


template <class PL>
size_t PartBunch<PL>::getLocalNum() const {
    return IpplParticleBase<PL>::getLocalNum();
}


template <class PL>
size_t PartBunch<PL>::getTotalNum() const {
    return IpplParticleBase<PL>::getTotalNum();
}

template <class PL>
Vector_t PartBunch<PL>::getR(int i) {
    return IpplParticleBase<PL>::R[i];
}


template <class PL>
inline double PartBunch<PL>::getQM(int i) {
    return qm[i];
}

template <class PL>
inline double PartBunch<PL>::getMass(int i) {
    return 0.0;
}


template <class PL>
inline Vector_t PartBunch<PL>::getP(int i) {
    return P[i];
}


template <class PL>
inline Vector_t PartBunch<PL>::getE(int i) {
    return E[i];
}


template <class PL>
inline Vector_t PartBunch<PL>::getB(int i) {
    return B[i];
}

template <class PL>
void PartBunch<PL>::setR(Vector_t pos, int i) {
    for (int d = 0; d < 3; ++d)
        IpplParticleBase<PL>::R[i](d) = pos(d);
}

template <class PL>
void PartBunch<PL>::setQM(double q, int i) {
    qm[i] = q;
}

template <class PL>
void PartBunch<PL>::setMass(double m, int i) {
    //TODO
}

template <class PL>
void PartBunch<PL>::setP(Vector_t v, int i) {
    for (int d = 0; d < 3; ++d)
        P(d) = v(d);
}

template <class PL>
void PartBunch<PL>::setE(Vector_t Ef, int i) {
    for (int d = 0; d < 3; ++d)
        E(d) = Ef(d);
}

template <class PL>
void PartBunch<PL>::setB(Vector_t Bf, int i) {
    for (int d = 0; d < 3; ++d)
        B(d) = Bf(d);
}


template <class PL>
double PartBunch<PL>::scatter() {
    Inform m("scatter ");
    double initialQ = sum(qm);

    scatterCIC_m();

    Field<double,Dim> tmpf;
    NDIndex<Dim> domain = getFieldLayout().getDomain();

    FieldLayout_t *FL  = new FieldLayout_t(getMesh(), decomp_m);
    tmpf.initialize(getMesh(), *FL);

    tmpf[domain] = EFDMag_m[domain];

    NDIndex<Dim> idx = getFieldLayout().getLocalNDIndex();
    NDIndex<Dim> idxdom = getFieldLayout().getDomain();
    NDIndex<Dim> elem;

    double Q = 0.0;
    for (int i=idx[0].min(); i<=idx[0].max(); ++i) {
        elem[0] = Index(i,i);
        for (int j=idx[1].min(); j<=idx[1].max(); ++j) {
            elem[1] = Index(j,j);
            for (int k=idx[2].min(); k<=idx[2].max(); ++k) {
                elem[2] = Index(k,k);
                Q +=  tmpf.localElement(elem);
            }
        }
    }

    reduce(Q,Q,OpAddAssign());
    m << "sum(qm)= " << initialQ << " sum(EFDMag)= " << sum(EFDMag_m) << endl;
    return initialQ-Q;
}


template <class PL>
void PartBunch<PL>::myUpdate() {
    if(fieldNotInitialized_m) {
        fieldNotInitialized_m=false;
        getMesh().set_meshSpacing(&(hr_m[0]));
        getMesh().set_origin(rmin_m);

        EFD_m.initialize(getMesh(), getFieldLayout(), GuardCellSizes<Dim>(1), vbc_m);
        EFDMag_m.initialize(getMesh(), getFieldLayout(), GuardCellSizes<Dim>(1), bc_m);
    }

    this->update();
}


template <class PL>
void PartBunch<PL>::gatherStatistics() {
    Inform m("gatherStatistics ");
    Inform m2a("gatherStatistics ", INFORM_ALL_NODES);

    double *partPerNode = new double[Ippl::getNodes()];
    double *globalPartPerNode = new double[Ippl::getNodes()];

    for (int i=0; i<Ippl::getNodes(); i++)
        partPerNode[i] = globalPartPerNode[i] = 0.0;

    partPerNode[Ippl::myNode()] = this->getLocalNum();

    reduce(partPerNode,partPerNode+Ippl::getNodes(),globalPartPerNode,OpAddAssign());

    for (int i = 0; i < Ippl::getNodes(); ++i)
        m << "Node " << i << " has "
          <<   globalPartPerNode[i]/this->getTotalNum()*100.0 << " \% of the total particles " << endl;
}


template <class PL>
void PartBunch<PL>::initFields() {
    Inform m("initFields ");

    NDIndex<Dim> domain = getFieldLayout().getDomain();

    for (unsigned int i=0; i<Dim; i++)
        nr_m[i] = domain[i].length();

    int nx = nr_m[0];
    int ny = nr_m[1];
    int nz = nr_m[2];

    double phi0 = 0.1*nx;

    m << "rmin= " << rmin_m << " rmax= " << rmax_m << " h= " << hr_m << " n= " << nr_m << endl;

    Index I(nx), J(ny), K(nz);

    assign(EFD_m[I][J][K](0), -2.0*pi*phi0/nx * cos(2.0*pi*(I+0.5)/nx) * cos(4.0*pi*(J+0.5)/ny) * cos(pi*(K+0.5)/nz));

    assign(EFD_m[I][J][K](1),  4.0*pi*phi0/ny * sin(2.0*pi*(I+0.5)/nx) * sin(4.0*pi*(J+0.5)/ny));

    assign(EFD_m[I][J][K](2),  4.0*pi*phi0/ny * sin(2.0*pi*(I+0.5)/nx) * sin(4.0*pi*(J+0.5)/ny));

    assign(EFDMag_m[I][J][K],
            EFD_m[I][J][K](0) * EFD_m[I][J][K](0) +
            EFD_m[I][J][K](1) * EFD_m[I][J][K](1) +
            EFD_m[I][J][K](2) * EFD_m[I][J][K](2));
}


template <class PL>
Vector_t PartBunch<PL>::getRMin() { return rmin_m;}


template <class PL>
Vector_t PartBunch<PL>::getRMax() { return rmax_m;}


template <class PL>
Vector_t PartBunch<PL>::getHr() { return hr_m;}

#endif