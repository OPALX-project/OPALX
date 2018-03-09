
// -*- C++ -*-
/**************************************************************************************************************************************
 *
 * The IPPL Framework
 *
 * This program was prepared by PSI.
 * All rights in the program are reserved by PSI.
 * Neither PSI nor the author(s)
 * makes any warranty, express or implied, or assumes any liability or
 * responsibility for the use of this software
 *

Example:
         grid size
            /|\  interaction radius
           / | \   /  
./PairTest 4 4 4 0.5--commlib mpi --info 9 | tee output

 *************************************************************************************************************************************/

#include "Ippl.h"
#include <string>
#include <vector>
#include "Particle/BoxParticleCachingPolicy.h"
#include "Particle/CellParticleCachingPolicy.h"
#include "Particle/PairBuilder/BasicPairBuilder.h"
#include "Particle/PairBuilder/SortingPairBuilder.h"
#include "Particle/PairBuilder/HashPairBuilder.h"
#include "Particle/PairBuilder/PairConditions.h"

// dimension of our positions
const unsigned Dim = 3;

// some typedefs
typedef UniformCartesian<Dim, double>                               Mesh_t;
typedef ParticleSpatialLayout < double, Dim, Mesh_t,
        BoxParticleCachingPolicy<double, Dim, Mesh_t> >             playout_t;
typedef playout_t::SingleParticlePos_t                              Vector_t;
typedef Cell                                                        Center_t;
typedef CenteredFieldLayout<Dim, Mesh_t, Center_t>                  FieldLayout_t;
typedef Field<double, Dim, Mesh_t, Center_t>                        Field_t;

const double qmmax = 1.0;       // maximum value for particle q/m
const double dt = 1.0;          // size of timestep

template<class PL>
class ChargedParticles : public IpplParticleBase<PL> {
public:
    ParticleAttrib<double>     qm;

    ChargedParticles(PL *pl, Vector_t nr, Vector_t hr, Vector_t rmin, e_dim_tag decomp[Dim], bool gCells = true) :
        IpplParticleBase<PL>(pl),
        nr_m(nr),
        hr_m(hr),
        rmin_m(rmin),
        withGuardCells_m(gCells) {
        this->addAttribute(qm);

        for(int i = 0; i < 2 * Dim; i++) {
            bc_m[i]  = new ParallelPeriodicFace<double, Dim, Mesh_t, Center_t>(i);
            this->getBConds()[i] =  ParticlePeriodicBCond;//ParticleNoBCond;//
        }
        for(int i = 0; i < Dim; i++)
            decomp_m[i] = decomp[i];

        getMesh().set_meshSpacing(&(hr_m[0]));
        getMesh().set_origin(rmin_m);
        if(withGuardCells_m)
            scalF_m.initialize(getMesh(), getFieldLayout(), GuardCellSizes<Dim>(1), bc_m);
        else
            scalF_m.initialize(getMesh(), getFieldLayout(), bc_m);
    }

    inline const Mesh_t &getMesh() const { return this->getLayout().getLayout().getMesh(); }

    inline Mesh_t &getMesh() { return this->getLayout().getLayout().getMesh(); }

    inline const FieldLayout_t &getFieldLayout() const {
        return dynamic_cast<FieldLayout_t &>(this->getLayout().getLayout().getFieldLayout());
    }

    inline FieldLayout_t &getFieldLayout() {
        return dynamic_cast<FieldLayout_t &>(this->getLayout().getLayout().getFieldLayout());
    }

    NDRegion<double, Dim> getLocalRegion() {
        return (*(this->getLayout().getLayout().begin_iv())).second->getDomain();
    }

    bool checkParticles() {
        Inform msg("CheckParticles", INFORM_ALL_NODES);
        bool ok = true;
        int i = 0;
        NDRegion<double, Dim> region = getLocalRegion();
        for(; i < this->getLocalNum(); ++i) {
            NDRegion<double, Dim> ppos;
            for(int d = 0; d < Dim; ++d)
                ppos[d] = PRegion<double>(this->R[i][d], this->R[i][d]);

            if(!region.contains(ppos)) {
                ok = false;
                msg << "Particle misplaced!\nPosition: "
                    << this->R[i] << "\nRegion: " << region << endl;
            }
        }
        for(; i < this->getLocalNum() + this->getGhostNum(); ++i) {
            NDRegion<double, Dim> ppos;
            for(int d = 0; d < Dim; ++d)
                ppos[d] = PRegion<double>(this->R[i][d], this->R[i][d]);

            if(region.contains(ppos)) {
                ok = false;
                msg << "Ghost Particle misplaced!\nPosition: "
                    << this->R[i] << "\nRegion: " << region << endl;
            }
        }
        return ok;
    }

    void printParticles() {
        Inform msg("IpplParticleBase", INFORM_ALL_NODES);
        Ippl::Comm->barrier();

        for(int i = 0; i < Ippl::getNodes(); ++i) {
            if(i == Ippl::myNode()) {
                msg << "local region: " << getLocalRegion() << '\n';
                msg << "local particles:\n";
                int i = 0;
                for(; i < this->getLocalNum(); ++i) {
                    msg << '\t' << this->R[i];
                }
                msg << "\nghost particles\n";
                for(; i < this->getLocalNum() + this->getGhostNum(); ++i) {
                    msg << '\t' << this->R[i];
                }
                msg << endl << endl;
            }
            Ippl::Comm->barrier();
        }
    }

private:
    Field<double, Dim> scalF_m;
    BConds<double, Dim, Mesh_t, Center_t> bc_m;

    Vektor<int, Dim> nr_m;
    Vector_t hr_m;
    Vector_t rmin_m;

    bool withGuardCells_m;
    e_dim_tag decomp_m[Dim];
};

struct PrintPair {
    void operator()(std::size_t i, std::size_t j, ChargedParticles<playout_t> &P) const {
        Inform msg("Pair");
        msg << P.R[i] << ' ' << P.R[j] << endl;
    }
};

struct CountPair {
    CountPair(unsigned &c) : count(c) { count = 0; }
    void operator()(std::size_t i, std::size_t j, ChargedParticles<playout_t> &P) const {
        count++;
    }
    unsigned &count;
};

int main(int argc, char *argv[]) {
    Ippl ippl(argc, argv);
    Inform msg(argv[0]);
    Inform msg2all(argv[0], INFORM_ALL_NODES);

    IpplTimings::TimerRef allTimer = IpplTimings::getTimer("AllTimer");
    IpplTimings::startTimer(allTimer);

    Vektor<int, Dim> nr;

    unsigned param = 1;

    if(Dim == 3) {
        nr = Vektor<int, Dim>(atoi(argv[param++]), atoi(argv[param++]), atoi(argv[param++]));
    } else {
        nr = Vektor<int, Dim>(atoi(argv[param++]), atoi(argv[param++]));

    }

    double interaction_radius = atof(argv[param++]);

    e_dim_tag decomp[Dim];
    Mesh_t *mesh;
    FieldLayout_t *FL;
    ChargedParticles<playout_t>  *P;

    NDIndex<Dim> domain;
    for(int i = 0; i < Dim; i++)
        domain[i] = domain[i] = Index(nr[i]);

    for(int d = 0; d < Dim; ++d)
        decomp[d] = PARALLEL;

    // create mesh and layout objects for this problem domain
    mesh          = new Mesh_t(domain);
    FL            = new FieldLayout_t(*mesh, decomp);
    playout_t *PL = new playout_t(*FL, *mesh);

    Vector_t hr(1.0);
    Vector_t rmin(0.0);

    PL->setAllCacheDimensions(interaction_radius);
    PL->enableCaching();
    //PL->setAllCacheDimensions(0);

    //PL->setAllCacheCellRanges(2);

    P = new ChargedParticles<playout_t>(PL, nr, hr, rmin, decomp, true);
    INFOMSG(P->getMesh() << endl);
    INFOMSG(P->getFieldLayout() << endl);
    msg << endl << endl;

    Ippl::Comm->barrier();

    {
        size_t size = size_t(nr[0] - 1) * size_t(nr[1] - 1) * size_t(nr[2] - 1);
        size_t begin = size * Ippl::myNode() / Ippl::getNodes();
        size_t end = size * (Ippl::myNode() + 1) / Ippl::getNodes();
        P->create(end - begin);
        size_t j = 0;
        for(size_t i = begin; i < end; ++i, ++j) {
            P->R[j](0) = double(0.5 + i % size_t(nr[0] - 1));
            P->R[j](1) = double(0.5 + (i / size_t(nr[0] - 1)) % size_t(nr[1] - 1));
            P->R[j](2) = double(0.5 + (i / size_t(nr[0] - 1)) / size_t(nr[1] - 1));
        }
    }

    P->qm = 1.0;
    P->update();

    unsigned basic_count;
    IpplTimings::TimerRef basicPairTimer = IpplTimings::getTimer("BasicPairTimer");
    IpplTimings::startTimer(basicPairTimer);
    BasicPairBuilder< ChargedParticles<playout_t> > BPB(*P);
    BPB.for_each(RadiusCondition<double, Dim>(interaction_radius), CountPair(basic_count));
    IpplTimings::stopTimer(basicPairTimer);
    msg << "BasicPairBuilder found " << basic_count << " pairs" << endl;

    unsigned sort_count;
    IpplTimings::TimerRef sortPairTimer = IpplTimings::getTimer("SortPairTimer");
    IpplTimings::startTimer(sortPairTimer);
    SortingPairBuilder< ChargedParticles<playout_t> > SPB(*P);
    SPB.for_each(RadiusCondition<double, Dim>(interaction_radius), CountPair(sort_count));
    IpplTimings::stopTimer(sortPairTimer);
    msg << "SortingPairBuilder found " << sort_count << " pairs" << endl;

    unsigned hash_count;
    IpplTimings::TimerRef hashPairTimer = IpplTimings::getTimer("HashPairTimer");
    IpplTimings::startTimer(hashPairTimer);
    HashPairBuilder< ChargedParticles<playout_t> > HPB(*P);
    HPB.for_each(RadiusCondition<double, Dim>(interaction_radius), CountPair(hash_count));
    IpplTimings::stopTimer(hashPairTimer);
    msg << "HashPairBuilder found " << hash_count << " pairs" << endl;

    IpplTimings::stopTimer(allTimer);
    IpplTimings::print();

    delete P;
    delete FL;
    delete mesh;
    return 0;
}

