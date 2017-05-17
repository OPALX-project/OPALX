
#include "Ippl.h"

#include "AmrParticle/AmrParticleBase.h"
#include "AmrParticle/ParticleAmrLayout.h"

#include <memory>


template <class T, unsigned Dim >
class BoxLibLayout : public ParticleAmrLayout<T, Dim>
{
public:
    typedef typename ParticleAmrLayout<T, Dim>::pair_t pair_t;
    typedef typename ParticleAmrLayout<T, Dim>::pair_iterator pair_iterator;
    typedef typename ParticleAmrLayout<T, Dim>::SingleParticlePos_t SingleParticlePos_t;
    typedef typename ParticleAmrLayout<T, Dim>::Index_t Index_t;
    
    typedef double AmrField_t;
    typedef double AmrFieldContainer_t;
    typedef typename ParticleAmrLayout<T, Dim>::ParticlePos_t ParticlePos_t;
    typedef ParticleAttrib<Index_t> ParticleIndex_t;
    
    BoxLibLayout() { }
    
    
    void update(AmrParticleBase< BoxLibLayout<T,Dim> >& PData, int lev_min = 0,
                const ParticleAttrib<char>* canSwap = 0)
    {
        std::cout << "BoxLibLayout::update()" << std::endl;
    }
    
    void update(IpplParticleBase< BoxLibLayout<T,Dim> >& PData, 
                const ParticleAttrib<char>* canSwap=0)
    {
        std::cout << "IpplBase update" << std::endl;
        //TODO: exit since we need AmrParticleBase with grids and levels for particles for this layout
        //if IpplParticleBase is used something went wrong
    }
};

template <class PLayout>
class BoxLibParticle : public AmrParticleBase<PLayout>
{
public:
    BoxLibParticle() {
        this->initializeAmr();
    }
    
};

class AbstractBunch {
    
public:
    typedef ParticleAttrib<ParticleLayout<double, 3>::SingleParticlePos_t> ParticlePos_t;
    typedef ParticleAttrib<ParticleLayout<double, 3>::Index_t>    ParticleIndex_t;
    
    AbstractBunch(ParticlePos_t &R, ParticleIndex_t& ID) : RR(R), IID(ID) {}
    
    virtual void print() = 0;
    
    virtual void createP(std::size_t n) = 0;
    
    ParticlePos_t& RR;
    ParticleIndex_t& IID;
    ParticleAttrib< int > Bin;
};

template <class Bunch>
class PartBunchBase : public AbstractBunch
{
public:
    
    PartBunchBase(ParticlePos_t& R, ParticleIndex_t& ID) : AbstractBunch(R, ID),
                                                           bunch_mp(static_cast<Bunch*>(this)) {
    }
    
    void createP(std::size_t n) {
        bunch_mp->create(n);
        
        srand(42);
        for (std::size_t i = 0; i < n; ++i) {
            bunch_mp->R[i][0] = (double)rand() / RAND_MAX;
            bunch_mp->R[i][1] = (double)rand() / RAND_MAX;
            bunch_mp->R[i][2] = (double)rand() / RAND_MAX;
            bunch_mp->Bin[i] = i;
        }
        bunch_mp->update();
    }
    
    
private:
    Bunch* bunch_mp;
};

class PartBunch : public PartBunchBase<PartBunch>,
                  public IpplParticleBase<ParticleSpatialLayout<double, 3> >
{
public:
    
    
    PartBunch() : PartBunchBase(this->R, this->ID) {
        this->initialize(new ParticleSpatialLayout<double, 3>()); // TODO Where is this done?
        this->addAttribute(Bin);
    }
    
    void print() {
        std::cout << "LocalNum:  " << this->getLocalNum() << std::endl
                  << "TotalNum:  " << this->getTotalNum() << std::endl
                  << "R[0] =     " << R[0] << std::endl
                  << "ID[0] =    " << ID[0] << std::endl
                  << "Bin[1] =   " << Bin[1] << std::endl;
    }
    
};

class AmrPartBunch : public PartBunchBase<AmrPartBunch>,
                     public BoxLibParticle<BoxLibLayout<double, 3> >
{
public:
    
//     ParticleAttrib< int > Bin;
    
    AmrPartBunch() : PartBunchBase(this->R, this->ID) {
        this->initialize(new BoxLibLayout<double, 3>()); // TODO Where is this done?
        this->addAttribute(Bin);
    }
    
    void print() {
        std::cout << "LocalNum:  " << this->getLocalNum() << std::endl
                  << "TotalNum:  " << this->getTotalNum() << std::endl
                  << "R[0] =     " << R[0] << std::endl
                  << "ID[0] =    " << ID[0] << std::endl
                  << "Level[0] = " << level[0] << std::endl
                  << "Grid[0]  = " << grid[0] << std::endl
                  << "Bin[1] =   " << Bin[1] << std::endl;
    }
    
};


int main(int argc, char** argv) {
    
    Ippl ippl(argc, argv);
    
    std::cout << "Init PartBunch" << std::endl
              << "--------------" << std::endl;
    
    std::unique_ptr<AbstractBunch> bunch(new PartBunch());
    
    
    bunch->createP(10);
    
    bunch->print();
    
    std::cout << std::endl << "Init AmrPartBunch" << std::endl
              << "-----------------" << std::endl;
    
    bunch = std::unique_ptr<AbstractBunch>(new AmrPartBunch());
    
    bunch->createP(10);
    
    bunch->print();
    
    std::cout << std::endl << "Modify Bin[1]" << std::endl
              << "-------------" << std::endl;
    
    bunch->Bin[1] = 100;
    
    bunch->print();
    
    std::cout << std::endl << "Modify RR[0]" << std::endl
              << "-------------" << std::endl;
    
    bunch->RR[0] = {1.0, 2.0, 3.0};
    
    bunch->print();
    
    return 0;
}