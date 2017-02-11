#ifndef PARTBUNCHAMR_H
#define PARTBUNCHAMR_H

#include "AmrParticleBase.h"

template<class PLayout>
class PartBunchAmr : public AmrParticleBase<PLayout>
{

public:
    
    ParticleAttrib<double> qm;
    typename PLayout::ParticlePos_t E;
    typename PLayout::ParticlePos_t P;
    ParticleAttrib<double> mass;
  ParticleAttrib<Vektor<double,2> > Rphase; //velocity of the particles
    
    PartBunchAmr() {
        this->addAttribute(qm);
        this->addAttribute(E);
        this->addAttribute(P);
        this->addAttribute(mass);
	this->addAttribute(Rphase);
    }
    
    ~PartBunchAmr() {}
    
    
    
    
    void gatherStatistics() {
        Inform msg("gatherStatistics ");
        
        std::unique_ptr<double[]> partPerNode( new double[Ippl::getNodes()] );
        std::unique_ptr<double[]> globalPartPerNode( new double[Ippl::getNodes()] );
    
        for (int i = 0; i < Ippl::getNodes(); ++i)
            partPerNode[i] = globalPartPerNode[i] = 0.0;
        
        partPerNode[Ippl::myNode()] = this->getLocalNum();
        
        reduce(partPerNode.get(),
            partPerNode.get() + Ippl::getNodes(),
            globalPartPerNode.get(),
            OpAddAssign());
        
        for (int i = 0; i < Ippl::getNodes(); ++i)
            msg << "Node " << i << " has "
                <<   globalPartPerNode[i] / this->getTotalNum() * 100.0
                << " \% of the total particles " << endl;
    }
    
    void dumpStatistics(const std::string& filename) {
        std::unique_ptr<double[]> partPerNode(new double[Ippl::getNodes()]);
        std::unique_ptr<double[]> globalPartPerNode(new double[Ippl::getNodes()]);
        
        for (int i = 0; i < Ippl::getNodes(); ++i)
            partPerNode[i] = globalPartPerNode[i] = 0.0;
        
        partPerNode[Ippl::myNode()] = this->getLocalNum();
        
        reduce(partPerNode.get(),
            partPerNode.get() + Ippl::getNodes(),
            globalPartPerNode.get(),
            OpAddAssign());
        
        std::size_t total = this->getTotalNum();
        
        if ( Ippl::myNode() == 0 ) {
            std::ofstream out(filename, std::ios::app);
            for (int i = 0; i < Ippl::getNodes(); ++i)
                out << globalPartPerNode[i] / double(total)*100.0 << " ";
            out << std::endl;
            out.close();
        }
    }
    
    void python_format(int step) {
        std::string st = std::to_string(step);
        Inform::WriteMode wm = Inform::OVERWRITE;
        for (int i = 0; i < Ippl::getNodes(); ++i) {
            if ( i == Ippl::myNode() ) {
                wm = (i == 0) ? wm : Inform::APPEND;
                
                PLayout* playout = &this->getLayout();
                const ParGDBBase* gdb = playout->GetParGDB();
                
                std::string grid_file = "pyplot_grids_" + st + ".dat";
                Inform msg("", grid_file.c_str(), wm, i);
                for (int l = 0; l < gdb->finestLevel() + 1; ++l) {
                    Geometry geom = gdb->Geom(l);
                    for (int g = 0; g < gdb->ParticleBoxArray(l).size(); ++g) {
                        msg << l << ' ';
                        RealBox loc = RealBox(gdb->boxArray(l)[g],geom.CellSize(),geom.ProbLo());
                        for (int n = 0; n < BL_SPACEDIM; n++)
                            msg << loc.lo(n) << ' ' << loc.hi(n) << ' ';
                        msg << endl;
                    }
                }
                
                std::string particle_file = "pyplot_particles_" + st + ".dat";
                Inform msg2all("", particle_file.c_str(), wm, i);
                for (size_t i = 0; i < this->getLocalNum(); ++i) {
                    msg2all << this->R[i](0) << " "
                            << this->R[i](1) << " "
                            << this->R[i](2) << " "
                            << this->P[i](0) << " "
                            << this->P[i](1) << " "
                            << this->P[i](2)
                            << endl;
                }
            }
            Ippl::Comm->barrier();
        }
    }

};


#endif
