#ifndef AMR_PART_BUNCH_H
#define AMR_PART_BUNCH_H

#include "Algorithms/PartBunchBase.h"
#include "Amr/AmrObject.h"

class AmrPartBunch : public PartBunchBase<double, 3>
{
public:
    typedef AmrParticle_t pbase_t;
    
public:
    
    AmrPartBunch(const PartData *ref);

    /// Conversion.
    AmrPartBunch(const std::vector<OpalParticle> &,
                 const PartData *ref);

    AmrPartBunch(const AmrPartBunch &);
    
    ~AmrPartBunch();
    
    pbase_t *getAmrParticleBase();
    
    const pbase_t *getAmrParticleBase() const;
    
//     pbase_t* clone();
    
    void initialize(FieldLayout_t *fLayout);
    
    VectorPair_t getEExtrema();
    
    double getRho(int x, int y, int z);
    
//     const Mesh_t &getMesh() const;
    
//     void setMesh(Mesh_t *mesh);
//     Mesh_t &getMesh();
    
//     void setFieldLayout(FieldLayout_t* fLayout);
    
    FieldLayout_t &getFieldLayout();
    
    
    void boundp();
    
    void computeSelfFields();
    
    void computeSelfFields(int bin);
    
    void computeSelfFields_cycl(double gamma);
    
    void computeSelfFields_cycl(int bin);
    
    void setSolver(FieldSolver *fs) {
        PartBunchBase<double, 3>::setSolver(fs);
        this->amrobj_mp = fs->getAmrObject();
    }
    
    // AmrPartBunch only
    PoissonSolver *getFieldSolver() {
        return fs_m->solver_m;
    }
    
    const PoissonSolver *getFieldSolver() const {
        return fs_m->solver_m;
    }
    
    void setBaseLevelMeshSpacing(const Vector_t& hr) {
        for (int i = 0; i < 3; ++i)
            hr_m[i] = hr[i];
    }
    
    
    //FIXME BCs
    void setBCAllPeriodic() {}
    void setBCAllOpen() {}
    void setBCForDCBeam() {}
    
    
    //BEGIN REMOVE
    void python_format(int step) {
        std::string st = std::to_string(step);
        Inform::WriteMode wm = Inform::OVERWRITE;
        for (int i = 0; i < Ippl::getNodes(); ++i) {
            if ( i == Ippl::myNode() ) {
                wm = (i == 0) ? wm : Inform::APPEND;
                
                AmrLayout_t* playout = static_cast<AmrLayout_t*>(&this->getLayout());
                
                std::string grid_file = "pyplot_grids_" + st + ".dat";
                Inform msg("", grid_file.c_str(), wm, i);
                for (int l = 0; l < playout->finestLevel() + 1; ++l) {
                    amrex::Geometry geom = playout->Geom(l);
                    for (int g = 0; g < playout->ParticleBoxArray(l).size(); ++g) {
                        msg << l << ' ';
                        amrex::RealBox loc = amrex::RealBox(playout->boxArray(l)[g],geom.CellSize(),geom.ProbLo());
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
    
    //END
    
    
private:
    void updateFieldContainers_m();
    
    void updateDomainLength(Vektor<int, 3>& grid);
    
    void updateFields(const Vector_t& hr, const Vector_t& origin);
    
private:
    
    /* pointer to AMR object that is part
     * of solver_m (AmrPoissonSolver) in src/Structure/FieldSolver.h
     */
    AmrObject *amrobj_mp;
    pbase_t *amrpbase_mp;
    
    /* We need this due to H5PartWrapper etc, but it's always nullptr.
     * Thus, don't use it.
     */
    FieldLayout_t* fieldlayout_m;
};

#endif
