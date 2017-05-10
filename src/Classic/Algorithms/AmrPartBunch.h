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
