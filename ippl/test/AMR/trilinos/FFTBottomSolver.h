#ifndef FFT_SOLVER_H
#define FFT_SOLVER_H

#include "BottomSolver.h"

#include "Solvers/FFTPoissonSolver.h"

#include <map>

#include "../AmrOpal.h"

#include "Amr/AmrDefs.h"

template <class Level>
class FFTBottomSolver : public BottomSolver<Teuchos::RCP<amr::matrix_t>,
                                            Teuchos::RCP<amr::multivector_t>,
                                            Level >,
                        public FFTPoissonSolver
{
public:
    typedef amr::matrix_t matrix_t;
    typedef amr::vector_t vector_t;
    typedef amr::scalar_t scalar_t;
    typedef amr::multivector_t mv_t;
    typedef amr::operator_t op_t;
    typedef amr::local_ordinal_t lo_t;
    typedef amr::global_ordinal_t go_t;
    typedef amr::AmrIntVect_t AmrIntVect_t;
    typedef amr::AmrProcMap_t AmrProcMap_t;
    typedef amr::AmrGrid_t AmrGrid_t;
    typedef amr::AmrGeometry_t AmrGeometry_t;
    typedef amr::AmrBox_t AmrBox_t;
    
public:
    
    FFTBottomSolver(Mesh_t *mesh,
                    FieldLayout_t *fl,
                    std::string greensFunction);
    
    static Mesh_t* initMesh(AmrOpal* amrobject_p);
    
    static FieldLayout_t* initFieldLayout(Mesh_t *mesh,
                                          AmrOpal* amrobject_p);
    
    void solve(const Teuchos::RCP<mv_t>& x,
               const Teuchos::RCP<mv_t>& b);
    
    void setOperator(const Teuchos::RCP<matrix_t>& A,
                      Level* level_p);
    
    std::size_t getNumIters() { return 1; }
    
private:
    void fillMap_m(Level* level_p);
    
    void field2vector_m(const Teuchos::RCP<mv_t>& vector);
    
    void vector2field_m(const Teuchos::RCP<mv_t>& vector);
    
private:
    std::map<lo_t, AmrIntVect_t> map_m;
    
    Field_t rho_m;
    BConds<double, 3, Mesh_t, Center_t> bc_m;
    
    Level* level_mp;
};

#endif
