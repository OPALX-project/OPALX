#ifndef BELOS_SOLVER_H
#define BELOS_SOLVER_H

#include "AmrMultiGridCore.h"

#include <Amesos2.hpp>
#include <Amesos2_TpetraMultiVecAdapter.hpp>

#include <string>

class AmesosBottomSolver : public BottomSolver<Teuchos::RCP<amr::matrix_t>,
                                               Teuchos::RCP<amr::vector_t> >
{
public:
    typedef amr::matrix_t matrix_t;
    typedef amr::vector_t vector_t;
//     typedef amr::multivector_t mv_t;
//     typedef amr::operator_t op_t;
    
    typedef Amesos2::Solver<matrix_t, vector_t> solver_t;
    
public:
    
    AmesosBottomSolver(std::string solvertype = "klu2");
    
    ~AmesosBottomSolver();
    
    void solve(const Teuchos::RCP<vector_t>& x,
               const Teuchos::RCP<vector_t>& b);
    
    void setOperator(const Teuchos::RCP<matrix_t>& A);
    
private:
    
    std::string solvertype_m;
    
    Teuchos::RCP<solver_t>  solver_mp;
};

#endif
