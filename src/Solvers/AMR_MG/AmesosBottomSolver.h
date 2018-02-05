#ifndef AMESOS_SOLVER_H
#define AMESOS_SOLVER_H

#include "AmrMultiGridCore.h"

#include <Amesos2.hpp>

#include <string>

class AmesosBottomSolver : public BottomSolver<Teuchos::RCP<amr::matrix_t>,
                                               Teuchos::RCP<amr::multivector_t> >
{
public:
    typedef amr::matrix_t matrix_t;
    typedef amr::multivector_t mv_t;
    
    typedef Amesos2::Solver<matrix_t, mv_t> solver_t;
    
public:
    
    AmesosBottomSolver(std::string solvertype = "klu2");
    
    ~AmesosBottomSolver();
    
    void solve(const Teuchos::RCP<mv_t>& x,
               const Teuchos::RCP<mv_t>& b);
    
    void setOperator(const Teuchos::RCP<matrix_t>& A);
    
private:
    
    std::string solvertype_m;
    
    Teuchos::RCP<solver_t> solver_mp;
};

#endif
