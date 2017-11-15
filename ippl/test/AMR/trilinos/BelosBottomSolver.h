#ifndef BELOS_SOLVER_H
#define BELOS_SOLVER_H

#include "AmrMultiGridCore.h"

#include <BelosLinearProblem.hpp>
#include <BelosSolverFactory_Tpetra.hpp>

#include <Ifpack2_Factory.hpp>

#include <string>

class BelosBottomSolver : public BottomSolver<Teuchos::RCP<amr::matrix_t>,
                                              Teuchos::RCP<amr::vector_t> >
{
public:
    typedef amr::matrix_t matrix_t;
    typedef amr::vector_t vector_t;
    typedef amr::scalar_t scalar_t;
    typedef amr::multivector_t mv_t;
    typedef amr::operator_t op_t;
    typedef amr::local_ordinal_t lo_t;
    typedef amr::global_ordinal_t go_t;
    typedef amr::node_t node_t;
    
    typedef Belos::SolverManager<scalar_t, mv_t, op_t> solver_t;
    typedef Belos::LinearProblem<scalar_t, mv_t, op_t> problem_t;
    
    typedef Ifpack2::Preconditioner<scalar_t, lo_t, go_t, node_t> precond_t;
    
public:
    
    BelosBottomSolver(std::string solvertype = "Pseudoblock CG",
                      Preconditioner precond = Preconditioner::NONE);
    
    ~BelosBottomSolver();
    
    void solve(const Teuchos::RCP<vector_t>& x,
               const Teuchos::RCP<vector_t>& b);
    
    void setOperator(const Teuchos::RCP<matrix_t>& A);
    
private:
    
    void initSolver_m(std::string solvertype);
    
    void initPreconditioner_m();
    
    void computePrecond_m(const Teuchos::RCP<const matrix_t>& A);
    
private:
    Teuchos::RCP<problem_t> problem_mp;
    Teuchos::RCP<Teuchos::ParameterList> params_mp;
    Teuchos::RCP<solver_t>  solver_mp;
    
    Preconditioner precond_m;
    std::string prectype_m;
    Teuchos::RCP<Teuchos::ParameterList> prec_mp;
    Teuchos::RCP<precond_t> P_mp;
    
    scalar_t reltol_m;
    
    int maxiter_m;
};

#endif
