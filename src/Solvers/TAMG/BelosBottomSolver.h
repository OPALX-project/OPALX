#ifndef BELOS_SOLVER_H
#define BELOS_SOLVER_H

#include "AmrMultiGridCore.h"

#include <BelosLinearProblem.hpp>
#include <BelosTpetraAdapter.hpp>
#include <BelosSolverFactory.hpp>
// #include <BelosSolverFactory_Tpetra.hpp>

#include <Ifpack2_Factory.hpp>

#include <string>

/// Interface to Belos solvers of the Trilinos package
class BelosBottomSolver : public BottomSolver<Teuchos::RCP<amr::matrix_t>,
                                              Teuchos::RCP<amr::multivector_t> >
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
    typedef amr::Preconditioner Preconditioner;
    
    typedef Belos::SolverManager<scalar_t, mv_t, op_t> solver_t;
    typedef Belos::LinearProblem<scalar_t, mv_t, op_t> problem_t;
    
    typedef Ifpack2::Preconditioner<scalar_t, lo_t, go_t, node_t> precond_t;
    
public:
    
    /*!
     * @param solvertype to use
     * @param precond preconditioner of matrix
     */
    BelosBottomSolver(std::string solvertype = "Pseudoblock CG",
                      Preconditioner precond = Preconditioner::NONE);
    
    ~BelosBottomSolver();
    
    void solve(const Teuchos::RCP<mv_t>& x,
               const Teuchos::RCP<mv_t>& b);
    
    void setOperator(const Teuchos::RCP<matrix_t>& A);
    
private:
    /*!
     * Create a solver instance
     * @param solvertype to create
     */
    void initSolver_m(std::string solvertype);
    
    /*!
     * Initialize preconditioner parameter list
     */
    void initPreconditioner_m();
    
    /*!
     * Create the preconditioner
     * @param A system matrix
     */
    void computePrecond_m(const Teuchos::RCP<const matrix_t>& A);
    
private:
    Teuchos::RCP<problem_t> problem_mp;             ///< represents linear problem
    Teuchos::RCP<Teuchos::ParameterList> params_mp; ///< parameter list of solver
    Teuchos::RCP<solver_t>  solver_mp;              ///< solver instance
    
    Preconditioner precond_m;                       ///< preconditioner type
    std::string prectype_m;                         ///< preconditioner type
    Teuchos::RCP<Teuchos::ParameterList> prec_mp;   ///< parameter list of preconditioner
    Teuchos::RCP<precond_t> P_mp;                   ///< preconditioner
    
    scalar_t reltol_m;                              ///< relative tolerance
    
    /// allowed number of steps for iterative solvers
    int maxiter_m;
};

#endif
