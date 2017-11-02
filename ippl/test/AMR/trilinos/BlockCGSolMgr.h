#ifndef TRILINOS_SOLVER_H
#define TRILINOS_SOLVER_H

#include "LinearSolver.h"

#include "AmrMultiGridCore.h"

#include <BelosTpetraAdapter.hpp>

#include <BelosLinearProblem.hpp>
#include <BelosEpetraAdapter.hpp>

#include <BelosBlockCGSolMgr.hpp>

class BlockCGSolMgr : public LinearSolver<Teuchos::RCP<amr::matrix_t>,
                                           Teuchos::RCP<amr::vector_t> >
{
public:
    typedef amr::matrix_t matrix_t;
    typedef amr::vector_t vector_t;
    typedef amr::scalar_t scalar_t;
    typedef amr::multivector_t mv_t;
    typedef amr::operator_t op_t;
    
    typedef Belos::BlockCGSolMgr<scalar_t, mv_t, op_t> solver_t;
    typedef Belos::LinearProblem<scalar_t, mv_t, op_t> problem_t;
    
public:
    
    BlockCGSolMgr() : problem_mp(Teuchos::null) {
        params_mp = Teuchos::rcp( new Teuchos::ParameterList );
    //     params->set( "Block Size", 4 );
        params_mp->set("Maximum Iterations", 1000);
        params_mp->set("Convergence Tolerance", 1.0e-1);
        params_mp->set("Block Size", 32);
        
        solver_mp = Teuchos::rcp( new solver_t() );
        solver_mp->setParameters(params_mp);
    }
    
    ~BlockCGSolMgr() {
        problem_mp = Teuchos::null;
        params_mp = Teuchos::null;
        solver_mp = Teuchos::null;
    };
    
    void solve(const Teuchos::RCP<matrix_t>& A,
               Teuchos::RCP<vector_t>& x,
               const Teuchos::RCP<vector_t>& b, double tol)
    {
        /*
         * solve linear system Ax = b
         */
        problem_mp = Teuchos::rcp( new problem_t(A, x, b) );
        
        
        if ( !problem_mp->setProblem() )
            throw std::runtime_error("Belos::LinearProblem failed to set up correctly!");
        
        solver_mp->setProblem(problem_mp);
        
        params_mp->set("Convergence Tolerance", tol);
        solver_mp->setParameters(params_mp);
        
        Belos::ReturnType ret = solver_mp->solve();
        
        // get the solution from the problem
        if ( ret == Belos::Converged ) {
            x->assign(*problem_mp->getLHS());
            
            // print stuff
//             if ( epetra_comm_m.MyPID() == 0 ) {
//                 std::cout << "Achieved tolerance: " << solver_mp->achievedTol() << std::endl
//                         << "Number of iterations: " << solver_mp->getNumIters() << std::endl;
//             }
            
        } else {
//             if ( epetra_comm_m.MyPID() == 0 ) {
                std::cout << "Not converged. Achieved tolerance after " << solver_mp->getNumIters() << " iterations is "
                        << solver_mp->achievedTol() << "." << std::endl;
//             }
        }
        
        
        problem_mp = Teuchos::null;
    }
    
    
    void residual(Teuchos::RCP<vector_t>& r, const Teuchos::RCP<vector_t>& x, const Teuchos::RCP<vector_t>& b) {
        problem_mp->computeCurrResVec(r.get(), x.get(), b.get());
    }
    
private:
    Teuchos::RCP<problem_t> problem_mp;
    Teuchos::RCP<Teuchos::ParameterList> params_mp;
    Teuchos::RCP<solver_t>  solver_mp;
};

#endif
