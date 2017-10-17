#ifndef TRILINOS_SOLVER_H
#define TRILINOS_SOLVER_H

#include "LinearSolver.h"

#include "AmrMultiGridCore.h"

#include <Epetra_LinearProblem.h>

#include <BelosLinearProblem.hpp>
#include <BelosEpetraAdapter.hpp>

#include <BelosBlockCGSolMgr.hpp>

class BlockCGSolMgr : public LinearSolver<Teuchos::RCP<amr::matrix_t>,
                                           Teuchos::RCP<amr::vector_t> >
{
public:
    typedef amr::matrix_t matrix_t;
    typedef amr::vector_t vector_t;
    
public:
    
    BlockCGSolMgr() : problem(Teuchos::null), params(Teuchos::null) { }
    
    ~BlockCGSolMgr() {
        problem = Teuchos::null;
        params = Teuchos::null;
    };
    
    void solve(const Teuchos::RCP<matrix_t>& A,
               Teuchos::RCP<vector_t>& x,
               const Teuchos::RCP<vector_t>& b)
    {
        /*
         * solve linear system Ax = b
         */
        problem = Teuchos::rcp( new Belos::LinearProblem<double, Epetra_MultiVector, Epetra_Operator>(A, x, b) );
        
        
        if ( !problem->setProblem() )
            throw std::runtime_error("Belos::LinearProblem failed to set up correctly!");
        
        params = Teuchos::rcp( new Teuchos::ParameterList );
        
    //     params->set( "Block Size", 4 );
        params->set("Maximum Iterations", 1000);
        params->set("Convergence Tolerance", 1.0e-9);
        params->set("Block Size", 32);
        
        Belos::BlockCGSolMgr<double, Epetra_MultiVector, Epetra_Operator> solver(problem, params);
        
        Belos::ReturnType ret = solver.solve();
        
        // get the solution from the problem
        if ( ret == Belos::Converged ) {
            x = Teuchos::null;
            x = Teuchos::rcp( new vector_t(Epetra_DataAccess::Copy, *problem->getLHS(), 0) );
            
            // print stuff
//             if ( epetra_comm_m.MyPID() == 0 ) {
                std::cout << "Achieved tolerance: " << solver.achievedTol() << std::endl
                        << "Number of iterations: " << solver.getNumIters() << std::endl;
//             }
            
        } else {
//             if ( epetra_comm_m.MyPID() == 0 ) {
                std::cout << "Not converged. Achieved tolerance after " << solver.getNumIters() << " iterations is "
                        << solver.achievedTol() << "." << std::endl;
//             }
        }
    }
    
    
    void residual(Teuchos::RCP<vector_t>& r, const Teuchos::RCP<vector_t>& x, const Teuchos::RCP<vector_t>& b) {
        problem->computeCurrResVec(r.get(), x.get(), b.get());
    }
    
private:
    Teuchos::RCP<Belos::LinearProblem<double, Epetra_MultiVector, Epetra_Operator> > problem;
    Teuchos::RCP<Teuchos::ParameterList> params;
};

#endif
