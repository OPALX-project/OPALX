#include "AmrSmoother.h"
#include <string>

AmrSmoother::AmrSmoother(const Teuchos::RCP<const matrix_t>& A,
                         const Smoother& smoother,
                         local_ordinal_t nSweeps)
{
    const std::string type = "RELAXATION";
    
    Ifpack2::Factory factory;
    prec_mp = factory.create(type, A);
    
    params_mp = Teuchos::rcp( new Teuchos::ParameterList );
    
    this->initParameter_m(smoother, nSweeps);
    
    
    prec_mp->setParameters(*params_mp);
    prec_mp->initialize();
    prec_mp->compute();
}


void AmrSmoother::smooth(const Teuchos::RCP<vector_t>& x,
                         const Teuchos::RCP<matrix_t>& A,
                         const Teuchos::RCP<vector_t>& b)
{
    Teuchos::RCP<vector_t> r = Teuchos::rcp( new vector_t(A->getDomainMap(), false) );
    A->apply(*x, *r);
    r->update(1.0, *x, -1.0);
    
    prec_mp->apply(*r, *x, Teuchos::NO_TRANS,
                   Teuchos::ScalarTraits<scalar_t>::one(),
                   Teuchos::ScalarTraits<scalar_t>::zero());
}


void AmrSmoother::initParameter_m(const Smoother& smoother,
                                  local_ordinal_t nSweeps)
{
    if ( params_mp == Teuchos::null )
        params_mp = Teuchos::rcp( new Teuchos::ParameterList );
    
    
    std::string type = "";
    double damping = 1.0;
    bool l1 = true;
    
    switch ( smoother ) {
        
        case GAUSS_SEIDEL:
        {
            type = "Gauss-Seidel";
            break;
        }
        case JACOBI:
        {
            type = "Jacobi";
            break;
        }
        default:
            break;
    };
    
    params_mp->set("relaxation: type", type);
    params_mp->set("relaxation: sweeps", nSweeps);
    params_mp->set("relaxation: zero starting solution", false);
    params_mp->set("relaxation: damping factor", damping);
    params_mp->set("relaxation: use l1", l1);
}
