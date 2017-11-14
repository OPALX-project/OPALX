#include "AmrSmoother.h"
#include <string>
#include <utility>

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
                         const Teuchos::RCP<vector_t>& b,
                         const Teuchos::RCP<matrix_t>& S)
{
//     Teuchos::RCP<vector_t> residual = Teuchos::rcp( new vector_t(A->getDomainMap(), false) );
//     A->apply(*x, *residual);
//     residual->update(1.0, *b, -1.0);
    
    
//     Teuchos::RCP<vector_t> correction = Teuchos::rcp( new vector_t(A->getDomainMap(), false) );
    
//     prec_mp->apply(*residual, *correction, Teuchos::NO_TRANS,
//                    Teuchos::ScalarTraits<scalar_t>::one(),
//                    Teuchos::ScalarTraits<scalar_t>::zero());
    
//     Teuchos::RCP<vector_t> accel = Teuchos::rcp( new vector_t(A->getDomainMap(), false) );
//     S->apply(*b, *accel);
    
//     x->update(1.0, *accel, 1.0);
    
//     Teuchos::RCP<vector_t> xcopy = Teuchos::rcp( new vector_t(A->getDomainMap(), false) );
//     Tpetra::deep_copy(*xcopy, *x);
    
    prec_mp->apply(*b, *x, Teuchos::NO_TRANS,
                   Teuchos::ScalarTraits<scalar_t>::one(),
                   Teuchos::ScalarTraits<scalar_t>::zero());
//     x->update(1.0, *xcopy, 1.0);
}


void AmrSmoother::initParameter_m(const Smoother& smoother,
                                  local_ordinal_t nSweeps)
{
    if ( params_mp == Teuchos::null )
        params_mp = Teuchos::rcp( new Teuchos::ParameterList );
    
    
    std::string type = "";
    double damping = 1.0;
    std::pair<bool, double> l1 = std::make_pair(true, 1.5);
    
    bool backward = false;
    std::pair<bool, double> fix = std::make_pair(true, 1.0e-5);
    bool check = true;
    
    switch ( smoother ) {
        
        case GAUSS_SEIDEL:
        {
            type = "Gauss-Seidel";
            backward = false;
            break;
        }
        case SGS:
        {
            type = "Symmetric Gauss-Seidel";
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
    params_mp->set("relaxation: use l1", l1.first);
    params_mp->set("relaxation: l1 eta", l1.second);
    params_mp->set("relaxation: backward mode", backward);
    params_mp->set("relaxation: fix tiny diagonal entries", fix.first);
    params_mp->set("relaxation: min diagonal value", fix.second);
    params_mp->set("relaxation: check diagonal entries", check);
}
