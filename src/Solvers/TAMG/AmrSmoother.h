#ifndef AMR_SMOOTHER_H
#define AMR_SMOOTHER_H

#include "AmrMultiGridCore.h"

#include "Ifpack2_Factory.hpp"

class AmrSmoother {
    
public:
    typedef amr::global_ordinal_t global_ordinal_t;
    typedef amr::local_ordinal_t local_ordinal_t;
    typedef amr::scalar_t scalar_t;
    typedef amr::node_t node_t;
    typedef amr::matrix_t matrix_t;
    typedef amr::vector_t vector_t;
    typedef Ifpack2::Preconditioner<scalar_t,
                                    local_ordinal_t,
                                    global_ordinal_t,
                                    node_t
                                    > preconditioner_t;
    
    enum Smoother {
        GAUSS_SEIDEL = 0,
        SGS,    // symmetric Gauss-Seidel
        JACOBI //,
//         SOR
    };
    
public:
    
    AmrSmoother(const Teuchos::RCP<const matrix_t>& A,
                const Smoother& smoother,
                local_ordinal_t nSweeps);
    
    ~AmrSmoother();
    
    void smooth(const Teuchos::RCP<vector_t>& x,
                const Teuchos::RCP<matrix_t>& A,
                const Teuchos::RCP<vector_t>& b);
    
    
private:
    
    void initParameter_m(const Smoother& smoother,
                         local_ordinal_t nSweeps);
    
    
private:
    Teuchos::RCP<preconditioner_t> prec_mp;
    
    Teuchos::RCP<Teuchos::ParameterList> params_mp;
};

#endif
