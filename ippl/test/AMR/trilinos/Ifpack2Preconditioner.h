#ifndef IFPACK_2_PRECONDITIONER_H
#define IFPACK_2_PRECONDITIONER_H

#include "AmrPreconditioner.h"

#include <Ifpack2_Factory.hpp>

class Ifpack2Preconditioner : public AmrPreconditioner<amr::matrix_t>
{
public:
    typedef amr::Preconditioner Preconditioner;
    typedef Ifpack2::Preconditioner<
        amr::scalar_t,
        amr::local_ordinal_t,
        amr::global_ordinal_t,
        amr::node_t
    > precond_t;
    
public:
    
    Ifpack2Preconditioner(Preconditioner prec);
    
    void create(const Teuchos::RCP<amr::matrix_t>& A);
    
    Teuchos::RCP<amr::operator_t> get();
    
private:
    /*!
     * Initialize preconditioner parameter list
     */
    void init_m(Preconditioner prec);
    
private:
    /// preconditioner type
    std::string prectype_m;
    
    /// parameter list of preconditioner
    Teuchos::RCP<Teuchos::ParameterList> params_mp;
    
    Teuchos::RCP<precond_t> prec_mp;
};

#endif
