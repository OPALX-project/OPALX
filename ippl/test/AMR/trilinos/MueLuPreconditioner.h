#ifndef MUELU_PRECONDITIONER_H
#define MUELU_PRECONDITIONER_H

#include "AmrPreconditioner.h"

#include <MueLu.hpp>
#include <MueLu_TpetraOperator.hpp>

class MueLuPreconditioner : public AmrPreconditioner<amr::matrix_t>
{
public:
    typedef MueLu::TpetraOperator<
        amr::scalar_t,
        amr::local_ordinal_t,
        amr::global_ordinal_t,
        amr::node_t
    > precond_t;
    
public:
    
    MueLuPreconditioner();
    
    void create(const Teuchos::RCP<amr::matrix_t>& A);
    
    Teuchos::RCP<amr::operator_t> get();
    
private:
    void init_m();
    
private:
    Teuchos::ParameterList params_m;
    
    Teuchos::RCP<precond_t> prec_mp;
};

#endif
