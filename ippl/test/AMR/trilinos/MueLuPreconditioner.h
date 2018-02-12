#ifndef MUELU_PRECONDITIONER_H
#define MUELU_PRECONDITIONER_H

#include "AmrPreconditioner.h"

#include "MueLu.hpp"

class MueLuPreconditioner : public AmrPreconditioner<amr::matrix_t>
{
public:
    typedef MueLu::TpetraOperator precond_t;
    
public:
    
    MueLuPreconditioner();
    
    void create(const Teuchos::RCP<const amr::matrix_t>& A);
    
    Teuchos::RCP<amr::operator_t> get();
    
private:
    void init_m();
    
private:
    Teuchos::ParameterList params_m;
    
    precond_t prec_mp;
};

#endif
