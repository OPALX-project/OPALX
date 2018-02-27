#ifndef MUELU_PRECONDITIONER_H
#define MUELU_PRECONDITIONER_H

#include "AmrPreconditioner.h"
#include "Amr/AmrDefs.h"

#include <MueLu.hpp>
#include <MueLu_TpetraOperator.hpp>

template <class Level>
class MueLuPreconditioner : public AmrPreconditioner<amr::matrix_t, Level>
{
public:
    typedef amr::Preconditioner Preconditioner;
    
    typedef amr::scalar_t scalar_t;
    typedef amr::local_ordinal_t lo_t;
    typedef amr::global_ordinal_t go_t;
    typedef amr::AmrBox_t AmrBox_t;
    
    typedef MueLu::TpetraOperator<
        scalar_t,
        lo_t,
        go_t,
        amr::node_t
    > precond_t;
    
    typedef amr::AmrIntVect_t AmrIntVect_t;
    
    typedef std::map<std::string, Preconditioner> map_t;
    
public:
    
    MueLuPreconditioner(const bool& rebalance);
    
    void create(const Teuchos::RCP<amr::matrix_t>& A, Level* level_p =  nullptr);
    
    Teuchos::RCP<amr::operator_t> get();
    
    static void fillMap(map_t& map);
    
private:
    void init_m();

private:
    Teuchos::ParameterList params_m;
    
    Teuchos::RCP<precond_t> prec_mp;

    Teuchos::RCP<amr::multivector_t> coords_mp;
    
    const bool rebalance_m;
};

#include "MueLuPreconditioner.hpp"

#endif
