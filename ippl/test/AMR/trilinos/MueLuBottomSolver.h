#ifndef MUELU_BOTTOM_SOLVER_H
#define MUELU_BOTTOM_SOLVER_H

#include "BottomSolver.h"

#include <MueLu.hpp>
#include <MueLu_Level.hpp>

class MueLuBottomSolver : public BottomSolver<Teuchos::RCP<amr::matrix_t>,
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
    
    typedef MueLu::Hierarchy<scalar_t, lo_t, go_t, node_t> hierarchy_t;
    typedef MueLu::Level level_t;
    
public:
    
    MueLuBottomSolver();
    
    void solve(const Teuchos::RCP<mv_t>& x,
               const Teuchos::RCP<mv_t>& b);
    
    void setOperator(const Teuchos::RCP<matrix_t>& A);
    
    std::size_t getNumIters();
    
    
private::
    Teuchos::RCP<hierarchy_t> hierarchy_mp;     ///< manages the multigrid hierarchy
    
    Teuchos:RCP<level_t> finest_mp;             ///< finest level of hierarchy
    
    hierarchy_t::ConvData conv_m;               ///< stopping criteria of multigrid iteration
    
};

#endif
