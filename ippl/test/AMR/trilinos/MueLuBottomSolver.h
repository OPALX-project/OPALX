#ifndef MUELU_BOTTOM_SOLVER_H
#define MUELU_BOTTOM_SOLVER_H

#include "BottomSolver.h"

#include <MueLu.hpp>
#include <MueLu_Level.hpp>
#include <MueLu_Utilities.hpp>

#include <MueLu_ParameterListInterpreter.hpp>

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
//    typedef amr::AmrGeometry_t AmrGeometry_t;
    
    typedef MueLu::Hierarchy<scalar_t, lo_t, go_t, node_t> hierarchy_t;
    typedef MueLu::Level level_t;
    typedef Xpetra::Matrix<scalar_t, lo_t, go_t, node_t> xmatrix_t;
    typedef Xpetra::MultiVector<scalar_t, lo_t, go_t, node_t> xmv_t;
    typedef MueLu::Utilities<scalar_t, lo_t, go_t, node_t> util_t;
    
    typedef MueLu::ParameterListInterpreter<scalar_t, lo_t, go_t, node_t> pListInterpreter_t;
    typedef MueLu::HierarchyManager<scalar_t, lo_t, go_t, node_t> manager_t;
        
public:
    
    MueLuBottomSolver(/*const AmrGeometry_t& geom*/);
    
    void solve(const Teuchos::RCP<mv_t>& x,
               const Teuchos::RCP<mv_t>& b);
    
    void setOperator(const Teuchos::RCP<matrix_t>& A);
    
    std::size_t getNumIters();

private:
    void initMueLuList_m();
    
private:
    Teuchos::RCP<hierarchy_t> hierarchy_mp;     ///< manages the multigrid hierarchy
    
    Teuchos::RCP<level_t> finest_mp;            ///< finest level of hierarchy

    Teuchos::RCP<xmatrix_t> A_mp;               ///< MueLu requires Xpetra

    lo_t nSweeps_m;                             ///< the number of multigrid iterations
    
    Teuchos::ParameterList mueluList_m;

//    const AmrGeometry_t& geom_mr;               ///< reference to bottom geometry
    
};

#endif
