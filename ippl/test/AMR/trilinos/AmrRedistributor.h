#ifndef AMR_REDISTRIBUTOR_H
#define AMR_REDISTRIBUTOR_H

#include "AmrMultiGridDefs.h"

#include <Zoltan2_PartitioningProblem.hpp>
#include <Zoltan2_XpetraCrsMatrixAdapter.hpp>
#include <Zoltan2_XpetraMultiVectorAdapter.hpp>

class AmrRedistributor {
    
public:
    typedef Zoltan2::XpetraCrsMatrixAdapter<amr::matrix_t> madapter_t;
    typedef Zoltan2::XpetraMultiVectorAdapter<amr::vector_t> vadapter_t;
    typedef Zoltan2::PartitioningProblem<madapter_t> partitioner_t;
    typedef Tpetra::Import<
        amr::local_ordinal_t,
        amr::global_ordinal_t,
        amr::node_t
    > import_t;
    typedef Tpetra::Export<
        amr::local_ordinal_t,
        amr::global_ordinal_t,
        amr::node_t
    > export_t;
    
public:
    
    AmrRedistributor();
    
    void solve(const Teuchos::RCP<amr::matrix_t>& matrix);
    
    void doImport(Teuchos::RCP<amr::matrix_t>& matrix);
    
    void doImport(Teuchos::RCP<amr::vector_t>& vector);
    
    void doExport(Teuchos::RCP<amr::matrix_t>& matrix);
    
    void doExport(Teuchos::RCP<amr::vector_t>& vector);
    
    void apply(Teuchos::RCP<amr::matrix_t>& matrix);
    
    void apply(Teuchos::RCP<amr::vector_t>& vector);
    
private:
    
    void initParameters_m();
    
private:
    Teuchos::RCP<import_t> importer_mp;
    Teuchos::RCP<export_t> exporter_mp;
    
    Teuchos::RCP<madapter_t> adapter_mp;        ///< input adapter for Tpetra matrix
    Teuchos::RCP<partitioner_t> problem_mp;     ///<
    Teuchos::ParameterList param_m;
};

#endif
