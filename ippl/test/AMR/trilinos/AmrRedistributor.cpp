#include "AmrRedistributor.h"

#include "Utilities/OpalException.h"

AmrRedistributor::AmrRedistributor()
    : importer_mp(Teuchos::null),
      exporter_mp(Teuchos::null)
{
    this->initParameters_m();
}


void AmrRedistributor::solve(const Teuchos::RCP<amr::matrix_t>& matrix) {
    adapter_mp = Teuchos::rcp(new madapter_t(matrix));
    
    problem_mp = Teuchos::rcp(new partitioner_t(adapter_mp.get(), &param_m));
    
    try {
        problem_mp->solve();
    } catch (const std::exception& ex) {
        // throw an OPAL exception instead
        throw OpalException("AmrRedistributor::solve()", ex.what());
    }
    
    Teuchos::RCP<amr::matrix_t> redistributed;
    adapter_mp->applyPartitioningSolution(*matrix, redistributed,
                                          problem_mp->getSolution());
    
    importer_mp = Teuchos::rcp( new import_t(matrix->getMap(),
                                             redistributed->getMap()) );
    
    exporter_mp = Teuchos::rcp( new export_t(redistributed->getMap(),
                                             matrix->getMap()) );
}


void AmrRedistributor::doImport(Teuchos::RCP<amr::matrix_t>& matrix) {
    Teuchos::RCP<amr::matrix_t> redistributed = Teuchos::rcp( new amr::matrix_t(importer_mp->getTargetMap(), 0) );
    redistributed->doImport(*matrix, *importer_mp, Tpetra::REPLACE);
    redistributed->fillComplete();
    matrix = redistributed->clone(matrix->getNode());
}


void AmrRedistributor::doImport(Teuchos::RCP<amr::vector_t>& vector) {
    Teuchos::RCP<amr::vector_t> redistributed = Teuchos::rcp( new amr::vector_t(importer_mp->getTargetMap(), 0) );
    redistributed->doImport(*vector, *importer_mp, Tpetra::REPLACE);
    Tpetra::deep_copy(*vector, *redistributed);
}


void AmrRedistributor::doExport(Teuchos::RCP<amr::matrix_t>& matrix) {
    std::cout << "export matrix" << std::endl;
    Teuchos::RCP<amr::matrix_t> redistributed = Teuchos::rcp( new amr::matrix_t(exporter_mp->getTargetMap(), 0) );
    redistributed->doExport(*matrix, *exporter_mp, Tpetra::REPLACE);
    std::cout << "HI" << std::endl;
    redistributed->fillComplete();
    std::cout << "Bye" << std::endl;
    matrix = redistributed->clone(matrix->getNode());
    std::cout << "Done" << std::endl;
}


void AmrRedistributor::doExport(Teuchos::RCP<amr::vector_t>& vector) {
    Teuchos::RCP<amr::vector_t> redistributed = Teuchos::rcp( new amr::vector_t(exporter_mp->getTargetMap(), 0) );
    redistributed->doExport(*vector, *exporter_mp, Tpetra::REPLACE);
    Tpetra::deep_copy(*vector, *redistributed);
}


void AmrRedistributor::apply(Teuchos::RCP<amr::matrix_t>& matrix) {
    Teuchos::RCP<amr::matrix_t> redistributed;
    adapter_mp->applyPartitioningSolution(*matrix, redistributed,
                                          problem_mp->getSolution());
    matrix = redistributed->clone(matrix->getNode());
}


void AmrRedistributor::apply(Teuchos::RCP<amr::vector_t>& vector) {
    Teuchos::RCP<amr::vector_t> redistributed;
    vadapter_t adapter(vector);
    adapter.applyPartitioningSolution(*vector, redistributed,
                                      problem_mp->getSolution());
    Tpetra::deep_copy(*vector, *redistributed);
}


void AmrRedistributor::initParameters_m() {
    param_m.set("partitioning_approach", "partition");
    param_m.set("algorithm", "parmetis");
}