#include "MueLuPreconditioner.h"

#include <MueLu_CreateTpetraPreconditioner.hpp>



MueLuPreconditioner::MueLuPreconditioner()
    : prec_mp(Teuchos::null)
{
    this->init_m();
}


void MueLuPreconditioner::create(Teuchos::RCP<amr::matrix_t>& A) {
    prec_mp = MueLu::CreateTpetraPreconditioner(A, params_m);
}


Teuchos::RCP<amr::operator_t> MueLuPreconditioner::get() {
    return prec_mp;
}


void MueLuPreconditioner::init_m() {
    params_m.set("problem: type", "Poisson-3D");
    params_m.set("problem: symmertric", true);
    params_m.set("verbosity", "extreme");
    params_m.set("number of equations", 1);
    params_m.set("max levels", 3);
    params_m.set("cycle type", "V");
    params_m.set("coarse: max size", 10);
    params_m.set("multigrid algorithm", "sa");
    params_m.set("repartition", "enable");
    
    params_m.set("smoother: type", "CHEBYSHEV");
    params_m.set("smoother: pre and post", "both");
    
    params_m.set("coarse: type", "KLU2");
}
