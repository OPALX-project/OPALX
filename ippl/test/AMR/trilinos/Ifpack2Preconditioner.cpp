#include "Ifpack2Preconditioner.h"

#include "Utilities/OpalException.h"

Ifpack2Preconditioner::Ifpack2Preconditioner(Preconditioner prec)
    : prec_mp(Teuchos::null)
{
    this->init_m(prec);
}


void Ifpack2Preconditioner::create(const Teuchos::RCP<amr::matrix_t>& A) {
    Ifpack2::Factory factory;
    
    prec_mp = factory.create(prectype_m, A.getConst());
    prec_mp->setParameters(*params_mp);
    prec_mp->initialize();
    prec_mp->compute();
}


Teuchos::RCP<amr::operator_t> Ifpack2Preconditioner::get() {
    return prec_mp;
}


void Ifpack2Preconditioner::fillMap(map_t& map) {
    map["ILUT"]         = Preconditioner::ILUT;
    map["CHEBYSHEV"]    = Preconditioner::CHEBYSHEV;
    map["RILUK"]        = Preconditioner::RILUK;
    map["JACOBI"]       = Preconditioner::JACOBI;
    map["BLOCK_JACOBI"] = Preconditioner::BLOCK_JACOBI;
    map["GS"]           = Preconditioner::GS;
    map["BLOCK_GS"]     = Preconditioner::BLOCK_GS;
}


void Ifpack2Preconditioner::init_m(Preconditioner prec)
{
    params_mp = Teuchos::parameterList();
    
    switch ( prec ) {
        case Preconditioner::ILUT:
            // inclomplete LU
            prectype_m = "ILUT";
            
            params_mp->set("fact: ilut level-of-fill", 5.0);
            params_mp->set("fact: drop tolerance", 0.0);
            
            break;
        case Preconditioner::CHEBYSHEV:
            prectype_m = "CHEBYSHEV";
            params_mp->set("chebyshev: degree", 1);
            break;
        case Preconditioner::RILUK:
            prectype_m = "RILUK";
            params_mp->set("fact: iluk level-of-fill", 0);
            params_mp->set("fact: relax value", 0.0);
            params_mp->set("fact: absolute threshold", 0.0);
            params_mp->set("fact: relative threshold", 1.0);
            break;
        case Preconditioner::JACOBI:
            prectype_m = "RELAXATION";
            params_mp->set("relaxation: type", "Jacobi");
            params_mp->set("relaxation: sweeps", 12);
            params_mp->set("relaxation: zero starting solution", false);
            params_mp->set("relaxation: damping factor", 6.0 / 7.0);
            params_mp->set("relaxation: use l1", true);
            params_mp->set("relaxation: l1 eta", 1.5);
            params_mp->set("relaxation: backward mode", false);
            params_mp->set("relaxation: fix tiny diagonal entries", true);
            params_mp->set("relaxation: min diagonal value", 1.0e-5);
            params_mp->set("relaxation: check diagonal entries", true);
            break;
        case Preconditioner::BLOCK_JACOBI:
            prectype_m = "BLOCK_RELAXATION";
            params_mp->set("relaxation: type", "Jacobi");
            params_mp->set("relaxation: sweeps", 12);
            params_mp->set("relaxation: zero starting solution", false);
            params_mp->set("relaxation: damping factor", 6.0 / 7.0);
            params_mp->set("relaxation: backward mode", false);
            
            params_mp->set("partitioner: type", "linear");
            params_mp->set("partitioner: overlap", 0);
            params_mp->set("partitioner: local parts", 1);
            
            break;
        case Preconditioner::GS:
            prectype_m = "RELAXATION";
            params_mp->set("relaxation: type", "Gauss-Seidel");
            params_mp->set("relaxation: sweeps", 12);
            params_mp->set("relaxation: zero starting solution", false);
            params_mp->set("relaxation: damping factor", 1.0);
            params_mp->set("relaxation: use l1", true);
            params_mp->set("relaxation: l1 eta", 1.5);
            params_mp->set("relaxation: backward mode", false);
            params_mp->set("relaxation: fix tiny diagonal entries", true);
            params_mp->set("relaxation: min diagonal value", 1.0e-5);
            params_mp->set("relaxation: check diagonal entries", true);
            break;
        case Preconditioner::BLOCK_GS:
            prectype_m = "BLOCK_RELAXATION";
            params_mp->set("relaxation: type", "Gauss-Seidel");
            params_mp->set("relaxation: sweeps", 12);
            params_mp->set("relaxation: zero starting solution", false);
            params_mp->set("relaxation: damping factor", 6.0 / 7.0);
            params_mp->set("relaxation: backward mode", false);
            
            params_mp->set("partitioner: type", "linear");
            params_mp->set("partitioner: overlap", 0);
            params_mp->set("partitioner: local parts", 1);
            break;
        case Preconditioner::NONE:
            prectype_m = "";
            break;
        default:
            throw OpalException("Ifpack2Preconditioner::init_m()",
                                "This type of Ifpack2 preconditioner not supported.");
    }
}
