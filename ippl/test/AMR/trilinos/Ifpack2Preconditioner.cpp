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


void Ifpack2Preconditioner::init_m(Preconditioner prec)
{
    params_mp = Teuchos::parameterList();
    
    switch ( prec ) {
        case Preconditioner::ILUT:
            // inclomplete LU
            prectype_m = "ILUT";
            
            params_mp->set("fact: ilut level-of-fill", 1.0);
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
        case Preconditioner::NONE:
            prectype_m = "";
            break;
        default:
            throw OpalException("Ifpack2Preconditioner::init_m()",
                                "This type of Ifpack2 preconditioner not supported.");
    }
}
