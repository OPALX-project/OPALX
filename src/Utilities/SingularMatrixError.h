// Copyright: see NOTICE and LICENSE.
#ifndef OPAL_SINGULAR_MATRIX_ERROR_H
#define OPAL_SINGULAR_MATRIX_ERROR_H

#include "Utilities/ArithmeticError.h"

/// A matrix is singular or a factorization pivot falls below the required floor.
class SingularMatrixError : public ArithmeticError {
public:
    explicit SingularMatrixError(const std::string& operation)
        : ArithmeticError(operation, "Singular matrix.") {}
};

#endif
