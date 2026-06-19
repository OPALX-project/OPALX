#ifndef OPALX_FIELD_STORAGE_HPP
#define OPALX_FIELD_STORAGE_HPP

#include "PartBunch/FieldContainer.hpp"

/**
 * @brief Runtime field storage used by `PartBunch` and Poisson backends.
 *
 * `FieldContainer` is kept as the implementation name for compatibility during the staged
 * refactor. New code should use `FieldStorage` at ownership boundaries so the final rename can be
 * mechanical.
 */
template <typename T, unsigned Dim = 3>
using FieldStorage = FieldContainer<T, Dim>;

#endif  // OPALX_FIELD_STORAGE_HPP
