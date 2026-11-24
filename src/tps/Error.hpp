// Copyright: see Copyright.readme (original OPAL/CLASSIC sources).
#ifndef TPS_ERROR_HPP
#define TPS_ERROR_HPP
#include <stdexcept>
#include <string>

namespace tps {
/** @brief Base for algebra failures; what() describes the failure, where() its operation.
 * @details Catch Error to handle any algebra failure, or a derived type to distinguish
 * a dimension/range failure from invalid arithmetic. Map metadata/preparation validation
 * uses standard exceptions; kernel evaluation returns EvaluationStatus instead.
 */
class Error : public std::runtime_error {
  public:
    Error(const std::string& operation, const std::string& message)
        : std::runtime_error(message), operation_(operation) {}
    [[nodiscard]] const std::string& where() const noexcept {
        return operation_;
    }

  private:
    std::string operation_;
};
/** @brief An index or order lies outside its supported range. */
class CLRangeError : public Error {
  public:
    using Error::Error;
};
/** @brief Array, matrix or polynomial dimensions are incompatible. */
class SizeError : public Error {
  public:
    using Error::Error;
};
/** @brief Serialized polynomial input has an invalid format. */
class FormatError : public Error {
  public:
    using Error::Error;
};
/** @brief An algebra operation violates its required preconditions. */
class LogicalError : public Error {
  public:
    using Error::Error;
};
/** @brief An algebra operation divides by zero. */
class DivideError : public Error {
  public:
    explicit DivideError(const std::string& operation) : Error(operation, "Divide error.") {}
};
/** @brief An elementary function is undefined at the expansion point. */
class DomainError : public Error {
  public:
    explicit DomainError(const std::string& operation) : Error(operation, "Domain error.") {}
};
/** @brief A linear solve encounters a singular matrix. */
class SingularMatrixError : public Error {
  public:
    explicit SingularMatrixError(const std::string& operation)
        : Error(operation, "Singular matrix.") {}
};
} // namespace tps
#endif
