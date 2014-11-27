#ifndef ERROR_H
#define ERROR_H

#include <cstdlib>
#include <iostream>
#include <string>

/// This class writes error messages to the terminal and stops execution
/*!
 * It serves as a base class for further extensions.
 */
class Error
{
public:
  /// Writes an error message to the terminal
  /*!
   * @param sender specifying the class where the error happened
   * @param error describing the error that occurred
   */
  void static message(std::string sender, std::string error) {
    std::cerr << "ERROR in " << sender << ": " << error << std::endl;
    std::exit(EXIT_FAILURE);
  }
  
  /// Range-access error message
  static std::string range;
  /// Undefined quantities
  static std::string notdefined;
  /// Impossible operation
  static std::string invalid;
  /// Wrong dimension
  static std::string dim;
  /// Wrong size
  static std::string size;
};

std::string Error::range = "VALUE OUT OF RANGE";
std::string Error::notdefined = "VALUE NOT DEFINED";
std::string Error::invalid = "INVALID OPERATION";
std::string Error::dim = "INVALID DIMENSION";
std::string Error::size = "INVALID SIZE";

#endif