#ifndef OPALX_FIELDMAP_ICC
#define OPALX_FIELDMAP_ICC


#include <string>
#include <iostream>
#include <fstream>
#include <typeindex>
#include "Fields/Fieldmap.h"

/**
 * @brief Parses a single value of type T from a line in the stream.
 * 
 * Reads a line using getLine(), then attempts to extract one value.
 * Validates that the extraction matches the type and that no unexpected
 * data remains on the line.
 * 
 * @param in Input file stream.
 * @param value Output variable reference.
 * @param file_length_known If false, EOF is acceptable.
 * @return true if successful, false if EOF (and permitted).
 */
template<class T>
bool Fieldmap::interpretLine(std::ifstream & in, T & value, const bool & file_length_known)
{
    bool read_all = true;
    std::string buffer;
    std::streampos start = in.tellg();

    getLine(in, buffer);
    if (in.eof() && !file_length_known) {
        return false;
    }
    std::istringstream interpreter(buffer, std::istringstream::in);

    interpreter >> value;

    if (interpreter.rdstate() ^ std::ios_base::eofbit) {
        in.seekg(start);
        if (in.eof()) {
            missingValuesWarning();
            return false;
        }
        std::string expecting(std::string(typeid(T).name()));
        interpretWarning((interpreter.rdstate() ^ std::ios_base::eofbit), read_all, expecting, buffer);
    }
    return (!(interpreter.rdstate() ^ std::ios_base::eofbit) && read_all);   // eof should not be an error but if not eof
}

/**
 * @brief Parses two values (types S, T) from a line.
 * @see interpretLine<T>
 */
template<class S, class T>
bool Fieldmap::interpretLine(std::ifstream & in, S & value1, T & value2, const bool & file_length_known)
{
    bool read_all = true;
    std::string buffer;
    std::streampos start = in.tellg();

    getLine(in, buffer);
    if (in.eof() && !file_length_known) {
        return false;
    }
    std::istringstream interpreter(buffer, std::istringstream::in);

    interpreter >> value1;
    if (interpreter.rdstate()) read_all = false;

    interpreter >> value2;

    if (interpreter.rdstate() ^ std::ios_base::eofbit) {
        in.seekg(start);
        if (in.eof()) {
            missingValuesWarning();
            return false;
        }
        std::string expecting(std::string(typeid(S).name()));
        expecting += std::string(" ") + std::string(typeid(T).name());
        interpretWarning((interpreter.rdstate() ^ std::ios_base::eofbit), read_all, expecting, buffer);
    }
    return (!(interpreter.rdstate() ^ std::ios_base::eofbit) && read_all);   // eof should not be an error but if not eof
}

/**
 * @brief Parses three values (types S, T, U) from a line.
 * @see interpretLine<T>
 */
template<class S, class T, class U>
bool Fieldmap::interpretLine(std::ifstream & in, S & value1, T & value2, U & value3, const bool & file_length_known)
{
    bool read_all = true;
    std::string buffer;
    std::streampos start = in.tellg();

    getLine(in, buffer);
    if (in.eof() && !file_length_known) {
        return false;
    }
    std::istringstream interpreter(buffer, std::istringstream::in);

    interpreter >> value1;

    interpreter >> value2;
    if (interpreter.rdstate()) read_all = false;

    interpreter >> value3;

    if (interpreter.rdstate() ^ std::ios_base::eofbit) {
        in.seekg(start);
        if (in.eof()) {
            missingValuesWarning();
            return false;
        }
        std::string expecting(std::string(typeid(S).name()) + std::string(" ") +
                              std::string(typeid(T).name()) + std::string(" ") +
                              std::string(typeid(U).name()));
        interpretWarning((interpreter.rdstate() ^ std::ios_base::eofbit), read_all, expecting, buffer);
    }
    return (!(interpreter.rdstate() ^ std::ios_base::eofbit) && read_all);   // eof should not be an error but if not eof
}

/**
 * @brief Parses four values (types S, T, U, V) from a line.
 * @see interpretLine<T>
 */
template<class S, class T, class U, class V>
bool Fieldmap::interpretLine(std::ifstream & in, S & value1, T & value2, U & value3, V & value4, const bool & file_length_known)
{
    bool read_all = true;
    std::string buffer;
    std::streampos start = in.tellg();

    getLine(in, buffer);
    if (in.eof() && !file_length_known) {
        return false;
    }
    std::istringstream interpreter(buffer, std::istringstream::in);

    interpreter >> value1;

    interpreter >> value2;

    interpreter >> value3;
    if (interpreter.rdstate()) read_all = false;

    interpreter >> value4;

    if (interpreter.rdstate() ^ std::ios_base::eofbit) {
        in.seekg(start);
        if (in.eof()) {
            missingValuesWarning();
            return false;
        }
        std::string expecting(std::string(typeid(S).name()) + std::string(" ") +
                              std::string(typeid(T).name()) + std::string(" ") +
                              std::string(typeid(U).name()) + std::string(" ") +
                              std::string(typeid(V).name()));
        interpretWarning((interpreter.rdstate() ^ std::ios_base::eofbit), read_all, expecting, buffer);
    }
    return (!(interpreter.rdstate() ^ std::ios_base::eofbit) && read_all);   // eof should not be an error but if not eof
}

/**
 * @brief Parses six values of the same type S from a line.
 * @see interpretLine<T>
 */
template<class S>
bool Fieldmap::interpretLine(std::ifstream & in, S & value1, S & value2, S & value3, S & value4, S & value5, S & value6, const bool & file_length_known)
{
    bool read_all = true;
    std::string buffer;
    std::streampos start = in.tellg();

    getLine(in, buffer);
    if (in.eof() && !file_length_known) {
        return false;
    }
    std::istringstream interpreter(buffer, std::istringstream::in);

    interpreter >> value1;

    interpreter >> value2;

    interpreter >> value3;

    interpreter >> value4;

    interpreter >> value5;
    if (interpreter.rdstate()) read_all = false;

    interpreter >> value6;

    if (interpreter.rdstate() ^ std::ios_base::eofbit) {
        in.seekg(start);
        if (in.eof()) {
            missingValuesWarning();
            return false;
        }
        std::string expecting(std::string(typeid(S).name()) + std::string(" ") +
                              std::string(typeid(S).name()) + std::string(" ") +
                              std::string(typeid(S).name()) + std::string(" ") +
                              std::string(typeid(S).name()) + std::string(" ") +
                              std::string(typeid(S).name()) + std::string(" ") +
                              std::string(typeid(S).name()));
        interpretWarning((interpreter.rdstate() ^ std::ios_base::eofbit), read_all, expecting, buffer);
    }
    return (!(interpreter.rdstate() ^ std::ios_base::eofbit) && read_all);   // eof should not be an error but if not eof
}


#endif