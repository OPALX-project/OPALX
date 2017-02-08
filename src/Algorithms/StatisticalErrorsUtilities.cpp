//
//  Copyright & License: See Copyright.readme in src directory
//

#include "Algorithms/StatisticalErrorsUtilities.h"
#include "Structure/SDDSParser/ast.hpp"

BinaryWriter::BinaryWriter(std::ofstream &out):
    out_m(out)
{ }

void BinaryWriter::operator()(const float &val) const {
    out_m.write(reinterpret_cast<const char*>(&val), sizeof(float));
}

void BinaryWriter::operator()(const double &val) const {
    out_m.write(reinterpret_cast<const char*>(&val), sizeof(double));
}

void BinaryWriter::operator()(const size_t &val) const {
    out_m.write(reinterpret_cast<const char*>(&val), sizeof(size_t));
}

void BinaryWriter::operator()(const int &val) const {
    out_m.write(reinterpret_cast<const char*>(&val), sizeof(int));
}

void BinaryWriter::operator()(const short &val) const {
    out_m.write(reinterpret_cast<const char*>(&val), sizeof(short));
}

void BinaryWriter::operator()(const long &val) const {
    out_m.write(reinterpret_cast<const char*>(&val), sizeof(long));
}

void BinaryWriter::operator()(const char &val) const
{ }

void BinaryWriter::operator()(const std::string & val) const
{ }


DataTypeWriter::DataTypeWriter(std::ofstream &out):
    out_m(out)
{ }

void DataTypeWriter::operator()(const float &) const {
    int dataType = SDDS::ast::FLOAT;
    out_m.write(reinterpret_cast<const char*>(&dataType), sizeof(int));
}

void DataTypeWriter::operator()(const double &) const {
    int dataType = SDDS::ast::DOUBLE;
    out_m.write(reinterpret_cast<const char*>(&dataType), sizeof(int));
}

void DataTypeWriter::operator()(const short &) const {
    int dataType = SDDS::ast::SHORT;
    out_m.write(reinterpret_cast<const char*>(&dataType), sizeof(int));
}

void DataTypeWriter::operator()(const long &) const {
    int dataType = SDDS::ast::LONG;
    out_m.write(reinterpret_cast<const char*>(&dataType), sizeof(int));
}

void DataTypeWriter::operator()(const char &) const {
    int dataType = SDDS::ast::CHARACTER;
    out_m.write(reinterpret_cast<const char*>(&dataType), sizeof(int));
}

void DataTypeWriter::operator()(const std::string &) const {
    int dataType = SDDS::ast::STRING;
    out_m.write(reinterpret_cast<const char*>(&dataType), sizeof(int));
}
