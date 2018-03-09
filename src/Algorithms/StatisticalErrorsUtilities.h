#ifndef STATISTICALERRORSUTILITIES_H
#define STATISTICALERRORSUTILITIES_H

//
//  Copyright & License: See Copyright.readme in src directory
//

/*!
  Class documentation
*/

#include <boost/variant/static_visitor.hpp>

#include <iostream>
#include <fstream>

class BinaryWriter: public boost::static_visitor<> {
public:
    BinaryWriter(std::ofstream &out);

    void operator()(const float &val) const;
    void operator()(const double &val) const;
    void operator()(const size_t &val) const;
    void operator()(const int &val) const;
    void operator()(const short &val) const;
    void operator()(const long &val) const;
    void operator()(const char &val) const;
    void operator()(const std::string &val) const;

private:
    std::ofstream &out_m;
};

class DataTypeWriter: public boost::static_visitor<> {
public:
    DataTypeWriter(std::ofstream &out);

    void operator()(const float &val) const;
    void operator()(const double &val) const;
    void operator()(const short &val) const;
    void operator()(const long &val) const;
    void operator()(const char &val) const;
    void operator()(const std::string &val) const;

private:
    std::ofstream &out_m;
};

#endif
