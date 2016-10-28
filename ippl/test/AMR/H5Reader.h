#ifndef H5READER_H
#define H5READER_H

#include <string>

#include <H5hut.h>

#include "Distribution.h"

#define READDATA(type, file, name, value) H5PartReadData##type(file, name, value);

class H5Reader {
    
public:
    
    H5Reader(const std::string& filename);
    
    void open(int step);
    void close();
    void read(Distribution::container_t& x,
              Distribution::container_t& px,
              Distribution::container_t& y,
              Distribution::container_t& py,
              Distribution::container_t& z,
              Distribution::container_t& pz,
              size_t firstParticle,
              size_t lastParticle);
    
    
    h5_ssize_t getNumParticles();
    
private:
    std::string filename_m;
    h5_file_t* file_m;
    
};

#endif