#ifndef H5READER_H
#define H5READER_H

#include <string>

#include <H5hut.h>

#include "Distribution.h"

#define READDATA(type, file, name, value) H5PartReadData##type(file, name, value);


/*!
 * @file H5Reader.h
 * @details This class is implemented according to
 * src/Structure/H5PartWrapper and
 * src/Structure/H5PartWrapperForPC.\n It is used in
 * Distribution.h to read in a particle distribution
 * from a H5 file (supports only opal flavour: opal-cycl).
 * @author Matthias Frey
 * @date 25. - 26. October 2016
 * @brief Read in a particle distribution from a H5 file
 */

/// Defines functions for reading a distribution from a cyclotron H5 file
class H5Reader {
    
public:
    /*!
     * @param filename specifies path and name of the file
     */
    H5Reader(const std::string& filename);
    
    /*!
     * Open the file and set the step to read
     * @param step to read in
     */
    void open(int step);
    
    /*!
     * Close the file and sets the pointer to NULL
     */
    void close();
    
    /*!
     * Copy the particle distribution to the containers of
     * the Distribution class. (parallel read)
     * It performs a shift in longitduinal direction due to the
     * box extent
     * @param x - coordinate
     * @param px - coordinate
     * @param y - coordinate
     * @param py - coordinate
     * @param z - coordinate
     * @param pz - coordinate
     * @param q is the particle charge
     * @param firstParticle to read (core specific)
     * @param lastParticle to read (core specific)
     */
    void read(Distribution::container_t& x,
              Distribution::container_t& px,
              Distribution::container_t& y,
              Distribution::container_t& py,
              Distribution::container_t& z,
              Distribution::container_t& pz,
              Distribution::container_t& q,
              size_t firstParticle,
              size_t lastParticle);
    
    /*!
     * @returns the number of particles
     */
    h5_ssize_t getNumParticles();
    
private:
    std::string filename_m;     ///< Path and filename
    h5_file_t* file_m;          ///< Pointer to the opened file
    
};

#endif