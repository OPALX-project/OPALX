#ifndef OPAL_LOSSOUTPUT_H_
#define OPAL_LOSSOUTPUT_H_

//////////////////////////////////////////////////////////////
#include "Utility/IpplInfo.h"
#include "Algorithms/Vektor.h"

#include <string>
#include <fstream>
#include <vector>

#include <hdf5.h>
#include "H5hut.h"

#include "config.h"

/*
  - In the destructor we do ALL the file handling
  - h5hut_mode_m defines h5hut or ASCII
 */
class LossDataSink {
 public:
    LossDataSink();

    LossDataSink(std::string elem, bool hdf5Save);

    LossDataSink(const LossDataSink &rsh);
    ~LossDataSink();

    bool inH5Mode() { return h5hut_mode_m;}

    void save();

    void addParticle(const Vector_t x, const Vector_t p, const size_t id);

    void addParticle(const Vector_t x, const Vector_t p, const size_t  id, const double time, const size_t turn);

private:

    void open() {
        if(Ippl::myNode() == 0) {
            os_m.open(fn_m.c_str(), std::ios::out);
        }
    }

    void append() {
        if(Ippl::myNode() == 0) {
            os_m.open(fn_m.c_str(), std::ios::app);
        }
    }

    void close() {
        if(Ippl::myNode() == 0)
            os_m.close();
    }

    void writeHeaderASCII() {
        if(Ippl::myNode() == 0) {
            if (time_m.size() != 0)
                os_m << "# Element " << element_m << " x (mm),  y (mm),  z (mm),  px ( ),  py ( ),  pz ( ), id,  turn,  time (ns) " << std::endl;
            else
                os_m << "# Element " << element_m << " x (mm),  y (mm),  z (mm),  px ( ),  py ( ),  pz ( ), id " << std::endl;
        }
    }

    bool hasNoParticlesToDump();

    bool hasTimeAttribute();

    void writeHeaderH5();

    void openH5();

    void saveH5();

    void saveASCII();

private:

    // filename without extension
    std::string fn_m;

    // used to write out the data
    std::ofstream os_m;

    std::string element_m;

    bool h5hut_mode_m;

    /// %Pointer to H5 file for particle data.
    h5_file_t *H5file_m;

    /// Current record, or time step, of H5 file.
    h5_int64_t H5call_m;

    std::vector<long>   id_m;
    std::vector<double>  x_m;
    std::vector<double>  y_m;
    std::vector<double>  z_m;
    std::vector<double> px_m;
    std::vector<double> py_m;
    std::vector<double> pz_m;

    std::vector<size_t> turn_m;
    std::vector<double> time_m;
};
#endif
