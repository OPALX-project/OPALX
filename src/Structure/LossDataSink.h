#ifndef OPAL_LOSSOUTPUT_H_
#define OPAL_LOSSOUTPUT_H_

//////////////////////////////////////////////////////////////
//#include "Ippl.h"
#include "Utility/IpplInfo.h"
#include <string>
//#include <iosfwd>
#include <fstream>
#include <vector>
//#include <iostream>
#include "AbstractObjects/OpalData.h"

#include <hdf5.h>
#include "H5hut.h"

#include "config.h"

#define MAXLOSS 1000000

//////////////////////////////////////////////////////////////
// specialization of LossDataSink to PPartBunch objects
class LossDataSink {
public:

    LossDataSink(size_t np, bool hdf5Save = false):
        element_m("NULL"),
        doHdf5Save_m(hdf5Save) {
        if(!doHdf5Save_m) {
            fn_m = OpalData::getInstance()->getInputFn();
            int pos = fn_m.find(std::string("."), 0);
            fn_m.erase(pos, fn_m.size() - pos);
            fn_m += std::string(".loss");
            open();
            writeHeader(np);
        }
        hdf5FileIsOpen_m = false;
        H5call_m = 0;
    }

    LossDataSink():
        element_m("NULL") {
        if(!doHdf5Save_m) {
            fn_m = OpalData::getInstance()->getInputFn();
            int pos = fn_m.find(std::string("."), 0);
            fn_m.erase(pos, fn_m.size() - pos);
            fn_m += std::string(".loss");
            append();
        }
        hdf5FileIsOpen_m = false;
        H5call_m = 0;
    }

    ~LossDataSink() {
        save(element_m);
        if(!doHdf5Save_m) {
            close();
        } else {
            if(hdf5FileIsOpen_m) {
                H5CloseFile(H5file_m);
                Ippl::Comm->barrier();
            }
        }

    }

    bool inH5Mode() { return doHdf5Save_m;}

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

    void writeHeader(size_t np) {
        if(Ippl::myNode() == 0) {
            os_m << "# Element x (mm),  y (mm),  z (mm),  px ( ),  py ( ),  pz ( ), id,  turn,  time (ns) " << std::endl;
        }
    }

    void writeH5FileAttributes();

public:

    void addParticle(const Vector_t x, const Vector_t p, const size_t id);
    void addParticle_time(const Vector_t x, const Vector_t p, const size_t  id, const double time, const size_t turn); 
    void openH5(std::string fn);
    void save(std::string element) ;
    void save_time(std::string element) ;

private:
    // filename without extension
    std::string fn_m;

    // used to write out the data
    std::ofstream os_m;

    std::string element_m;

    bool doHdf5Save_m;

    bool hdf5FileIsOpen_m;

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



