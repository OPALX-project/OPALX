#ifndef OPAL_LOSSOUTPUT_H_
#define OPAL_LOSSOUTPUT_H_

//////////////////////////////////////////////////////////////
#include "Ippl.h"
#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include <set>
#include "AbstractObjects/OpalData.h"

#include <hdf5.h>
#include "H5Part.h"

#include "config.h"

#define MAXLOSS 1000000

using namespace std;

//////////////////////////////////////////////////////////////
// specialization of LossDataSink to PPartBunch objects
class LossDataSink {
public:

    LossDataSink(size_t np, bool hdf5Save = false):
        element_m("NULL"),
        doHdf5Save_m(hdf5Save) {
        if(!doHdf5Save_m) {
            fn_m = OPAL.getInputFn();
            int pos = fn_m.find(string("."), 0);
            fn_m.erase(pos, fn_m.size() - pos);
            fn_m += string(".loss");
            open();
            writeHeader(np);
        }
        hdf5FileIsOpen_m = false;
        H5call_m = 0;
    }

    LossDataSink():
        element_m("NULL") {
        if(!doHdf5Save_m) {
            fn_m = OPAL.getInputFn();
            int pos = fn_m.find(string("."), 0);
            fn_m.erase(pos, fn_m.size() - pos);
            fn_m += string(".loss");
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
            H5PartCloseFile(H5file_m);
            Ippl::Comm->barrier();
        }

    }

    bool inH5Mode() { return doHdf5Save_m;}

private:

    void open() {
        if(Ippl::myNode() == 0) {
            os_m.open(fn_m.c_str(), ios::out);
        }
    }

    void append() {
        if(Ippl::myNode() == 0) {
            os_m.open(fn_m.c_str(), ios::app);
        }
    }

    void close() {
        if(Ippl::myNode() == 0)
            os_m.close();
    }

    void writeHeader(size_t np) {
        if(Ippl::myNode() == 0) {
            os_m << "# Element x y z px py pz id | np = " << np << endl;
        }
    }

    void writeH5FileAttributes();

public:

    void addParticle(Vektor<double, 3> x, Vektor<double, 3> p, size_t id);
    void openH5(string fn);
    void save(string element) ;

private:
    // filename without extension
    string fn_m;

    // used to write out the data
    ofstream os_m;

    //  string element_m;
    string element_m;

    bool doHdf5Save_m;

    bool hdf5FileIsOpen_m;

    /// %Pointer to H5 file for particle data.
    H5PartFile *H5file_m;

    /// Current record, or time step, of H5 file.
    h5part_int64_t H5call_m;

    vector<long>   id_m;
    vector<double>  x_m;
    vector<double>  y_m;
    vector<double>  z_m;
    vector<double> px_m;
    vector<double> py_m;
    vector<double> pz_m;
};
#endif



