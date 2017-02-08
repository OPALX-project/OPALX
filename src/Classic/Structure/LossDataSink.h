#ifndef OPAL_LOSSOUTPUT_H_
#define OPAL_LOSSOUTPUT_H_

//////////////////////////////////////////////////////////////
#include "Utility/IpplInfo.h"
#include "Algorithms/Vektor.h"
#include "AbsBeamline/ElementBase.h"

#include <string>
#include <fstream>
#include <vector>

#include <hdf5.h>
#include "H5hut.h"

/*
  - In the destructor we do ALL the file handling
  - h5hut_mode_m defines h5hut or ASCII
 */
class LossDataSink {
 public:
    LossDataSink();

    LossDataSink(std::string elem, bool hdf5Save, ElementBase::ElementType type = ElementBase::ANY);

    LossDataSink(const LossDataSink &rsh);
    ~LossDataSink();

    bool inH5Mode() { return h5hut_mode_m;}

    void save(unsigned int numSets = 1);

    void addReferenceParticle(const Vector_t &x,
                              const Vector_t &p,
                              double time,
                              double spos,
                              long long globalTrackStep);

    void addParticle(const Vector_t &x, const Vector_t &p, const size_t id);

    void addParticle(const Vector_t &x, const Vector_t &p, const size_t  id, const double time, const size_t turn);

    size_t size() const;

    static void writeStatistics();

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

    void openH5(h5_int32_t mode = H5_O_WRONLY);

    void saveH5(unsigned int setIdx);

    void saveASCII();

    void reportOnError(h5_int64_t rc, const char* file, int line);

    void splitSets(unsigned int numSets);
    void saveStatistics(unsigned int numSets);

    // filename without extension
    std::string fn_m;

    // used to write out the data
    std::ofstream os_m;

    std::string element_m;

    bool h5hut_mode_m;

    /// H5 file for particle data.
#if defined (USE_H5HUT2)
    h5_file_t H5file_m;
#else
    h5_file_t *H5file_m;
#endif

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

    std::vector<Vector_t> RefPartR_m;
    std::vector<Vector_t> RefPartP_m;
    std::vector<h5_int64_t> globalTrackStep_m;
    std::vector<double> refTime_m;
    std::vector<double> spos_m;

    std::vector<unsigned long> startSet_m;

    ElementBase::ElementType type_m;
    static std::map<double, std::string> statFileEntries_s;
};

inline
size_t LossDataSink::size() const {
    return x_m.size();
}

inline
void LossDataSink::reportOnError(h5_int64_t rc, const char* file, int line) {
    if (rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << file << " @ line " << line << endl);
}

#endif