//
//  Copyright & License: See Copyright.readme in src directory
//

#include "Structure/DataSink.h"

#include "config.h"
#include "Algorithms/bet/EnvelopeBunch.h"
#include "AbstractObjects/OpalData.h"
#include "Utilities/OpalOptions.h"
#include "Utilities/Options.h"
#include "Fields/Fieldmap.hh"
#include "Structure/BoundaryGeometry.h"
#include "Utilities/Timer.h"

#include "H5hut.h"

#include <queue>

extern Inform *gmsg;

using namespace std;

//using Physics::m_e;

// backward compatibility with 1.99.5
#if defined(H5_O_FLUSHSTEP)
#define H5_FLUSH_STEP H5_O_FLUSHSTEP
#endif

DataSink::DataSink() :
    lossWrCounter_m(0),
    doHDF5_m(true),
    H5file_m(NULL)
{
    /// Constructor steps:
    /// Get timers from IPPL.
    H5PartTimer_m = IpplTimings::getTimer("H5PartTimer");
    H5BlockTimer_m = IpplTimings::getTimer("H5BlockTimer");
    StatMarkerTimer_m = IpplTimings::getTimer("StatMarkerTimer");

    /// Set file write flags to true. These will be set to false after first
    /// write operation.
    firstWriteToStat_m = true;
    firstWriteH5part_m = true;
    firstWriteH5Surface_m = true;
    /// Define file names.
    string fn = OpalData::getInstance()->getInputBasename();
    surfaceLossFileName_m = fn + string(".SurfaceLoss.h5");
    statFileName_m = fn + string(".stat");
    lBalFileName_m = fn + string(".lbal");

    fn += string(".h5");

    doHDF5_m = Options::enableHDF5;

    if (doHDF5_m) {
        /// Open H5 file. Check that it opens correctly.
#ifdef PARALLEL_IO
        H5file_m = H5OpenFile(fn.c_str(), H5_FLUSH_STEP | H5_O_WRONLY, Ippl::getComm());
#else
        H5file_m = H5OpenFile(fn.c_str(), H5_FLUSH_STEP | H5_O_WRONLY, 0);
#endif
        if (H5file_m == (void*)H5_ERR) {
            throw OpalException("DataSink::DataSink()",
                                "failed to open h5 file '" + fn + "'");
        }
        /// Write file attributes.
        writeH5FileAttributes();
    }
    else {
        H5file_m = NULL;
    }
    /// Set current record/time step to 0.
    H5call_m = 0;
}

DataSink::DataSink(int restartStep) :
    lossWrCounter_m(0),
    H5file_m(NULL)
{
    doHDF5_m = Options::enableHDF5;
    if (!doHDF5_m) {
        throw OpalException("DataSink::DataSink(int)",
                            "Can not restart when HDF5 is disabled");
    }
    /// Constructor steps:
    h5_int64_t rc;
    /// Get timers from IPPL.
    H5PartTimer_m = IpplTimings::getTimer("H5PartTimer");
    H5BlockTimer_m = IpplTimings::getTimer("H5BlockTimer");
    StatMarkerTimer_m = IpplTimings::getTimer("StatMarkerTimer");

    /// Set file write flags. Since this is a restart, assume H5 file is old.
    firstWriteToStat_m = true;
    firstWriteH5part_m = false;
    /// Get file name root.
    string fn = OpalData::getInstance()->getInputBasename();

    statFileName_m = fn + string(".stat");
    lBalFileName_m = fn + string(".lbal");
    /// Test if .stat and .lbal files exist. If they do, new data will be appended.
    ofstream statFile(statFileName_m.c_str(), ios::in);
    ofstream lBalFile(lBalFileName_m.c_str(), ios::in);

    if(statFile.is_open()) {
        // File exists so we append data to end.
        firstWriteToStat_m = false;
        statFile.close();
        *gmsg << "Appending statistical data to existing data file: " << statFileName_m << endl;
    } else {
        statFile.clear();
        *gmsg << "Creating new file for statistical data: " << statFileName_m << endl;
    }

    if(lBalFile.is_open()) {
        // File exists so we append data to end.
        lBalFile.close();
        *gmsg << "Appending load balance data to existing data file: " << lBalFileName_m << endl;
    } else {
        lBalFile.clear();
        *gmsg << "Creating new file for load balance data: " << lBalFileName_m << endl;
    }

    // Define file name.
    fn += string(".h5");

    H5file_m = H5OpenFile(fn.c_str(), H5_FLUSH_STEP | H5_O_RDWR, Ippl::getComm());

    *gmsg << "Will append to " << fn << endl;

    if(H5file_m == (void*)H5_ERR) {
        throw OpalException("DataSink::DataSink(int)",
                            "failed to open h5 file '" + fn + "'");
    }

    int numStepsInFile = H5GetNumSteps(H5file_m);

    *gmsg << "numStepsInFile " << numStepsInFile << endl;

    if(restartStep == -1) {
        restartStep = numStepsInFile;
    }

    // Use same dump frequency.
    h5_int64_t dumpfreq = 0;
    rc = H5ReadFileAttribInt64(H5file_m, "dump frequency", &dumpfreq);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    OpalData::getInstance()->setRestartDumpFreq(dumpfreq);

    // Set current record/time step to restart step.
    H5call_m = restartStep;

}

DataSink::~DataSink() {
    h5_int64_t rc;
    if(H5file_m) {
        rc = H5CloseFile(H5file_m);
	H5file_m = NULL;
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    }
    Ippl::Comm->barrier();
}

void DataSink::writeH5FileAttributes() {

    if (!doHDF5_m) return;

    h5_int64_t rc;
    /// Function steps:

    /// Write file attributes to describe phase space to H5 file.
    stringstream OPAL_version;
    OPAL_version << PACKAGE_NAME << " " << PACKAGE_VERSION << " svn rev. " << GIT_VERSION;
    rc = H5WriteFileAttribString(H5file_m, "OPAL_version", OPAL_version.str().c_str());
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteFileAttribString(H5file_m, "tUnit", "s");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteFileAttribString(H5file_m, "xUnit", "m");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteFileAttribString(H5file_m, "yUnit", "m");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteFileAttribString(H5file_m, "zUnit", "m");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteFileAttribString(H5file_m, "pxUnit", "#beta#gamma");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteFileAttribString(H5file_m, "pyUnit", "#beta#gamma");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteFileAttribString(H5file_m, "pzUnit", "#beta#gamma");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteFileAttribString(H5file_m, "qUnit", "Cb");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteFileAttribString(H5file_m, "idUnit", "1");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteFileAttribString(H5file_m, "ptype", "1");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteFileAttribString(H5file_m, "lastsection", "1");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteFileAttribString(H5file_m, "SPOSUnit", "m");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteFileAttribString(H5file_m, "TIMEUnit", "s");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteFileAttribString(H5file_m, "#gammaUnit", "1");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteFileAttribString(H5file_m, "ENERGYUnit", "MeV");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteFileAttribString(H5file_m, "#varepsilonUnit", "m rad");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteFileAttribString(H5file_m, "#varepsilonrUnit", "m rad");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteFileAttribString(H5file_m, "#varepsilon-geomUnit", "m rad");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteFileAttribString(H5file_m, "#sigmaUnit", " ");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteFileAttribString(H5file_m, "RMSXUnit", "m");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteFileAttribString(H5file_m, "RMSRUnit", "m");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteFileAttribString(H5file_m, "RMSPUnit", "#beta#gamma");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteFileAttribString(H5file_m, "maxdEUnit", "MeV");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteFileAttribString(H5file_m, "max#phiUnit", "deg");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteFileAttribString(H5file_m, "phizUnit", "deg");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteFileAttribString(H5file_m, "dEUnit", "MeV");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteFileAttribString(H5file_m, "mpart", "GeV");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteFileAttribString(H5file_m, "qi", "C");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    /// Write file attributes to describe fields of head/ref particle/tail.
    rc = H5WriteFileAttribString(H5file_m, "spos-headUnit", "m");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteFileAttribString(H5file_m, "E-headUnit", "MV/m");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteFileAttribString(H5file_m, "B-headUnit", "T");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteFileAttribString(H5file_m, "spos-refUnit", "m");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteFileAttribString(H5file_m, "E-refUnit", "MV/m");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteFileAttribString(H5file_m, "B-refUnit", "T");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteFileAttribString(H5file_m, "spos-tailUnit", "m");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteFileAttribString(H5file_m, "E-tailUnit", "MV/m");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteFileAttribString(H5file_m, "B-tailUnit", "T");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteFileAttribString(H5file_m, "StepUnit", " ");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteFileAttribString(H5file_m, "LocalTrackStepUnit", " ");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteFileAttribString(H5file_m, "GlobalTrackStepUnit", " ");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteFileAttribString(H5file_m, "NumBunchUnit", " ");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteFileAttribString(H5file_m, "NumPartUnit", " ");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteFileAttribString(H5file_m, "RefPartRUnit", "m");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteFileAttribString(H5file_m, "RefPartPUnit", "#beta#gamma");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteFileAttribString(H5file_m, "SteptoLastInjUnit", "");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    /// Write file dump frequency.
    h5_int64_t dumpfreq = Options::psDumpFreq;
    rc = H5WriteFileAttribInt64(H5file_m, "dump frequency", &dumpfreq, 1);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    /// Write global phase change
    h5_float64_t dphi = OpalData::getInstance()->getGlobalPhaseShift();
    rc = H5WriteFileAttribFloat64(H5file_m, "dPhiGlobal", &dphi, 1);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    /// Reset first write flag.
    firstWriteH5part_m = false;

    /*
      stringstream inputFileContent;
      hsize_t ContentLength;
      hsize_t write_length;
      hsize_t length;
      hsize_t start = 0;
      hsize_t dmax = H5S_UNLIMITED;

      herr_t herr;

      hid_t group_id;
      hid_t shape;
      hid_t dataset_id;
      hid_t diskshape;
      hid_t memshape;

      char group_name[] = "INPUT";
      char dataset_name[] = "InputFile";

      char *FileContent = NULL;

      if(H5file_m->timegroup >= 0) {
      herr = H5Gclose(H5file_m->timegroup);
      H5file_m->timegroup = -1;
      }

      if(Ippl::myNode() == 0) {
      struct stat st;
      off_t fsize;
      if(stat(OpalData::getInstance()->getInputFn().c_str(), &st) == 0) {
      fsize = st.st_size;
      }
      ContentLength = fsize / sizeof(char);
      FileContent = new char[ContentLength];

      filebuf inputFileBuffer;
      inputFileBuffer.open(OpalData::getInstance()->getInputFn().c_str(), ios::in);
      istream inputFile(&inputFileBuffer);

      inputFile.get(FileContent, ContentLength, '\0');

      inputFileBuffer.close();
      write_length = ContentLength;

      } else {
      FileContent = new char[1];
      //        FileContent[0] = '.';
      write_length = 0;
      }
      //     long n = static_cast<long>(floor( 0.5 + ContentLength / Ippl::getNodes() ) );
      //     int N = n * Ippl::getNodes() - ContentLength;

      //     int signN = N > 0 ? 1 : -1;
      //     if (Ippl::myNode() < signN * N) {
      //         length = n - signN;
      //         start = Ippl::myNode() * length;
      //     } else {
      //         length = n;
      //         start = Ippl::myNode() * length - N;
      //     }

      MPI_Bcast(&ContentLength,
      1,
      MPI_LONG_LONG_INT,
      0,
      Ippl::getComm());

      herr = H5Gget_objinfo(H5file_m->file, group_name, 1, NULL);
      if(herr >= 0) {  // there exists a group 'INPUT'
      delete[] FileContent;
      return;
      }

      group_id = H5Gcreate(H5file_m->file, group_name, 0);

      shape = H5Screate_simple(1, &ContentLength, &ContentLength);
      dataset_id = H5Dcreate(group_id,
      dataset_name,
      H5T_NATIVE_CHAR,
      shape,
      H5P_DEFAULT);
      H5Sclose(shape);

      diskshape = H5Dget_space(dataset_id);
      H5Sselect_hyperslab(diskshape,
      H5S_SELECT_SET,
      &start,
      NULL,
      &write_length,
      NULL);

      memshape = H5Screate_simple(1, &write_length, &dmax);

      herr = H5Dwrite(dataset_id,
      H5T_NATIVE_CHAR,
      memshape,
      diskshape,
      H5file_m->xfer_prop,
      FileContent);

      H5Sclose(memshape);
      H5Dclose(dataset_id);
      H5Sclose(diskshape);
      H5Gclose(group_id);

      delete[] FileContent;

    */
}

void DataSink::retriveCavityInformation(string fn) {

    h5_int64_t rc;
    h5_int64_t nAutoPhaseCavities = 0;
    rc = H5ReadFileAttribInt64(H5file_m, "nAutoPhaseCavities", &nAutoPhaseCavities);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    for(long i = 1; i <= nAutoPhaseCavities; i++) {
        stringstream is;
        is << i;
        string elName = string("Cav-") + is.str() + string("-name");
        string elVal  = string("Cav-") + is.str() + string("-value");
        char name[128];
        h5_float64_t phi = 0;

        rc = H5ReadFileAttribString(H5file_m, elName.c_str(), name);
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        rc = H5ReadFileAttribFloat64(H5file_m, elVal.c_str(), &phi);
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        OpalData::getInstance()->setMaxPhase(string(name), double(phi));
    }
}

void DataSink::storeCavityInformation() {
    if (!doHDF5_m) return;
    /// Write number of Cavities with autophase information
    h5_int64_t nAutopPhaseCavities = OpalData::getInstance()->getNumberOfMaxPhases();

    h5_int64_t rc;

    rc = H5WriteFileAttribInt64(H5file_m, "nAutoPhaseCavities", &nAutopPhaseCavities, 1);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    int elNum = 1;
    for(vector<MaxPhasesT>::iterator it = OpalData::getInstance()->getFirstMaxPhases(); it < OpalData::getInstance()->getLastMaxPhases(); it++) {
        stringstream is;
        is << elNum;
        string elName = string("Cav-") + is.str() + string("-name");
        string elVal  = string("Cav-") + is.str() + string("-value");

        h5_float64_t phi = (*it).second;
        string name   = (*it).first;

        rc = H5WriteFileAttribString(H5file_m, elName.c_str(), name.c_str());
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        rc = H5WriteFileAttribFloat64(H5file_m, elVal.c_str(), &phi, 1);
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        INFOMSG("Saved phases in the h5 file: " << elName << "->" << name <<   " --- " << elVal << " ->" << phi << endl);
        elNum++;
    }
}


int DataSink::storeFieldmaps() {

    vector<string> fieldmap_list = Fieldmap::getListFieldmapNames();
    vector<file_size_name> fieldmap_list2;
    off_t fsize;
    struct stat st;
    int Nf = 0;
    int Np = Ippl::getNodes();
    if(Np > 20) Np = 20;  // limit the number of processors that write data

    for(vector<string>::const_iterator it = fieldmap_list.begin(); it != fieldmap_list.end(); ++ it) {
        if(stat((*it).c_str(), &st) == 0) {
            fsize = st.st_size;
            ++ Nf;
        } else {
            continue;
        }
        fieldmap_list2.push_back(file_size_name((*it), fsize));
    }
    sort(fieldmap_list2.begin(), fieldmap_list2.end(), file_size_name::SortAsc);



    if(Ippl::myNode() < Np) {
        for(int i = 0; i < Nf / Np; ++ i) {
            int lid = i * Np + Ippl::myNode();
            string filename = fieldmap_list2[lid].file_name_m;
            int ContentLength = fieldmap_list2[lid].file_size_m / sizeof(char);

            filebuf inputFileBuffer;
            inputFileBuffer.open(filename.c_str(), ios::in);
            istream inputFile(&inputFileBuffer);

            std::unique_ptr<char[]> FileContent(new char[ContentLength]);

            inputFile.get(FileContent.get(), ContentLength, '\0');

            inputFileBuffer.close();
        }
    }
    //group_id = H5Gcreate ( H5file_m->file, group_name, 0 );
    return Nf;

    return 1;
}

void DataSink::writePhaseSpace(PartBunch &beam, Vector_t FDext[], double sposHead, double sposRef, double sposTail) {

    if (!doHDF5_m) return;

    h5_int64_t rc;
    /// Function steps:

    /// Start timer.
    IpplTimings::startTimer(H5PartTimer_m);

    /// Calculate beam statistical parameters etc. Put them in the right format
    /// for H5 file write.
    beam.calcBeamParameters();

    double               actPos   = beam.get_sPos();
    double               t        = beam.getT();
    Vektor< double, 3 >  rmin     = beam.get_origin();
    Vektor< double, 3 >  rmax     = beam.get_maxExtend();
    Vektor< double, 3 >  centroid = beam.get_centroid();

    size_t nLoc                   = beam.getLocalNum();

    Vektor< double, 3 >  maxP(0.0);
    Vektor< double, 3 >  minP(0.0);

    Vektor<double, 3 > xsigma = beam.get_rrms();
    Vektor<double, 3 > psigma = beam.get_prms();
    Vektor<double, 3 > vareps = beam.get_norm_emit();
    Vektor<double, 3 > geomvareps = beam.get_emit();
    Vektor<double, 3 > RefPartR = beam.RefPart_R;
    Vektor<double, 3 > RefPartP = beam.RefPart_P;
    Vektor<double, 3>  pmean = beam.get_pmean();

    double meanEnergy = beam.get_meanEnergy();
    double energySpread = beam.getdE();

    double sigma = ((xsigma[0] * xsigma[0]) + (xsigma[1] * xsigma[1])) /
        (2.0 * beam.get_gamma() * 17.0e3 * ((geomvareps[0] * geomvareps[0]) + (geomvareps[1] * geomvareps[1])));

    beam.get_PBounds(minP, maxP);

    h5_int64_t localTrackStep = (h5_int64_t)beam.getLocalTrackStep();
    h5_int64_t globalTrackStep = (h5_int64_t)beam.getGlobalTrackStep();

    std::unique_ptr<char[]> varray(new char[(nLoc)*sizeof(double)]);
    double *farray = reinterpret_cast<double *>(varray.get());
    h5_int64_t *larray = reinterpret_cast<h5_int64_t *>(varray.get());

    ///Get the particle decomposition from all the compute nodes.
    std::unique_ptr<size_t[]> locN(new size_t[Ippl::getNodes()]);
    std::unique_ptr<size_t[]> globN(new size_t[Ippl::getNodes()]);

    for(int i = 0; i < Ippl::getNodes(); i++) {
        globN[i] = locN[i] = 0;
    }
    locN[Ippl::myNode()] = nLoc;
    reduce(locN.get(), locN.get() + Ippl::getNodes(), globN.get(), OpAddAssign());

    /* ------------------------------------------------------------------------ */

    rc = H5SetStep(H5file_m, H5call_m);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5PartSetNumParticles(H5file_m, nLoc);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribFloat64(H5file_m, "RMSX", (h5_float64_t *)&xsigma, 3);    //sigma
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribFloat64(H5file_m, "#varepsilon", (h5_float64_t *)&vareps, 3);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribFloat64(H5file_m, "RMSP", (h5_float64_t *)&psigma, 3);    //sigma
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribFloat64(H5file_m, "MEANP", (h5_float64_t *)&pmean, 3);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribFloat64(H5file_m, "centroid", (h5_float64_t *)&centroid, 3);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribFloat64(H5file_m, "maxX", (h5_float64_t *)&rmax, 3);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribFloat64(H5file_m, "minX", (h5_float64_t *)&rmin, 3);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribFloat64(H5file_m, "maxP", (h5_float64_t *)&maxP, 3);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribFloat64(H5file_m, "minP", (h5_float64_t *)&minP, 3);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribFloat64(H5file_m, "#varepsilon-geom", (h5_float64_t *)&geomvareps, 3);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribFloat64(H5file_m, "B-ref", (h5_float64_t *)&FDext[2], 3);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribFloat64(H5file_m, "E-ref", (h5_float64_t *)&FDext[3], 3);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribFloat64(H5file_m, "B-head", (h5_float64_t *)&FDext[0], 3);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribFloat64(H5file_m, "E-head", (h5_float64_t *)&FDext[1], 3);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribFloat64(H5file_m, "B-tail", (h5_float64_t *)&FDext[4], 3);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribFloat64(H5file_m, "E-tail", (h5_float64_t *)&FDext[5], 3);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribFloat64(H5file_m, "RefPartR", (h5_float64_t *)&RefPartR, 3);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribFloat64(H5file_m, "RefPartP", (h5_float64_t *)&RefPartP, 3);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    // scalar values

    rc = H5WriteStepAttribFloat64(H5file_m, "SPOS",     &actPos, 1);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribInt64(H5file_m, "Step",       &H5call_m, 1);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribFloat64(H5file_m, "#sigma",   &sigma, 1);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribFloat64(H5file_m, "TIME",     &t, 1);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribFloat64(H5file_m, "ENERGY",   &meanEnergy, 1);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribFloat64(H5file_m, "dE",       &energySpread, 1);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribInt64(H5file_m, "LocalTrackStep",        &localTrackStep, 1);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribInt64(H5file_m, "GlobalTrackStep",        &globalTrackStep, 1);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);


    /// Write particle mass and charge per particle. (Consider making these file attributes.)
    double mpart = 1.0e-9 * beam.getM();
    double     Q = beam.getCharge();

    rc = H5WriteStepAttribFloat64(H5file_m, "mpart", &mpart, 1);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribFloat64(H5file_m, "Q", &Q, 1);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribFloat64(H5file_m, "spos-head", &sposHead, 1);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteStepAttribFloat64(H5file_m, "spos-ref",  &sposRef, 1);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteStepAttribFloat64(H5file_m, "spos-tail", &sposTail, 1);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    /*
      Attributes originally not in OPAL-t

    */

    h5_int64_t numBunch = 1;
    h5_int64_t SteptoLastInj = 0;

    rc = H5WriteStepAttribInt64(H5file_m, "NumBunch",         &numBunch, 1);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribInt64(H5file_m, "SteptoLastInj",    &SteptoLastInj, 1);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    /*
      Done attributes originally not in OPAL-t

    */

    setOPALt();

    /// Write beam phase space.
    for(size_t i = 0; i < nLoc; i++)
        farray[i] =  beam.R[i](0);
    rc = H5PartWriteDataFloat64(H5file_m, "x", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    for(size_t i = 0; i < nLoc; i++)
        farray[i] =  beam.R[i](1);
    rc = H5PartWriteDataFloat64(H5file_m, "y", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    for(size_t i = 0; i < nLoc; i++)
        farray[i] =  beam.R[i](2);
    rc = H5PartWriteDataFloat64(H5file_m, "z", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    for(size_t i = 0; i < nLoc; i++)
        farray[i] =  beam.P[i](0);
    rc = H5PartWriteDataFloat64(H5file_m, "px", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    for(size_t i = 0; i < nLoc; i++)
        farray[i] =  beam.P[i](1);
    rc = H5PartWriteDataFloat64(H5file_m, "py", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    for(size_t i = 0; i < nLoc; i++)
        farray[i] =  beam.P[i](2);
    rc = H5PartWriteDataFloat64(H5file_m, "pz", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    for(size_t i = 0; i < nLoc; i++)
        farray[i] =  beam.Q[i];
    rc = H5PartWriteDataFloat64(H5file_m, "q", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    for(size_t i = 0; i < nLoc; i++)
        larray[i] =  beam.ID[i];
    rc = H5PartWriteDataInt64(H5file_m, "id", larray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    for(size_t i = 0; i < nLoc; i++)
        larray[i] = (h5_int64_t) beam.PType[i];
    rc = H5PartWriteDataInt64(H5file_m, "ptype", larray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    for(size_t i = 0; i < nLoc; i++)
        larray[i] =  beam.LastSection[i];
    rc = H5PartWriteDataInt64(H5file_m, "lastsection", larray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    // Need this because rc = H5root does not yet work with
    // a vaiable number of data

    if(Options::ebDump) {
        for(size_t i = 0; i < nLoc; i++)
            farray[i] =  beam.Ef[i](0);
        rc = H5PartWriteDataFloat64(H5file_m, "Ex", farray);
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

        for(size_t i = 0; i < nLoc; i++)
            farray[i] =  beam.Ef[i](1);
        rc = H5PartWriteDataFloat64(H5file_m, "Ey", farray);
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

        for(size_t i = 0; i < nLoc; i++)
            farray[i] =  beam.Ef[i](2);
        rc = H5PartWriteDataFloat64(H5file_m, "Ez", farray);
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

        for(size_t i = 0; i < nLoc; i++)
            farray[i] =  beam.Bf[i](0);
        rc = H5PartWriteDataFloat64(H5file_m, "Bx", farray);
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

        for(size_t i = 0; i < nLoc; i++)
            farray[i] =  beam.Bf[i](1);
        rc = H5PartWriteDataFloat64(H5file_m, "By", farray);
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

        for(size_t i = 0; i < nLoc; i++)
            farray[i] =  beam.Bf[i](2);
        rc = H5PartWriteDataFloat64(H5file_m, "Bz", farray);
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    }

    /// Write space charge field map if asked for.
    if(Options::rhoDump) {

        IpplTimings::startTimer(H5BlockTimer_m);

        NDIndex<3> idx = beam.getFieldLayout().getLocalNDIndex();
        NDIndex<3> elem;
        h5_err_t herr = H5Block3dSetView(
                                         H5file_m,
                                         idx[0].min(), idx[0].max(),
                                         idx[1].min(), idx[1].max(),
                                         idx[2].min(), idx[2].max());

        if(herr < 0)
            *gmsg << "H5Block3dSetView err " << herr << endl;

        std::unique_ptr<h5_float64_t[]> data(new h5_float64_t[(idx[0].max() + 1)  * (idx[1].max() + 1) * (idx[2].max() + 1)]);

        int ii = 0;
        // h5block uses the fortran convention of storing data:
        // INTEGER, DIMENSION(2,3) :: a
        // => {a(1,1), a(2,1), a(1,2), a(2,2), a(1,3), a(2,3)}
        for(int i = idx[2].min(); i <= idx[2].max(); ++i) {
            for(int j = idx[1].min(); j <= idx[1].max(); ++j) {
                for(int k = idx[0].min(); k <= idx[0].max(); ++k) {
                    data[ii] = beam.getRho(k, j, i);
                    ii++;
                }
            }
        }
        herr = H5Block3dWriteScalarFieldFloat64(H5file_m, "rho", data.get());
        if(herr < 0)
            *gmsg << "H5Block3dWriteScalarField err " << herr << endl;

        /// Need this to align particles and fields when writing space charge map.
        herr = H5Block3dSetFieldOrigin(H5file_m, "rho",
                                       (h5_float64_t)beam.get_origin()(0),
                                       (h5_float64_t)beam.get_origin()(1),
                                       (h5_float64_t)beam.get_origin()(2));

        herr = H5Block3dSetFieldSpacing(H5file_m, "rho",
                                        (h5_float64_t)beam.get_hr()(0),
                                        (h5_float64_t)beam.get_hr()(1),
                                        (h5_float64_t)beam.get_hr()(2));

    }
    /// Step record/time step index.
    H5call_m++;

    /// %Stop timer.
    IpplTimings::stopTimer(H5PartTimer_m);
    if(Options::rhoDump)
        IpplTimings::stopTimer(H5BlockTimer_m);
}



int DataSink::writePhaseSpace_cycl(PartBunch &beam, Vector_t FDext[], double meanEnergy,
                                   double refPr, double refPt, double refPz,
                                   double refR, double refTheta, double refZ,
                                   double azimuth, double elevation, bool local) {

    if (!doHDF5_m) return -1;
    //if (beam.getLocalNum() == 0) return -1; //TEMP for testing -DW

    h5_int64_t rc;
    IpplTimings::startTimer(H5PartTimer_m);

    beam.calcBeamParameters_cycl();

    double               t        = beam.getT();

    Vektor< double, 3 >  rmin     = beam.get_origin();
    Vektor< double, 3 >  rmax     = beam.get_maxExtend();
    Vektor< double, 3 >  centroid = beam.get_centroid();
    size_t nLoc                   = beam.getLocalNum();
    Vektor<double, 3 > xsigma = beam.get_rrms();
    Vektor<double, 3 > psigma = beam.get_prms();
    Vektor<double, 3 > geomvareps = beam.get_emit();
    Vektor<double, 3 > vareps = beam.get_norm_emit();

    Vektor<double, 3 > RefPartR = beam.RefPart_R;
    Vektor<double, 3 > RefPartP = beam.RefPart_P;

    double energySpread = beam.getdE();

    double sigma = ((xsigma[0] * xsigma[0]) + (xsigma[1] * xsigma[1])) /
        (2.0 * beam.get_gamma() * 17.0e3 * ((geomvareps[0] * geomvareps[0]) + (geomvareps[1] * geomvareps[1])));

    Vektor< double, 3 >  maxP(0.0);
    Vektor< double, 3 >  minP(0.0);
    beam.get_PBounds(minP, maxP);

    std::unique_ptr<char[]> varray(new char[(nLoc)*sizeof(double)]);
    double *farray = reinterpret_cast<double *>(varray.get());
    h5_int64_t *larray = reinterpret_cast<h5_int64_t *>(varray.get());

    double  pathLength = beam.getLPath();

    h5_int64_t localTrackStep = (h5_int64_t)beam.getLocalTrackStep();
    h5_int64_t globalTrackStep = (h5_int64_t)beam.getGlobalTrackStep();
    h5_int64_t numBunch = (h5_int64_t)beam.getNumBunch();
    h5_int64_t SteptoLastInj = (h5_int64_t)beam.getSteptoLastInj();

    h5_int64_t localFrame = 0;
    if (local) localFrame = 1;

    ///Get the particle decomposition from all the compute nodes.
    std::unique_ptr<size_t[]> locN(new size_t[Ippl::getNodes()]);
    std::unique_ptr<size_t[]> globN(new size_t[Ippl::getNodes()]);

    for(int i = 0; i < Ippl::getNodes(); i++) {
        globN[i] = locN[i] = 0;
    }
    locN[Ippl::myNode()] = nLoc;
    reduce(locN.get(), locN.get() + Ippl::getNodes(), globN.get(), OpAddAssign());

    /* ------------------------------------------------------------------------ */

    rc = H5SetStep(H5file_m, H5call_m);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5PartSetNumParticles(H5file_m, nLoc);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribInt64(H5file_m, "LocalTrackStep", &localTrackStep, 1);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribInt64(H5file_m, "GlobalTrackStep", &globalTrackStep, 1);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribFloat64(H5file_m, "SPOS",     &pathLength, 1);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribInt64(H5file_m, "Step",       &H5call_m, 1);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribFloat64(H5file_m, "#sigma",   &sigma, 1);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribFloat64(H5file_m, "RMSX", (h5_float64_t *)&xsigma, 3);    //sigma
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribFloat64(H5file_m, "RMSP", (h5_float64_t *)&psigma, 3);    //sigma
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribFloat64(H5file_m, "maxX", (h5_float64_t *)&rmax, 3);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribFloat64(H5file_m, "minX", (h5_float64_t *)&rmin, 3);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribFloat64(H5file_m, "maxP", (h5_float64_t *)&maxP, 3);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribFloat64(H5file_m, "minP", (h5_float64_t *)&minP, 3);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribFloat64(H5file_m, "centroid", (h5_float64_t *)&centroid, 3);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribFloat64(H5file_m, "TIME",     &t, 1);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribFloat64(H5file_m, "ENERGY",   &meanEnergy, 1);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribFloat64(H5file_m, "#varepsilon", (h5_float64_t *)&vareps, 3);     //unnormalized
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribFloat64(H5file_m, "#varepsilon-geom", (h5_float64_t *)&geomvareps, 3); //normalized
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribFloat64(H5file_m, "dE",       &energySpread, 1);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribFloat64(H5file_m, "B-ref", (h5_float64_t *)&FDext[0], 3);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribFloat64(H5file_m, "E-ref", (h5_float64_t *)&FDext[1], 3);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribFloat64(H5file_m, "B-head", (h5_float64_t *)&FDext[0], 3);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribFloat64(H5file_m, "E-head", (h5_float64_t *)&FDext[1], 3);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribFloat64(H5file_m, "B-tail", (h5_float64_t *)&FDext[4], 3);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribFloat64(H5file_m, "E-tail", (h5_float64_t *)&FDext[5], 3);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribFloat64(H5file_m, "RefPartR", (h5_float64_t *)&RefPartR, 3);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribFloat64(H5file_m, "RefPartP", (h5_float64_t *)&RefPartP, 3);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribFloat64(H5file_m, "REFPR", (h5_float64_t *)&refPr, 1);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribFloat64(H5file_m, "REFPT", (h5_float64_t *)&refPt, 1);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribFloat64(H5file_m, "REFPZ", (h5_float64_t *)&refPz, 1);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribFloat64(H5file_m, "REFR", (h5_float64_t *)&refR, 1);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribFloat64(H5file_m, "REFTHETA", (h5_float64_t *)&refTheta, 1);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribFloat64(H5file_m, "REFZ", (h5_float64_t *)&refZ, 1);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribFloat64(H5file_m, "AZIMUTH", (h5_float64_t *)&azimuth, 1);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribFloat64(H5file_m, "ELEVATION", (h5_float64_t *)&elevation, 1);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribInt64(H5file_m, "LOCAL", &localFrame, 1);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    /// Write particle mass and charge per particle. (Consider making these file attributes.)
    double mpart = 1.0e-9 * beam.getM();
    double     Q = beam.getCharge();

    rc = H5WriteStepAttribFloat64(H5file_m, "mpart", &mpart, 1);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribFloat64(H5file_m, "Q", &Q, 1);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);


    rc = H5WriteStepAttribInt64(H5file_m, "NumBunch",         &numBunch, 1);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribInt64(H5file_m, "SteptoLastInj",    &SteptoLastInj, 1);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);


    /*
      Attributes originally not in OPAL-cycl
    */

    double sposHead = 0.0;
    double sposRef = 0.0;
    double sposTail = 0.0;
    rc = H5WriteStepAttribFloat64(H5file_m, "spos-head", &sposHead, 1);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteStepAttribFloat64(H5file_m, "spos-ref",  &sposRef, 1);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteStepAttribFloat64(H5file_m, "spos-tail", &sposTail, 1);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    /*
      Done attributes originally not in OPAL-cycl
    */

    setOPALcycl();

    for(size_t i = 0; i < nLoc; i++)
        farray[i] =  beam.R[i](0);
    rc = H5PartWriteDataFloat64(H5file_m, "x", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    for(size_t i = 0; i < nLoc; i++)
        farray[i] =  beam.R[i](1);
    rc = H5PartWriteDataFloat64(H5file_m, "y", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    for(size_t i = 0; i < nLoc; i++)
        farray[i] =  beam.R[i](2);
    rc = H5PartWriteDataFloat64(H5file_m, "z", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    for(size_t i = 0; i < nLoc; i++)
        farray[i] =  beam.P[i](0);
    rc = H5PartWriteDataFloat64(H5file_m, "px", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    for(size_t i = 0; i < nLoc; i++)
        farray[i] =  beam.P[i](1);
    rc = H5PartWriteDataFloat64(H5file_m, "py", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    for(size_t i = 0; i < nLoc; i++)
        farray[i] =  beam.P[i](2);
    rc = H5PartWriteDataFloat64(H5file_m, "pz", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    for(size_t i = 0; i < nLoc; i++)
        farray[i] =  beam.Q[i];
    rc = H5PartWriteDataFloat64(H5file_m, "q", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    for(size_t i = 0; i < nLoc; i++)
        larray[i] =  beam.ID[i];
    rc = H5PartWriteDataInt64(H5file_m, "id", larray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    for(size_t i = 0; i < nLoc; i++)
        larray[i] =  beam.PType[i];
    rc = H5PartWriteDataInt64(H5file_m, "ptype", larray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    for(size_t i = 0; i < nLoc; i++)
        larray[i] =  beam.LastSection[i];
    rc = H5PartWriteDataInt64(H5file_m, "lastsection", larray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    // Need this because rc = H5root does not yet work with
    // avaiable number of data
    if(Options::ebDump) {

        for(size_t i = 0; i < nLoc; i++)
            farray[i] =  beam.Ef[i](0);
        rc = H5PartWriteDataFloat64(H5file_m, "Ex", farray);
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

        for(size_t i = 0; i < nLoc; i++)
            farray[i] =  beam.Ef[i](1);
        rc = H5PartWriteDataFloat64(H5file_m, "Ey", farray);
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

        for(size_t i = 0; i < nLoc; i++)
            farray[i] =  beam.Ef[i](2);
        rc = H5PartWriteDataFloat64(H5file_m, "Ez", farray);
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

        for(size_t i = 0; i < nLoc; i++)
            farray[i] =  beam.Bf[i](0);
        rc = H5PartWriteDataFloat64(H5file_m, "Bx", farray);
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

        for(size_t i = 0; i < nLoc; i++)
            farray[i] =  beam.Bf[i](1);
        rc = H5PartWriteDataFloat64(H5file_m, "By", farray);
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

        for(size_t i = 0; i < nLoc; i++)
            farray[i] =  beam.Bf[i](2);
        rc = H5PartWriteDataFloat64(H5file_m, "Bz", farray);
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    }

    /// Write space charge field map if asked for.
    if(Options::rhoDump) {

        IpplTimings::startTimer(H5BlockTimer_m);

        NDIndex<3> idx = beam.getFieldLayout().getLocalNDIndex();
        NDIndex<3> elem;
        h5_err_t herr = H5Block3dSetView(
                                         H5file_m,
                                         idx[0].min(), idx[0].max(),
                                         idx[1].min(), idx[1].max(),
                                         idx[2].min(), idx[2].max());

        if(herr < 0)
            *gmsg << "H5Block3dSetView err " << herr << endl;

        std::unique_ptr<h5_float64_t[]> data(new h5_float64_t[(idx[0].max() + 1)  * (idx[1].max() + 1) * (idx[2].max() + 1)]);

        int ii = 0;
        // h5block uses the fortran convention of storing data:
        // INTEGER, DIMENSION(2,3) :: a
        // => {a(1,1), a(2,1), a(1,2), a(2,2), a(1,3), a(2,3)}
        for(int i = idx[2].min(); i <= idx[2].max(); ++i) {
            for(int j = idx[1].min(); j <= idx[1].max(); ++j) {
                for(int k = idx[0].min(); k <= idx[0].max(); ++k) {
                    data[ii] = beam.getRho(k, j, i);
                    ii++;
                }
            }
        }
        herr = H5Block3dWriteScalarFieldFloat64(H5file_m, "rho", data.get());
        if(herr < 0)
            *gmsg << "H5Block3dWriteScalarField err " << herr << endl;

        /// Need this to align particles and fields when writing space charge map.
        herr = H5Block3dSetFieldOrigin(H5file_m, "rho",
                                       (h5_float64_t)beam.get_origin()(0),
                                       (h5_float64_t)beam.get_origin()(1),
                                       (h5_float64_t)beam.get_origin()(2));

        herr = H5Block3dSetFieldSpacing(H5file_m, "rho",
                                        (h5_float64_t)beam.get_hr()(0),
                                        (h5_float64_t)beam.get_hr()(1),
                                        (h5_float64_t)beam.get_hr()(2));

    }
    /// Step record/time step index.
    H5call_m++;

    IpplTimings::stopTimer(H5PartTimer_m);
    if(Options::rhoDump)
        IpplTimings::stopTimer(H5BlockTimer_m);
    return (int)(H5call_m - 1);
}

void DataSink::writePhaseSpaceEnvelope(EnvelopeBunch &beam, Vector_t FDext[], double sposHead, double sposRef, double sposTail) {

    if (!doHDF5_m) return;

    h5_int64_t rc;
    /// Function steps:

    /// Start timer.
    IpplTimings::startTimer(H5PartTimer_m);

    /// Calculate beam statistical parameters etc. Put them in the right format
    /// for H5 file write.
    beam.calcBeamParameters();

    /*
      "dEdt", "double","MeV/ps", "beam energy");
      "dE",   "double","MeV",    "rms energy spread");
      "Imax", "double","A",      "max bunch current");
      "Irms", "double","A",      "rms bunch current");

      "tau",  "double","ps",     "rms bunch length");
      "Rx",   "double","m",      "beam radius x");
      "Ry",   "double","m",      "beam radius y");
      "Px",   "double","mrad",   "beam divergence x");
      "Py",   "double","mrad",   "beam divergence y");
    */

    //TODO:
    Vektor<double, 3> RefPartR = Vector_t(0.0); //beam.RefPart_R;
    Vektor<double, 3> RefPartP = Vector_t(0.0); //beam.RefPart_P;
    Vektor<double, 3> centroid = Vector_t(0.0);
    Vektor<double, 3> geomvareps = Vector_t(0.0); //beam.emtn();

    double actPos   = beam.get_sPos();
    double t        = beam.getT();

    size_t nLoc = beam.getLocalNum();

    Vektor<double, 3> xsigma = beam.sigmax();
    Vektor<double, 3> psigma = beam.sigmap();
    Vektor<double, 3> vareps = beam.get_norm_emit();
    // min beam radius
    Vektor<double, 3> rmin     = beam.minX();
    // max beam radius
    Vektor<double, 3> rmax     = beam.maxX();
    Vektor<double, 3> maxP = beam.maxP();
    Vektor<double, 3> minP = beam.minP();

    //in MeV
    double meanEnergy = beam.get_meanEnergy() * 1e-6;

    double sigma = ((xsigma[0] * xsigma[0]) + (xsigma[1] * xsigma[1])) /
        (2.0 * beam.get_gamma() * 17.0e3 * ((geomvareps[0] * geomvareps[0]) + (geomvareps[1] * geomvareps[1])));

    //beam.get_PBounds(minP,maxP);

    std::unique_ptr<char[]> varray(new char[(nLoc)*sizeof(double)]);
    double *farray = reinterpret_cast<double *>(varray.get());

    ///Get the particle decomposition from all the compute nodes.
    std::unique_ptr<size_t[]> locN(new size_t[Ippl::getNodes()]);
    std::unique_ptr<size_t[]> globN(new size_t[Ippl::getNodes()]);

    for(int i = 0; i < Ippl::getNodes(); i++) {
        globN[i] = locN[i] = 0;
    }
    locN[Ippl::myNode()] = nLoc;
    reduce(locN.get(), locN.get() + Ippl::getNodes(), globN.get(), OpAddAssign());

    /// Set current record/time step.
    rc = H5SetStep(H5file_m, H5call_m);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5PartSetNumParticles(H5file_m, nLoc);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    /// Write statistical data.
    rc = H5WriteStepAttribFloat64(H5file_m, "SPOS",     &actPos, 1);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteStepAttribFloat64(H5file_m, "#sigma",   &sigma, 1);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteStepAttribFloat64(H5file_m, "RMSX", (h5_float64_t *)&xsigma, 3);    //sigma
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteStepAttribFloat64(H5file_m, "RMSP", (h5_float64_t *)&psigma, 3);    //sigma
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteStepAttribFloat64(H5file_m, "maxX", (h5_float64_t *)&rmax, 3);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteStepAttribFloat64(H5file_m, "minX", (h5_float64_t *)&rmin, 3);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteStepAttribFloat64(H5file_m, "maxP", (h5_float64_t *)&maxP, 3);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteStepAttribFloat64(H5file_m, "minP", (h5_float64_t *)&minP, 3);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteStepAttribFloat64(H5file_m, "centroid", (h5_float64_t *)&centroid, 3);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteStepAttribFloat64(H5file_m, "TIME",     &t, 1);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteStepAttribFloat64(H5file_m, "ENERGY",   &meanEnergy, 1);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteStepAttribFloat64(H5file_m, "RefPartR", (h5_float64_t *)&RefPartR, 3);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteStepAttribFloat64(H5file_m, "RefPartP", (h5_float64_t *)&RefPartP, 3);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    /// Write head/reference particle/tail field information.

    rc = H5WriteStepAttribFloat64(H5file_m, "B-head", (h5_float64_t *)&FDext[0], 3);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteStepAttribFloat64(H5file_m, "E-head", (h5_float64_t *)&FDext[1], 3);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteStepAttribFloat64(H5file_m, "B-ref", (h5_float64_t *)&FDext[2], 3);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteStepAttribFloat64(H5file_m, "E-ref", (h5_float64_t *)&FDext[3], 3);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteStepAttribFloat64(H5file_m, "B-tail", (h5_float64_t *)&FDext[4], 3);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteStepAttribFloat64(H5file_m, "E-tail", (h5_float64_t *)&FDext[5], 3);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribFloat64(H5file_m, "spos-head", &sposHead, 1);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteStepAttribFloat64(H5file_m, "spos-ref",  &sposRef, 1);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteStepAttribFloat64(H5file_m, "spos-tail", &sposTail, 1);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    /// Write number of compute nodes.
    //    H5WriteStepAttribInt64(H5file_m, "nloc", globN, Ippl::getNodes());

    /// Write particle mass and charge per particle. (Consider making these
    /// file attributes.)
    double mpart = 1.0e-9 * beam.getM();
    double    qi = beam.getChargePerParticle();

    rc = H5WriteStepAttribFloat64(H5file_m, "mpart", &mpart, 1);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteStepAttribFloat64(H5file_m, "qi", &qi, 1);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    /// Write normalized emittance.
    rc = H5WriteStepAttribFloat64(H5file_m, "#varepsilon", (h5_float64_t *)&vareps, 3);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    /// Write geometric emittance.
    rc = H5WriteStepAttribFloat64(H5file_m, "#varepsilon-geom", (h5_float64_t *)&geomvareps, 3);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    /// Write beam phase space.
    for(size_t i = 0; i < nLoc; i++)
        farray[i] =  beam.getX(i);
    rc = H5PartWriteDataFloat64(H5file_m, "x", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    for(size_t i = 0; i < nLoc; i++)
        farray[i] =  beam.getY(i);
    rc = H5PartWriteDataFloat64(H5file_m, "y", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    for(size_t i = 0; i < nLoc; i++)
        farray[i] =  beam.getZ(i);
    rc = H5PartWriteDataFloat64(H5file_m, "z", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    for(size_t i = 0; i < nLoc; i++)
        farray[i] =  beam.getPx(i);
    rc = H5PartWriteDataFloat64(H5file_m, "px", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    for(size_t i = 0; i < nLoc; i++)
        farray[i] =  beam.getPy(i);
    rc = H5PartWriteDataFloat64(H5file_m, "py", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    for(size_t i = 0; i < nLoc; i++)
        farray[i] =  beam.getPz(i);
    rc = H5PartWriteDataFloat64(H5file_m, "pz", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    for(size_t i = 0; i < nLoc; i++)
        farray[i] =  beam.getBeta(i);
    rc = H5PartWriteDataFloat64(H5file_m, "beta", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    for(size_t i = 0; i < nLoc; i++)
        farray[i] =  beam.getX0(i);
    rc = H5PartWriteDataFloat64(H5file_m, "X0", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    for(size_t i = 0; i < nLoc; i++)
        farray[i] =  beam.getPx0(i);
    rc = H5PartWriteDataFloat64(H5file_m, "pX0", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    for(size_t i = 0; i < nLoc; i++)
        farray[i] =  beam.getY0(i);
    rc = H5PartWriteDataFloat64(H5file_m, "Y0", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    for(size_t i = 0; i < nLoc; i++)
        farray[i] =  beam.getPy0(i);
    rc = H5PartWriteDataFloat64(H5file_m, "pY0", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    /// Step record/time step index.
    H5call_m++;

    /// %Stop timer.
    IpplTimings::stopTimer(H5PartTimer_m);
}

void DataSink::stashPhaseSpaceEnvelope(EnvelopeBunch &beam, Vector_t FDext[], double sposHead, double sposRef, double sposTail) {

    if (!doHDF5_m) return;

    /// Start timer.
    IpplTimings::startTimer(H5PartTimer_m);

    /// Calculate beam statistical parameters etc. Put them in the right format
    /// for H5 file write.
    beam.calcBeamParameters();

    //TODO:
    stash_RefPartR.push_back(Vector_t(0.0)); //beam.RefPart_R;
    stash_RefPartP.push_back(Vector_t(0.0)); //beam.RefPart_P;
    stash_centroid.push_back(Vector_t(0.0));

    stash_actPos.push_back(beam.get_sPos());
    stash_t.push_back(beam.getT());

    stash_nTot.push_back(beam.getTotalNum());
    size_t nLoc = beam.getLocalNum();
    stash_nLoc.push_back(nLoc);

    stash_xsigma.push_back(beam.sigmax());
    stash_psigma.push_back(beam.sigmap());
    stash_vareps.push_back(beam.get_norm_emit());
    stash_geomvareps.push_back(beam.emtn());
    stash_rmin.push_back(beam.minX());
    stash_rmax.push_back(beam.maxX());
    stash_maxP.push_back(beam.maxP());
    stash_minP.push_back(beam.minP());

    //in MeV
    stash_meanEnergy.push_back(beam.get_meanEnergy() * 1e-6);

    //stash_sigma.push_back(((xsigma[0]*xsigma[0])+(xsigma[1]*xsigma[1])) /
    //    (2.0*beam.get_gamma()*17.0e3*((geomvareps[0]*geomvareps[0]) + (geomvareps[1]*geomvareps[1]))));

    //beam.get_PBounds(minP,maxP);

    ///Get the particle decomposition from all the compute nodes.
    std::unique_ptr<size_t[]> locN(new size_t[Ippl::getNodes()]);
    std::unique_ptr<size_t[]> globN(new size_t[Ippl::getNodes()]);

    for(int i = 0; i < Ippl::getNodes(); i++) {
        globN[i] = locN[i] = 0;
    }
    locN[Ippl::myNode()] = nLoc;
    reduce(locN.get(), locN.get() + Ippl::getNodes(), globN.get(), OpAddAssign());

    stash_sposHead.push_back(sposHead);
    stash_sposRef.push_back(sposRef);
    stash_sposTail.push_back(sposTail);
    stash_Bhead.push_back(FDext[0]);
    stash_Ehead.push_back(FDext[1]);
    stash_Bref.push_back(FDext[2]);
    stash_Eref.push_back(FDext[3]);
    stash_Btail.push_back(FDext[4]);
    stash_Etail.push_back(FDext[5]);

    /// Write particle mass and charge per particle. (Consider making these
    /// file attributes.)
    stash_mpart.push_back(1.0e-9 * beam.getM());
    stash_qi.push_back(beam.getChargePerParticle());

    /// Write beam phase space.
    //for (size_t i=0; i<nLoc;i++)
    //farray[i] =  beam.getX(i);

    //for (size_t i=0; i<nLoc;i++)
    //farray[i] =  beam.getY(i);

    //for (size_t i=0; i<nLoc;i++)
    //farray[i] =  beam.getZ(i);

    //for (size_t i=0; i<nLoc;i++)
    //farray[i] =  beam.getPx(i);

    //for (size_t i=0; i<nLoc;i++)
    //farray[i] =  beam.getPy(i);

    //for (size_t i=0; i<nLoc;i++)
    //farray[i] =  beam.getPz(i);

    /// Step record/time step index.
    H5call_m++;

    /// %Stop timer.
    IpplTimings::stopTimer(H5PartTimer_m);
}

void DataSink::dumpStashedPhaseSpaceEnvelope() {

    if (!doHDF5_m) return;

    h5_int64_t rc;
    /// Start timer.
    IpplTimings::startTimer(H5PartTimer_m);

    size_t nLoc = stash_nLoc[0];
    std::unique_ptr<char[]> varray(new char[(nLoc)*sizeof(double)]);
    double *farray = reinterpret_cast<double *>(varray.get());

    //FIXME: restart step
    for(int step = 0; step < H5call_m; step++) {

        /// Set current record/time step.
        rc = H5SetStep(H5file_m, step);
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        nLoc = stash_nLoc[step];
        rc = H5PartSetNumParticles(H5file_m, nLoc);
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

        /// Write statistical data.
        rc = H5WriteStepAttribFloat64(H5file_m, "SPOS",     &stash_actPos[step], 1);
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        //rc = H5WriteStepAttribFloat64(H5file_m,"#sigma",  &stash_sigma[step],1);
        rc = H5WriteStepAttribFloat64(H5file_m, "RMSX", (h5_float64_t *)&stash_xsigma[step], 3);    //sigma
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        rc = H5WriteStepAttribFloat64(H5file_m, "RMSP", (h5_float64_t *)&stash_psigma[step], 3);    //sigma
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        rc = H5WriteStepAttribFloat64(H5file_m, "maxX", (h5_float64_t *)&stash_rmax[step], 3);
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        rc = H5WriteStepAttribFloat64(H5file_m, "minX", (h5_float64_t *)&stash_rmin[step], 3);
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        rc = H5WriteStepAttribFloat64(H5file_m, "maxP", (h5_float64_t *)&stash_maxP[step], 3);
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        rc = H5WriteStepAttribFloat64(H5file_m, "minP", (h5_float64_t *)&stash_minP[step], 3);
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        rc = H5WriteStepAttribFloat64(H5file_m, "centroid", (h5_float64_t *)&stash_centroid[step], 3);
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        rc = H5WriteStepAttribFloat64(H5file_m, "TIME", (h5_float64_t *)&stash_t[step], 1);
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        rc = H5WriteStepAttribFloat64(H5file_m, "ENERGY", (h5_float64_t *)&stash_meanEnergy[step], 1);
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        rc = H5WriteStepAttribFloat64(H5file_m, "RefPartR", (h5_float64_t *)&stash_RefPartR[step], 3);
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        rc = H5WriteStepAttribFloat64(H5file_m, "RefPartP", (h5_float64_t *)&stash_RefPartP[step], 3);
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        /// Write head/reference particle/tail field information.
        rc = H5WriteStepAttribFloat64(H5file_m, "B-head", (h5_float64_t *)&stash_Bhead[step], 3);
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        rc = H5WriteStepAttribFloat64(H5file_m, "E-head", (h5_float64_t *)&stash_Ehead[step], 3);
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        rc = H5WriteStepAttribFloat64(H5file_m, "B-ref", (h5_float64_t *)&stash_Bref[step], 3);
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        rc = H5WriteStepAttribFloat64(H5file_m, "E-ref", (h5_float64_t *)&stash_Eref[step], 3);
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        rc = H5WriteStepAttribFloat64(H5file_m, "B-tail", (h5_float64_t *)&stash_Btail[step], 3);
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        rc = H5WriteStepAttribFloat64(H5file_m, "E-tail", (h5_float64_t *)&stash_Etail[step], 3);
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

        rc = H5WriteStepAttribFloat64(H5file_m, "spos-head", (h5_float64_t *)&stash_sposHead[step], 1);
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        rc = H5WriteStepAttribFloat64(H5file_m, "spos-ref", (h5_float64_t *)&stash_sposRef[step], 1);
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        rc = H5WriteStepAttribFloat64(H5file_m, "spos-tail", (h5_float64_t *)&stash_sposTail[step], 1);
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

        /// Write number of compute nodes.
        //H5WriteStepAttribFloat64(H5file_m,"nloc", globN, Ippl::getNodes());

        /// Write particle mass and charge per particle. (Consider making these
        /// file attributes.)
        rc = H5WriteStepAttribFloat64(H5file_m, "mpart", &stash_mpart[step], 1);
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        rc = H5WriteStepAttribFloat64(H5file_m, "qi", &stash_qi[step], 1);
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        /// Write normalized emittance.
        rc = H5WriteStepAttribFloat64(H5file_m, "#varepsilon", (h5_float64_t *)&stash_vareps[step], 3);
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        /// Write geometric emittance.
        rc = H5WriteStepAttribFloat64(H5file_m, "#varepsilon-geom", (h5_float64_t *)&stash_geomvareps[step], 3);
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        /// Write beam phase space.
        for(size_t i = 0; i < nLoc; i++)
            farray[i] =  0.0; //beam.getX(i);
        rc = H5PartWriteDataFloat64(H5file_m, "x", farray);
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        for(size_t i = 0; i < nLoc; i++)
            farray[i] =  0.0; //beam.getY(i);
        rc = H5PartWriteDataFloat64(H5file_m, "y", farray);
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        for(size_t i = 0; i < nLoc; i++)
            farray[i] =  0.0; //beam.getZ(i);
        rc = H5PartWriteDataFloat64(H5file_m, "z", farray);
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        for(size_t i = 0; i < nLoc; i++)
            farray[i] =  0.0; //beam.getPx(i);
        rc = H5PartWriteDataFloat64(H5file_m, "px", farray);
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        for(size_t i = 0; i < nLoc; i++)
            farray[i] =  0.0; //beam.getPy(i);
        rc = H5PartWriteDataFloat64(H5file_m, "py", farray);
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        for(size_t i = 0; i < nLoc; i++)
            farray[i] =  0.0; //beam.getPz(i);
        rc = H5PartWriteDataFloat64(H5file_m, "pz", farray);
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        //rc = H5Fflush(H5file_m->file, rc = H5F_SCOPE_GLOBAL);

    }

    /// %Stop timer.
    IpplTimings::stopTimer(H5PartTimer_m);
}

void DataSink::writeStatData(PartBunch &beam, Vector_t FDext[], double sposHead, double sposRef, double sposTail, double E) {
    doWriteStatData(beam, FDext, sposHead, sposRef, sposTail, E, std::vector<std::pair<std::string, unsigned int> >());
}

void DataSink::writeStatData(PartBunch &beam, Vector_t FDext[], double sposHead, double sposRef, double sposTail, const std::vector<std::pair<std::string, unsigned int> >& losses) {
    doWriteStatData(beam, FDext, sposHead, sposRef, sposTail, beam.get_meanEnergy(), losses);
}


void DataSink::doWriteStatData(PartBunch &beam, Vector_t FDext[], double sposHead, double sposRef, double sposTail, double E, const std::vector<std::pair<std::string, unsigned int> > &losses) {

    /// Function steps:

    /// Start timer.
    IpplTimings::startTimer(StatMarkerTimer_m);

    /// Set width of write fields in output files.
    unsigned int pwi = 10;

    /// Calculate beam statistics and gather load balance statistics.
    beam.calcBeamParameters();
    beam.gatherLoadBalanceStatistics();

    double  pathLength = 0.0;
    if (OpalData::getInstance()->isInOPALCyclMode())
        pathLength = beam.getLPath();
    else
        pathLength = sposRef;

    /// Write data to files. If this is the first write to the beam statistics file, write SDDS
    /// header information.
    ofstream os_statData;
    ofstream os_lBalData;
    double Q = beam.getCharge();

    if(Ippl::myNode() == 0) {
        if(firstWriteToStat_m) {
            os_statData.open(statFileName_m.c_str(), ios::out);
            os_statData.precision(15);
            os_statData.setf(ios::scientific, ios::floatfield);
            writeSDDSHeader(os_statData, losses);

            os_lBalData.open(lBalFileName_m.c_str(), ios::out);
            os_lBalData.precision(15);
            os_lBalData.setf(ios::scientific, ios::floatfield);
            os_lBalData << "# " << Ippl::getNodes() << endl;

            firstWriteToStat_m = false;
        } else {
            os_statData.open(statFileName_m.c_str(), ios::app);
            os_statData.precision(15);
            os_statData.setf(ios::scientific, ios::floatfield);

            os_lBalData.open(lBalFileName_m.c_str(), ios::app);
            os_lBalData.precision(15);
            os_lBalData.setf(ios::scientific, ios::floatfield);
        }

        os_statData << beam.getT() << setw(pwi) << "\t"                                       // 1
                    << pathLength << setw(pwi) << "\t"                                        // 2

                    << beam.getTotalNum() << setw(pwi) << "\t"                                // 3
                    << Q << setw(pwi) << "\t"                                                 // 4

                    << E << setw(pwi) << "\t"                                                 // 5

                    << beam.get_rrms()(0) << setw(pwi) << "\t"                                // 6
                    << beam.get_rrms()(1) << setw(pwi) << "\t"                                // 7
                    << beam.get_rrms()(2) << setw(pwi) << "\t"                                // 8

                    << beam.get_prms()(0) << setw(pwi) << "\t"                                // 9
                    << beam.get_prms()(1) << setw(pwi) << "\t"                                // 10
                    << beam.get_prms()(2) << setw(pwi) << "\t"                                // 11

                    << beam.get_norm_emit()(0) << setw(pwi) << "\t"                           // 12
                    << beam.get_norm_emit()(1) << setw(pwi) << "\t"                           // 13
                    << beam.get_norm_emit()(2) << setw(pwi) << "\t"                           // 14

                    << beam.get_rmean()(0)  << setw(pwi) << "\t"                              // 15
                    << beam.get_rmean()(1)  << setw(pwi) << "\t"                              // 16
                    << beam.get_rmean()(2)  << setw(pwi) << "\t"                              // 17

                    << beam.get_maxExtend()(0) << setw(pwi) << "\t"                           // 18
                    << beam.get_maxExtend()(1) << setw(pwi) << "\t"                           // 19
                    << beam.get_maxExtend()(2) << setw(pwi) << "\t"                           // 20

            // Write out Courant Snyder parameters.
                    << beam.get_rprrms()(0) << setw(pwi) << "\t"                              // 21
                    << beam.get_rprrms()(1) << setw(pwi) << "\t"                              // 22
                    << beam.get_rprrms()(2) << setw(pwi) << "\t"                              // 23

                    << 0.0 << setw(pwi) << "\t"                                                // 24
                    << 0.0 << setw(pwi) << "\t"                                                // 25

            // Write out dispersion.
                    << beam.get_Dx() << setw(pwi) << "\t"                                      // 26
                    << beam.get_DDx() << setw(pwi) << "\t"                                     // 27
                    << beam.get_Dy() << setw(pwi) << "\t"                                      // 28
                    << beam.get_DDy() << setw(pwi) << "\t"                                     // 29


            // Write head/reference particle/tail field information.
                    << FDext[0](0) << setw(pwi) << "\t"                                         // 30 B-head x
                    << FDext[0](1) << setw(pwi) << "\t"                                         // 31 B-head y
                    << FDext[0](2) << setw(pwi) << "\t"                                         // 32 B-head z

                    << FDext[1](0) << setw(pwi) << "\t"                                         // 33 E-head x
                    << FDext[1](1) << setw(pwi) << "\t"                                         // 34 E-head y
                    << FDext[1](2) << setw(pwi) << "\t"                                         // 35 E-head z

                    << FDext[2](0) << setw(pwi) << "\t"                                         // 36 B-ref x
                    << FDext[2](1) << setw(pwi) << "\t"                                         // 37 B-ref y
                    << FDext[2](2) << setw(pwi) << "\t"                                         // 38 B-ref z

                    << FDext[3](0) << setw(pwi) << "\t"                                         // 39 E-ref x
                    << FDext[3](1) << setw(pwi) << "\t"                                         // 40 E-ref y
                    << FDext[3](2) << setw(pwi) << "\t"                                         // 41 E-ref z

                    << FDext[4](0) << setw(pwi) << "\t"                                         // 42 B-tail x
                    << FDext[4](1) << setw(pwi) << "\t"                                         // 43 B-tail y
                    << FDext[4](2) << setw(pwi) << "\t"                                         // 44 B-tail z

                    << FDext[5](0) << setw(pwi) << "\t"                                         // 45 E-tail x
                    << FDext[5](1) << setw(pwi) << "\t"                                         // 46 E-tail y
                    << FDext[5](2) << setw(pwi) << "\t"                                         // 47 E-tail z

                    << beam.getdE() << setw(pwi) << "\t";                                       // 48 dE energy spread

        if(Ippl::getNodes() == 1) {
            os_statData << beam.R[0](0) << setw(pwi) << "\t";                                    // 49 R0_x
            os_statData << beam.R[0](1) << setw(pwi) << "\t";                                    // 50 R0_y
            os_statData << beam.R[0](2) << setw(pwi) << "\t";                                    // 51 R0_z
            os_statData << beam.P[0](0) << setw(pwi) << "\t";                                    // 52 P0_x
            os_statData << beam.P[0](1) << setw(pwi) << "\t";                                    // 53 P0_y
            os_statData << beam.P[0](2) << setw(pwi) << "\t";                                    // 54 P0_z
        }

        for(size_t i = 0; i < losses.size(); ++ i) {
            os_statData << losses[i].second << setw(pwi) << "\t";
        }
        os_statData   << endl;

        for(int p = 0; p < Ippl::getNodes(); p++)
            os_lBalData << beam.getLoadBalance(p)  << setw(pwi) << "\t";
        os_lBalData << endl;

        os_statData.close();
        os_lBalData.close();
    }

    /// %Stop timer.
    IpplTimings::stopTimer(StatMarkerTimer_m);
}

void DataSink::writeStatData(EnvelopeBunch &beam, Vector_t FDext[], double sposHead, double sposRef, double sposTail) {
    /// Function steps:

    /// Start timer.
    IpplTimings::startTimer(StatMarkerTimer_m);

    /// Set width of write fields in output files.
    unsigned int pwi = 10;

    /// Calculate beam statistics and gather load balance statistics.
    beam.calcBeamParameters();
    beam.gatherLoadBalanceStatistics();

    /// Write data to files. If this is the first write to the beam statistics file, write SDDS
    /// header information.
    ofstream os_statData;
    ofstream os_lBalData;
    double en = beam.get_meanEnergy() * 1e-6;

    if(Ippl::myNode() == 0) {
        if(firstWriteToStat_m) {
            os_statData.open(statFileName_m.c_str(), ios::out);
            os_statData.precision(15);
            os_statData.setf(ios::scientific, ios::floatfield);
            writeSDDSHeader(os_statData);

            os_lBalData.open(lBalFileName_m.c_str(), ios::out);
            os_lBalData.precision(15);
            os_lBalData.setf(ios::scientific, ios::floatfield);
            os_lBalData << "# " << Ippl::getNodes() << endl;

            firstWriteToStat_m = false;
        } else {
            os_statData.open(statFileName_m.c_str(), ios::app);
            os_statData.precision(15);
            os_statData.setf(ios::scientific, ios::floatfield);

            os_lBalData.open(lBalFileName_m.c_str(), ios::app);
            os_lBalData.precision(15);
            os_lBalData.setf(ios::scientific, ios::floatfield);
        }


        os_statData << beam.getT() << setw(pwi) << "\t"                                       // 1
                    << sposRef << setw(pwi) << "\t"                                           // 2

                    << beam.getTotalNum() << setw(pwi) << "\t"                                // 3
                    << beam.getTotalNum() * beam.getChargePerParticle() << setw(pwi) << "\t"  // 4
                    << en << setw(pwi) << "\t"                                                // 5

                    << beam.get_rrms()(0) << setw(pwi) << "\t"                                // 6
                    << beam.get_rrms()(1) << setw(pwi) << "\t"                                // 7
                    << beam.get_rrms()(2) << setw(pwi) << "\t"                                // 8

                    << beam.get_prms()(0) << setw(pwi) << "\t"                                // 9
                    << beam.get_prms()(1) << setw(pwi) << "\t"                                // 10
                    << beam.get_prms()(2) << setw(pwi) << "\t"                                // 11

                    << beam.get_norm_emit()(0) << setw(pwi) << "\t"                           // 12
                    << beam.get_norm_emit()(1) << setw(pwi) << "\t"                           // 13
                    << beam.get_norm_emit()(2) << setw(pwi) << "\t"                           // 14

                    << beam.get_rmean()(0)  << setw(pwi) << "\t"                              // 15
                    << beam.get_rmean()(1)  << setw(pwi) << "\t"                              // 16
                    << beam.get_rmean()(2)  << setw(pwi) << "\t"                              // 17

                    << beam.get_maxExtend()(0) << setw(pwi) << "\t"                           // 18
                    << beam.get_maxExtend()(1) << setw(pwi) << "\t"                           // 19
                    << beam.get_maxExtend()(2) << setw(pwi) << "\t"                           // 20

            // Write out Courant Snyder parameters.
                    << 0.0  << setw(pwi) << "\t"                                              // 21
                    << 0.0  << setw(pwi) << "\t"                                              // 22

                    << 0.0 << setw(pwi) << "\t"                                               // 23
                    << 0.0 << setw(pwi) << "\t"                                               // 24

            // Write out dispersion.
                    << beam.get_Dx() << setw(pwi) << "\t"                                     // 25
                    << beam.get_DDx() << setw(pwi) << "\t"                                    // 26
                    << beam.get_Dy() << setw(pwi) << "\t"                                     // 27
                    << beam.get_DDy() << setw(pwi) << "\t"                                    // 28


            // Write head/reference particle/tail field information.
                    << FDext[0](0) << setw(pwi) << "\t"                                       // 29 B-head x
                    << FDext[0](1) << setw(pwi) << "\t"                                       // 30 B-head y
                    << FDext[0](2) << setw(pwi) << "\t"                                       // 31 B-head z

                    << FDext[1](0) << setw(pwi) << "\t"                                       // 32 E-head x
                    << FDext[1](1) << setw(pwi) << "\t"                                       // 33 E-head y
                    << FDext[1](2) << setw(pwi) << "\t"                                       // 34 E-head z

                    << FDext[2](0) << setw(pwi) << "\t"                                       // 35 B-ref x
                    << FDext[2](1) << setw(pwi) << "\t"                                       // 36 B-ref y
                    << FDext[2](2) << setw(pwi) << "\t"                                       // 37 B-ref z

                    << FDext[3](0) << setw(pwi) << "\t"                                       // 38 E-ref x
                    << FDext[3](1) << setw(pwi) << "\t"                                       // 39 E-ref y
                    << FDext[3](2) << setw(pwi) << "\t"                                       // 40 E-ref z

                    << FDext[4](0) << setw(pwi) << "\t"                                       // 41 B-tail x
                    << FDext[4](1) << setw(pwi) << "\t"                                       // 42 B-tail y
                    << FDext[4](2) << setw(pwi) << "\t"                                       // 43 B-tail z

                    << FDext[5](0) << setw(pwi) << "\t"                                       // 44 E-tail x
                    << FDext[5](1) << setw(pwi) << "\t"                                       // 45 E-tail y
                    << FDext[5](2) << setw(pwi) << "\t"                                       // 46 E-tail z

                    << beam.get_dEdt() << setw(pwi) << "\t"                                   // 47 dE energy spread
                    << endl;


        for(int p = 0; p < Ippl::getNodes(); p++)
            os_lBalData << beam.getLoadBalance(p)  << setw(pwi) << "\t";
        os_lBalData << endl;

        os_statData.close();
        os_lBalData.close();
    }

    /// %Stop timer.
    IpplTimings::stopTimer(StatMarkerTimer_m);
}

void DataSink::writeSDDSHeader(ofstream &outputFile) {
    writeSDDSHeader(outputFile,
                    std::vector<std::pair<std::string, unsigned int> >());
}
void DataSink::writeSDDSHeader(ofstream &outputFile,
                               const std::vector<std::pair<std::string, unsigned int> > &losses) {
    OPALTimer::Timer simtimer;

    string dateStr(simtimer.date());
    string timeStr(simtimer.time());

    outputFile << "SDDS1" << endl;
    outputFile << "&description text=\"Statistics data " << OpalData::getInstance()->getInputFn() << " " << dateStr << " " << timeStr << "\" " << endl;
    outputFile << ", contents=\"stat parameters\" &end" << endl;

    outputFile << "&parameter name=processors, type=long, ";
    outputFile << "description=\"Number of Cores used\" &end" << endl;

    outputFile << "&parameter name=revision, type=string, "
               << "description=\"svn revision of opal\" &end\n";

    outputFile << "&column name=t, type=double, units=s, ";
    outputFile << "description=\"1 Time\" &end" << endl;
    outputFile << "&column name=s, type=double, units=m, ";
    outputFile << "description=\"2 Average Longitudinal Position\" &end" << endl;

    outputFile << "&column name=numParticles, type=long, units=1, ";
    outputFile << "description=\"3 Number of Macro Particles\" &end" << endl;
    outputFile << "&column name=charge, type=double, units=1, ";
    outputFile << "description=\"4 Bunch Charge\" &end" << endl;

    outputFile << "&column name=energy, type=double, units=MeV, ";
    outputFile << "description=\"5 Mean Energy\" &end" << endl;

    outputFile << "&column name=rms_x, type=double, units=m , ";
    outputFile << "description=\"6 RMS Beamsize in x  \" &end" << endl;
    outputFile << "&column name=rms_y, type=double, units=m , ";
    outputFile << "description=\"7 RMS Beamsize in y  \" &end" << endl;
    outputFile << "&column name=rms_s, type=double, units=m , ";
    outputFile << "description=\"8 RMS Beamsize in s  \" &end" << endl;

    outputFile << "&column name=rms_px, type=double, units=1 , ";
    outputFile << "description=\"9 RMS Momenta in x  \" &end" << endl;
    outputFile << "&column name=rms_py, type=double, units=1 , ";
    outputFile << "description=\"10 RMS Momenta in y  \" &end" << endl;
    outputFile << "&column name=rms_ps, type=double, units=1 , ";
    outputFile << "description=\"11 RMS Momenta in s  \" &end" << endl;

    outputFile << "&column name=emit_x, type=double, units=m , ";
    outputFile << "description=\"12 Normalized Emittance x  \" &end" << endl;
    outputFile << "&column name=emit_y, type=double, units=m , ";
    outputFile << "description=\"13 Normalized Emittance y  \" &end" << endl;
    outputFile << "&column name=emit_s, type=double, units=m , ";
    outputFile << "description=\"14 Normalized Emittance s  \" &end" << endl;

    outputFile << "&column name=mean_x, type=double, units=m , ";
    outputFile << "description=\"15 Mean Beam Position in x  \" &end" << endl;
    outputFile << "&column name=mean_y, type=double, units=m , ";
    outputFile << "description=\"16 Mean Beam Position in y  \" &end" << endl;
    outputFile << "&column name=mean_s, type=double, units=m , ";
    outputFile << "description=\"17 Mean Beam Position in s  \" &end" << endl;

    outputFile << "&column name=max_x, type=double, units=m , ";
    outputFile << "description=\"18 Max Beamsize in x  \" &end" << endl;
    outputFile << "&column name=max_y, type=double, units=m , ";
    outputFile << "description=\"19 Max Beamsize in y  \" &end" << endl;
    outputFile << "&column name=max_s, type=double, units=m , ";
    outputFile << "description=\"20 Max Beamsize in s  \" &end" << endl;

    outputFile << "&column name=xpx, type=double, units=1 , ";
    outputFile << "description=\"21 Correlation xpx  \" &end" << endl;
    outputFile << "&column name=ypy, type=double, units=1 , ";
    outputFile << "description=\"22 Correlation ypyy  \" &end" << endl;

    outputFile << "&column name=zpz, type=double, units=1 , ";
    outputFile << "description=\"23 Correlation zpz  \" &end" << endl;

    outputFile << "&column name=notused1, type=double, units=1 , ";
    outputFile << "description=\"24 notused1 in y  \" &end" << endl;

    outputFile << "&column name=notused2, type=double, units=1 , ";
    outputFile << "description=\"25 notused2 in y  \" &end" << endl;

    outputFile << "&column name=Dx, type=double, units=m , ";
    outputFile << "description=\"26 Dispersion in x  \" &end" << endl;
    outputFile << "&column name=DDx, type=double, units=1 , ";
    outputFile << "description=\"27 Derivative of dispersion in x  \" &end" << endl;

    outputFile << "&column name=Dy, type=double, units=m , ";
    outputFile << "description=\"28 Dispersion in y  \" &end" << endl;
    outputFile << "&column name=DDy, type=double, units=1 , ";
    outputFile << "description=\"29 Derivative of dispersion in y  \" &end" << endl;

    outputFile << "&column name=Bx_head, type=double, units=T , ";
    outputFile << "description=\"30 Bx-Field component of head particle  \" &end" << endl;
    outputFile << "&column name=By_head, type=double, units=T , ";
    outputFile << "description=\"31 By-Field component of head particle  \" &end" << endl;
    outputFile << "&column name=Bz_head, type=double, units=T , ";
    outputFile << "description=\"32 Bz-Field component of head particle  \" &end" << endl;

    outputFile << "&column name=Ex_head, type=double, units=MV/m , ";
    outputFile << "description=\"33 Ex-Field component of head particle  \" &end" << endl;
    outputFile << "&column name=Ey_head, type=double, units=MV/m , ";
    outputFile << "description=\"34 Ey-Field component of head particle  \" &end" << endl;
    outputFile << "&column name=Ez_head, type=double, units=MV/m , ";
    outputFile << "description=\"35 Ez-Field component of head particle  \" &end" << endl;

    outputFile << "&column name=Bx_ref, type=double, units=T , ";
    outputFile << "description=\"36 Bx-Field component of ref particle  \" &end" << endl;
    outputFile << "&column name=By_ref, type=double, units=T , ";
    outputFile << "description=\"37 By-Field component of ref particle  \" &end" << endl;
    outputFile << "&column name=Bz_ref, type=double, units=T , ";
    outputFile << "description=\"38 Bz-Field component of ref particle  \" &end" << endl;

    outputFile << "&column name=Ex_ref, type=double, units=MV/m , ";
    outputFile << "description=\"39 Ex-Field component of ref particle  \" &end" << endl;
    outputFile << "&column name=Ey_ref, type=double, units=MV/m , ";
    outputFile << "description=\"40 Ey-Field component of ref particle  \" &end" << endl;
    outputFile << "&column name=Ez_ref, type=double, units=MV/m , ";
    outputFile << "description=\"41 Ez-Field component of ref particle  \" &end" << endl;

    outputFile << "&column name=Bx_tail, type=double, units=T , ";
    outputFile << "description=\"42 Bx-Field component of tail particle  \" &end" << endl;
    outputFile << "&column name=By_tail, type=double, units=T , ";
    outputFile << "description=\"43 By-Field component of tail particle  \" &end" << endl;
    outputFile << "&column name=Bz_tail, type=double, units=T , ";
    outputFile << "description=\"44 Bz-Field component of tail particle  \" &end" << endl;

    outputFile << "&column name=Ex_tail, type=double, units=MV/m , ";
    outputFile << "description=\"45 Ex-Field component of tail particle  \" &end" << endl;
    outputFile << "&column name=Ey_tail, type=double, units=MV/m , ";
    outputFile << "description=\"46 Ey-Field component of tail particle  \" &end" << endl;
    outputFile << "&column name=Ez_tail, type=double, units=MV/m , ";
    outputFile << "description=\"47 Ez-Field component of tail particle  \" &end" << endl;

    outputFile << "&column name=dE, type=double, units=MeV , ";
    outputFile << "description=\"48 energy spread of the beam  \" &end" << endl;

    unsigned int columnStart = 49;
    if(Ippl::getNodes() == 1) {
        outputFile << "&column name=R0_x, type=double, units=m , ";
        outputFile << "description=\"49 R0 Particle Position in x  \" &end" << endl;
        outputFile << "&column name=R0_y, type=double, units=m , ";
        outputFile << "description=\"50 R0 Particle Position in y  \" &end" << endl;
        outputFile << "&column name=R0_s, type=double, units=m , ";
        outputFile << "description=\"51 R0 Particle Position in s  \" &end" << endl;
        outputFile << "&column name=P0_x, type=double, units=1 , ";
        outputFile << "description=\"52 R0 Particle Position in x  \" &end" << endl;
        outputFile << "&column name=P0_y, type=double, units=1 , ";
        outputFile << "description=\"53 R0 Particle Position in y  \" &end" << endl;
        outputFile << "&column name=P0_s, type=double, units=1 , ";
        outputFile << "description=\"54 R0 Particle Position in s  \" &end" << endl;
        columnStart = 55;
    }

    for (size_t i = 0; i < losses.size(); ++ i) {
        outputFile << "&column name=" << losses[i].first << ", type=long, units=1, ";
        outputFile << "description=\"" << columnStart++ << " " << losses[i].second << "\" &end" << endl;
    }
    outputFile << "&data mode=ascii, no_row_counts=1 &end" << endl;

    outputFile << Ippl::getNodes() << endl;
    outputFile << PACKAGE_NAME << " " << PACKAGE_VERSION << " svn rev. " << GIT_VERSION << endl;
}


void DataSink::writePartlossZASCII(PartBunch &beam, BoundaryGeometry &bg, string fn) {

    size_t temp = lossWrCounter_m ;

    string ffn = fn + convert2Int(temp) + string("Z.dat");
    std::unique_ptr<Inform> ofp(new Inform(NULL, ffn.c_str(), Inform::OVERWRITE, 0));
    Inform &fid = *ofp;
    setInform(fid);
    fid.precision(6);

    string ftrn =  fn + string("triangle") + convert2Int(temp) + string(".dat");
    std::unique_ptr<Inform> oftr(new Inform(NULL, ftrn.c_str(), Inform::OVERWRITE, 0));
    Inform &fidtr = *oftr;
    setInform(fidtr);
    fidtr.precision(6);

    Vector_t Geo_nr = bg.getnr();
    Vector_t Geo_hr = bg.gethr();
    Vector_t Geo_mincoords = bg.getmincoords();
    double t = beam.getT();
    double t_step = t * 1.0e9;
    double* prPartLossZ = new double[bg.getnr() (2)];
    double* sePartLossZ = new double[bg.getnr() (2)];
    double* fePartLossZ = new double[bg.getnr() (2)];
    fidtr << "# Time/ns" << std::setw(18) << "Triangle_ID" << std::setw(18)
          << "Xcoordinates (m)" << std::setw(18)
          << "Ycoordinates (m)" << std::setw(18)
          << "Zcoordinates (m)" << std::setw(18)
          << "Primary part. charge (C)" << std::setw(40)
          << "Field emit. part. charge (C)" << std::setw(40)
          << "Secondary emit. part. charge (C)" << std::setw(40) << endl ;
    for(int i = 0; i < Geo_nr(2) ; i++) {
        prPartLossZ[i] = 0;
        sePartLossZ[i] = 0;
        fePartLossZ[i] = 0;
        for(int j = 0; j < bg.getNumBFaces(); j++) {
            if(((Geo_mincoords[2] + Geo_hr(2)*i) < bg.TriBarycenters_m[j](2))
               && (bg.TriBarycenters_m[j](2) < (Geo_hr(2)*i + Geo_hr(2) + Geo_mincoords[2]))) {
                prPartLossZ[i] += bg.TriPrPartloss_m[j];
                sePartLossZ[i] += bg.TriSePartloss_m[j];
                fePartLossZ[i] += bg.TriFEPartloss_m[j];
            }

        }
    }
    for(int j = 0; j < bg.getNumBFaces(); j++) {
        fidtr << t_step << std::setw(18) << j << std::setw(18)// fixme: maybe gether particle loss data, i.e., do a reduce() for each triangle in each node befor write to file.
              << bg.TriBarycenters_m[j](0) << std::setw(18)
              << bg.TriBarycenters_m[j](1) << std::setw(18)
              << bg.TriBarycenters_m[j](2) <<  std::setw(40)
              << -bg.TriPrPartloss_m[j] << std::setw(40)
              << -bg.TriFEPartloss_m[j] <<  std::setw(40)
              << -bg.TriSePartloss_m[j] << endl;
    }
    fid << "# Delta_Z/m" << std::setw(18)
        << "Zcoordinates (m)" << std::setw(18)
        << "Primary part. charge (C)" << std::setw(40)
        << "Field emit. part. charge (C)" << std::setw(40)
        << "Secondary emit. part. charge (C)" << std::setw(40) << "t" << endl ;


    for(int i = 0; i < Geo_nr(2) ; i++) {
        double primaryPLoss = -prPartLossZ[i];
        double secondaryPLoss = -sePartLossZ[i];
        double fieldemissionPLoss = -fePartLossZ[i];
        reduce(primaryPLoss, primaryPLoss, OpAddAssign());
        reduce(secondaryPLoss, secondaryPLoss, OpAddAssign());
        reduce(fieldemissionPLoss, fieldemissionPLoss, OpAddAssign());
        fid << Geo_hr(2) << std::setw(18)
            << Geo_mincoords[2] + Geo_hr(2)*i << std::setw(18)
            << primaryPLoss << std::setw(40)
            << fieldemissionPLoss << std::setw(40)
            << secondaryPLoss << std::setw(40) << t << endl;
    }
    lossWrCounter_m++;
}

void DataSink::writeSurfaceInteraction(PartBunch &beam, long long &step, BoundaryGeometry &bg, string fn) {

    if (!doHDF5_m) return;

    h5_int64_t rc;
    /// Start timer.
    IpplTimings::startTimer(H5PartTimer_m);
    if(firstWriteH5Surface_m) {
        firstWriteH5Surface_m = false;
#ifdef PARALLEL_IO
        H5fileS_m = H5OpenFile(surfaceLossFileName_m.c_str(), H5_FLUSH_STEP | H5_O_WRONLY, Ippl::getComm());
#else
        H5fileS_m = H5OpenFile(surfaceLossFileName_m.c_str(), H5_FLUSH_STEP | H5_O_WRONLY, 0);
#endif

        if(H5fileS_m == (void*)H5_ERR) {
            throw OpalException("DataSink::writeSurfaceInteraction",
                                "failed to open h5 file '" + surfaceLossFileName_m + "' for surface loss");
        }


    }
    int nTot = bg.getNumBFaces();

    int N_mean = static_cast<int>(floor(nTot / Ippl::getNodes()));
    int N_extra = static_cast<int>(nTot - N_mean * Ippl::getNodes());
    int pc = 0;
    int count = 0;
    if(Ippl::myNode() == 0) {
        N_mean += N_extra;
    }
    std::unique_ptr<char[]> varray(new char[(N_mean)*sizeof(double)]);
    double *farray = reinterpret_cast<double *>(varray.get());
    h5_int64_t *larray = reinterpret_cast<h5_int64_t *>(varray.get());


    rc = H5SetStep(H5fileS_m, step);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5PartSetNumParticles(H5fileS_m, N_mean);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    double    qi = beam.getChargePerParticle();
    rc = H5WriteStepAttribFloat64(H5fileS_m, "qi", &qi, 1);

    std::unique_ptr<double[]> tmploss(new double[nTot]);
    for(int i = 0; i < nTot; i++)
        tmploss[i] = 0.0;
    //memset( tmploss, 0.0, nTot * sizeof(double));
    reduce(bg.TriPrPartloss_m, bg.TriPrPartloss_m + nTot, tmploss.get(), OpAddAssign()); // may be removed if we have parallelized the geometry .

    for(int i = 0; i < nTot; i++) {
        if(pc == Ippl::myNode()) {
            if(count < N_mean) {
                if(pc != 0) {
                    //farray[count] =  bg.TriPrPartloss_m[Ippl::myNode()*N_mean+count];
                    if(((bg.TriBGphysicstag_m[pc * N_mean + count + N_extra] & (BGphysics::Absorption)) == (BGphysics::Absorption)) && (((bg.TriBGphysicstag_m[pc * N_mean + count + N_extra] & (BGphysics::FNEmission)) != (BGphysics::FNEmission)) && ((bg.TriBGphysicstag_m[pc * N_mean + count + N_extra] & (BGphysics::SecondaryEmission)) != (BGphysics::SecondaryEmission)))) {
                        farray[count] = 0.0;
                    } else {
                        farray[count] =  tmploss[pc * N_mean + count + N_extra];
                    }
                    count ++;
                } else {
                    if(((bg.TriBGphysicstag_m[count] & (BGphysics::Absorption)) == (BGphysics::Absorption)) && (((bg.TriBGphysicstag_m[count] & (BGphysics::FNEmission)) != (BGphysics::FNEmission)) && ((bg.TriBGphysicstag_m[count] & (BGphysics::SecondaryEmission)) != (BGphysics::SecondaryEmission)))) {
                        farray[count] = 0.0;
                    } else {
                        farray[count] =  tmploss[count];
                    }
                    count ++;

                }
            }
        }
        pc++;
        if(pc == Ippl::getNodes())
            pc = 0;
    }
    rc = H5PartWriteDataFloat64(H5fileS_m, "PrimaryLoss", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    for(int i = 0; i < nTot; i++)
        tmploss[i] = 0.0;
    reduce(bg.TriSePartloss_m, bg.TriSePartloss_m + nTot, tmploss.get(), OpAddAssign()); // may be removed if we parallelize the geometry as well.
    count = 0;
    pc = 0;
    for(int i = 0; i < nTot; i++) {
        if(pc == Ippl::myNode()) {
            if(count < N_mean) {

                if(pc != 0) {

                    if(((bg.TriBGphysicstag_m[pc * N_mean + count + N_extra] & (BGphysics::Absorption)) == (BGphysics::Absorption)) && (((bg.TriBGphysicstag_m[pc * N_mean + count + N_extra] & (BGphysics::FNEmission)) != (BGphysics::FNEmission)) && ((bg.TriBGphysicstag_m[pc * N_mean + count + N_extra] & (BGphysics::SecondaryEmission)) != (BGphysics::SecondaryEmission)))) {
                        farray[count] = 0.0;
                    } else {
                        farray[count] =  tmploss[pc * N_mean + count + N_extra];
                    }
                    count ++;
                } else {
                    if(((bg.TriBGphysicstag_m[count] & (BGphysics::Absorption)) == (BGphysics::Absorption)) && (((bg.TriBGphysicstag_m[count] & (BGphysics::FNEmission)) != (BGphysics::FNEmission)) && ((bg.TriBGphysicstag_m[count] & (BGphysics::SecondaryEmission)) != (BGphysics::SecondaryEmission)))) {
                        farray[count] = 0.0;
                    } else {
                        farray[count] =  tmploss[count];
                    }
                    count ++;

                }
            }
        }
        pc++;
        if(pc == Ippl::getNodes())
            pc = 0;
    }
    rc = H5PartWriteDataFloat64(H5fileS_m, "SecondaryLoss", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    for(int i = 0; i < nTot; i++)
        tmploss[i] = 0.0;

    reduce(bg.TriFEPartloss_m, bg.TriFEPartloss_m + nTot, tmploss.get(), OpAddAssign()); // may be removed if we parallelize the geometry as well.
    count = 0;
    pc = 0;
    for(int i = 0; i < nTot; i++) {
        if(pc == Ippl::myNode()) {
            if(count < N_mean) {

                if(pc != 0) {

                    if(((bg.TriBGphysicstag_m[pc * N_mean + count + N_extra] & (BGphysics::Absorption)) == (BGphysics::Absorption)) && (((bg.TriBGphysicstag_m[pc * N_mean + count + N_extra] & (BGphysics::FNEmission)) != (BGphysics::FNEmission)) && ((bg.TriBGphysicstag_m[pc * N_mean + count + N_extra] & (BGphysics::SecondaryEmission)) != (BGphysics::SecondaryEmission)))) {
                        farray[count] = 0.0;
                    } else {
                        farray[count] =  tmploss[pc * N_mean + count + N_extra];
                    }
                    count ++;
                } else {
                    if(((bg.TriBGphysicstag_m[count] & (BGphysics::Absorption)) == (BGphysics::Absorption)) && (((bg.TriBGphysicstag_m[count] & (BGphysics::FNEmission)) != (BGphysics::FNEmission)) && ((bg.TriBGphysicstag_m[count] & (BGphysics::SecondaryEmission)) != (BGphysics::SecondaryEmission)))) {
                        farray[count] = 0.0;
                    } else {
                        farray[count] =  tmploss[count];
                    }
                    count ++;

                }

            }
        }
        pc++;
        if(pc == Ippl::getNodes())
            pc = 0;
    }
    rc = H5PartWriteDataFloat64(H5fileS_m, "FNEmissionLoss", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    count = 0;
    pc = 0;
    for(int i = 0; i < nTot; i++) {
        if(pc == Ippl::myNode()) {
            if(count < N_mean) {
                if(pc != 0)
                    larray[count] =  Ippl::myNode() * N_mean + count + N_extra; // node 0 will be 0*N_mean+count+N_extra also correct.
                else
                    larray[count] = count;
                count ++;
            }
        }
        pc++;
        if(pc == Ippl::getNodes())
            pc = 0;
    }
    rc = H5PartWriteDataInt64(H5fileS_m, "TriangleID", larray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    /// %Stop timer.
    IpplTimings::stopTimer(H5PartTimer_m);

}

void DataSink::writeImpactStatistics(PartBunch &beam, long long &step, size_t &impact, double &sey_num,
                                     size_t numberOfFieldEmittedParticles, bool nEmissionMode, string fn) {

    double charge = 0.0;
    size_t Npart = 0;
    double Npart_d = 0.0;
    if(!nEmissionMode) {
        charge = -1.0 * beam.getCharge();
        //reduce(charge, charge, OpAddAssign());
        Npart_d = -1.0 * charge / beam.getChargePerParticle();
    } else {
        Npart = beam.getTotalNum();
    }
    if(Ippl::myNode() == 0) {
        string ffn = fn + string(".dat");

        std::unique_ptr<Inform> ofp(new Inform(NULL, ffn.c_str(), Inform::APPEND, 0));
        Inform &fid = *ofp;
        setInform(fid);

        fid.precision(6);
        fid << setiosflags(ios::scientific);
        double t = beam.getT() * 1.0e9;
        if(!nEmissionMode) {

            if(step == 0) {
                fid << "#Time/ns"  << std::setw(18) << "#Geometry impacts" << std::setw(18) << "tot_sey" << std::setw(18)
                    << "TotalCharge" << std::setw(18) << "PartNum" << " numberOfFieldEmittedParticles " << endl ;
            }
            fid << t << std::setw(18) << impact << std::setw(18) << sey_num << std::setw(18) << charge
                << std::setw(18) << Npart_d << std::setw(18) << numberOfFieldEmittedParticles << endl ;
        } else {

            if(step == 0) {
                fid << "#Time/ns"  << std::setw(18) << "#Geometry impacts" << std::setw(18) << "tot_sey" << std::setw(18)
                    << "ParticleNumber" << " numberOfFieldEmittedParticles " << endl;
            }
            fid << t << std::setw(18) << impact << std::setw(18) << sey_num
                << std::setw(18) << double(Npart) << std::setw(18) << numberOfFieldEmittedParticles << endl ;
        }
    }
}

void DataSink::writeGeomToVtk(BoundaryGeometry &bg, string fn) {
    if(Ippl::myNode() == 0) {
        bg.writeGeomToVtk (fn);
    }
}


void DataSink::storeOneBunch(const PartBunch &beam, const string fn_appendix) {

    if (!doHDF5_m) return;

    h5_int64_t rc;
    /// Define file names.
    string fn = OpalData::getInstance()->getInputBasename();
    fn = fn + fn_appendix + string(".h5");
    h5_file_t *H5file;

    const size_t nLoc = beam.getLocalNum();

    std::unique_ptr<char[]> varray(new char[(nLoc)*sizeof(double)]);
    double *farray = reinterpret_cast<double *>(varray.get());
    h5_int64_t *larray = reinterpret_cast<h5_int64_t *>(varray.get());


#ifdef PARALLEL_IO
    H5file = H5OpenFile(fn.c_str(), H5_FLUSH_STEP | H5_O_WRONLY, Ippl::getComm());
#else
    H5file = H5OpenFile(fn.c_str(), H5_FLUSH_STEP | H5_O_WRONLY, 0);
#endif

    if(H5file == (void*)H5_ERR) {
        throw OpalException("DataSink::storeOneBunch",
                            "failed to open h5 file '" + fn + "' for backup bunch in OPAL-CYCL");
    }

    h5_int64_t H5call = 0;
    /// Set current record/time step.
    rc = H5SetStep(H5file, H5call);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5PartSetNumParticles(H5file, nLoc);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    stringstream OPAL_version;
    OPAL_version << PACKAGE_NAME << " " << PACKAGE_VERSION << " svn rev. " << GIT_VERSION;
    rc = H5WriteFileAttribString(H5file, "OPAL_version", OPAL_version.str().c_str());
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteFileAttribString(H5file, "StepUnit", " ");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteFileAttribString(H5file, "xUnit", "m");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteFileAttribString(H5file, "yUnit", "m");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteFileAttribString(H5file, "zUnit", "m");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteFileAttribString(H5file, "pxUnit", "#beta#gamma");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteFileAttribString(H5file, "pyUnit", "#beta#gamma");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteFileAttribString(H5file, "pzUnit", "#beta#gamma");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteFileAttribString(H5file, "massUnit", "GeV");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteFileAttribString(H5file, "chargeUnit", "C");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteFileAttribString(H5file, "ptype", "1");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteStepAttribInt64(H5file, "Step", &H5call, 1);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    for(size_t i = 0; i < nLoc; i++)
        farray[i] =  beam.R[i](0) / 1000.0;
    rc = H5PartWriteDataFloat64(H5file, "x", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    for(size_t i = 0; i < nLoc; i++)
        farray[i] =  beam.R[i](1) / 1000.0;
    rc = H5PartWriteDataFloat64(H5file, "y", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    for(size_t i = 0; i < nLoc; i++)
        farray[i] =  beam.R[i](2) / 1000.0;
    rc = H5PartWriteDataFloat64(H5file, "z", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    for(size_t i = 0; i < nLoc; i++)
        farray[i] =  beam.P[i](0);
    rc = H5PartWriteDataFloat64(H5file, "px", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    for(size_t i = 0; i < nLoc; i++)
        farray[i] =  beam.P[i](1);
    rc = H5PartWriteDataFloat64(H5file, "py", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    for(size_t i = 0; i < nLoc; i++)
        farray[i] =  beam.P[i](2);
    rc = H5PartWriteDataFloat64(H5file, "pz", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    for(size_t i = 0; i < nLoc; i++)
        farray[i] =  beam.M[i];
    rc = H5PartWriteDataFloat64(H5file, "mass", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    for(size_t i = 0; i < nLoc; i++)
        farray[i] =  beam.Q[i];
    rc = H5PartWriteDataFloat64(H5file, "charge", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    for(size_t i = 0; i < nLoc; i++)
        larray[i] =  beam.PType[i];
    rc = H5PartWriteDataInt64(H5file, "ptype", larray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5CloseFile(H5file);
    H5file = NULL;
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    Ippl::Comm->barrier();
}

bool DataSink::readOneBunch(PartBunch &beam, const string fn_appendix, const size_t BinID) {

    if (!doHDF5_m) return false;

    h5_int64_t rc;
    /// Define file names.
    string fn = OpalData::getInstance()->getInputBasename();
    fn = fn + fn_appendix + string(".h5");
    h5_file_t *H5file;

#ifdef PARALLEL_IO
    H5file = H5OpenFile(fn.c_str(), H5_O_RDONLY, Ippl::getComm());
#else
    H5file = H5OpenFile(fn.c_str(), H5_O_RDONLY);
#endif
    if(H5file == (void*)H5_ERR) {
        throw OpalException("DataSink::readOneBunch",
                            "failed to open h5 file '" + fn + "'");
    }

    h5_int64_t H5call = 0;
    rc = H5SetStep(H5file, H5call);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    const size_t globalN = (size_t)H5PartGetNumParticles(H5file);

    const size_t myNode = (size_t)Ippl::myNode();
    size_t numberOfParticlesPerNode = (size_t) floor((double) globalN / Ippl::getNodes());
    long long starti = myNode * numberOfParticlesPerNode;
    long long endi = 0;
    // ensure that we dont miss any particle in the end
    if(myNode + 1 == (size_t) Ippl::getNodes())
        endi = -1;
    else
        endi = starti + numberOfParticlesPerNode;

    rc = H5PartSetView(H5file, starti, endi);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    const size_t InjectN = (size_t)H5PartGetNumParticles(H5file);


    std::unique_ptr<char[]> varray(new char[(InjectN*sizeof(double))]);
    double *farray = reinterpret_cast<double *>(varray.get());

    const size_t LocalNum = beam.getLocalNum();
    const size_t NewLocalNum = LocalNum + InjectN;

    beam.create(InjectN);

    rc = H5PartReadDataFloat64(H5file, "x", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    for(size_t ii = LocalNum; ii < NewLocalNum; ii++) {
        beam.R[ii](0) = farray[ii - LocalNum] * 1000.0; // m --> mm
        // unlike the process in of restart, here we set bin index forcely,
        // because this new bunch should  be a new bin with lowest energy.
        beam.Bin[ii] = BinID;
    }
    rc = H5PartReadDataFloat64(H5file, "y", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    for(size_t ii = LocalNum; ii < NewLocalNum; ii++)
        beam.R[ii](1) = farray[ii - LocalNum] * 1000.0; // m --> mm;

    rc = H5PartReadDataFloat64(H5file, "z", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    for(size_t ii = LocalNum; ii < NewLocalNum; ii++)
        beam.R[ii](2) = farray[ii - LocalNum] * 1000.0; // m --> mm;

    rc = H5PartReadDataFloat64(H5file, "px", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    for(size_t ii = LocalNum; ii < NewLocalNum; ii++)
        beam.P[ii](0) = farray[ii - LocalNum];

    rc = H5PartReadDataFloat64(H5file, "py", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    for(size_t ii = LocalNum; ii < NewLocalNum; ii++)
        beam.P[ii](1) = farray[ii - LocalNum];

    rc = H5PartReadDataFloat64(H5file, "pz", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    for(size_t ii = LocalNum; ii < NewLocalNum; ii++)
        beam.P[ii](2) = farray[ii - LocalNum];

    rc = H5PartReadDataFloat64(H5file, "mass", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    for(size_t ii = LocalNum; ii < NewLocalNum; ii++)
        beam.M[ii] = farray[ii - LocalNum];

    rc = H5PartReadDataFloat64(H5file, "charge", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    for(size_t ii = LocalNum; ii < NewLocalNum; ii++)
        beam.Q[ii] = farray[ii - LocalNum];

    for(size_t ii = LocalNum; ii < NewLocalNum; ii++)
        beam.PType[ii] = 0;

    // update the bin status
    if(beam.weHaveBins())
        beam.pbin_m->updateStatus(BinID + 1, InjectN);

    Ippl::Comm->barrier();
    rc = H5CloseFile(H5file);
    H5file = NULL;
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    // update statistics parameters of PartBunch
    // allocate ID for new particles
    // do not call boundp() here

    return true;

}

/** \brief Find out which if we write HDF5 or not
 *
 *
 *
 */
bool DataSink::doHDF5() {
    return doHDF5_m;
}

/** \brief Find out which flavor has written the data of
 *   the h5 file.
 *
 *
 */
bool DataSink::isOPALt() {

    if (!doHDF5_m) return false;

    char opalFlavour[128];
    h5_int64_t rc = H5ReadStepAttribString(H5file_m, "OPAL_flavour", opalFlavour);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    return (std::string(opalFlavour)==std::string("opal-t"));
}

/** \brief
 *
 *
 *
 */
void DataSink::setOPALcycl() {

    if (!doHDF5_m) return;

    string OPALFlavour("opal-cycl");
    h5_int64_t rc = H5WriteStepAttribString(H5file_m, "OPAL_flavour", OPALFlavour.c_str());
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
}

/** \brief
 *
 *
 *
 */
void DataSink::setOPALt() {

    if (!doHDF5_m) return;

    string OPALFlavour("opal-t");
    h5_int64_t rc = H5WriteStepAttribString(H5file_m, "OPAL_flavour", OPALFlavour.c_str());
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
}

/** \brief
 *  delete the last 'numberOfLines' lines of the file 'fileName'
 */
void DataSink::rewindLines(const std::string &fileName, size_t numberOfLines) const {
    std::string line;
    std::queue<std::string> allLines;
    std::fstream fs;

    fs.open (fileName.c_str(), std::fstream::in);

    if (!fs.is_open()) return;

    while (getline(fs, line)) {
        allLines.push(line);
    }
    fs.close();


    fs.open (fileName.c_str(), std::fstream::out);

    if (!fs.is_open()) return;

    while (allLines.size() > numberOfLines) {
        fs << allLines.front() << "\n";
        allLines.pop();
    }
    fs.close();
}

/***************************************************************************
 * $RCSfile: DataSink.cpp,v $   $Author: adelmann $
 * $Revision: 1.3 $   $Date: 2004/06/02 19:38:54 $
 ***************************************************************************/