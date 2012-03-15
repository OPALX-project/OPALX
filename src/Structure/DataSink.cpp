// ------------------------------------------------------------------------
// $RCSfile: DataSink.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.3 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: DataSink
//
// ------------------------------------------------------------------------
//
// Revision History:
// $Date: 2004/06/02 19:38:54 $
// $Author: adelmann $
// $Log: DataSink.cpp,v $
// Revision 1.3  2004/06/02 19:38:54  adelmann
// add tabs in the stat file
//
// Revision 1.2  2003/12/12 10:46:39  adelmann
// Remove old stuff
//
// Revision 1.1.1.2  2003/07/04 12:52:00  adelmann
// July 4 2003
//
// Revision 3.2  2001/09/08 16:34:08  adelmann
// Add writeout of: moments, dispersion
//
// Revision 3.1  2001/08/30 11:28:26  adelmann
// Remove .madprogress and write to gnuplot
// add courant snyder parametres to statfile
//
// Revision 3.0  2001/08/22 14:41:33  adelmann
// The stable Version
//
// Revision 2.17  2001/08/11 05:30:22  adelmann
// Production Version
//
// Revision 1.6  2001/08/11 05:21:47  adelmann
// August 11 2001 V2.17
//
// Add #-tag to string in statfile
//
// Revision 1.5  2001/07/05 20:52:53  adelmann
// Mad9p V2.12
//
// Used reduction on one processor only -> program hangs
//
// Revision 1.4  2001/02/23 12:28:34  adelmann
// Add space between the column of the STAMARKER
//
// Revision 1.3  2001/02/11 17:14:27  adelmann
// Remove .mad9pprogress
//
// Revision 1.2  2001/01/26 12:42:49  adelmann
// This is MAD9p Version: Tue Jan 25 2001 V1.28
//
// Revision 1.1.1.1  2000/11/30 20:29:48  adelmann
// g++ and KCPP
//
// Revision 1.5  2000/11/09 15:24:24  adelmann
// - no output in .progress file anymore
//   the function is NOT removed so one can use
//   this to dump some debug information
//
// - collect all reference data from PartData
//
// Revision 1.4  2000/11/05 04:24:47  adelmann
// Add mean values of momentum to the statistics file LANL Nov 2000
//
// Revision 1.3  2000/10/26 05:17:58  adelmann
// Remove DX stuff and add timestamp and title
//
// Revision 1.2  2000/08/10 10:56:34  adelmann
// Some cleanup and add a option dx (data explorer) !
// The new printall script extracts data from this stat file format
//
// Revision 1.1.1.1  2000/07/14 07:20:54  adelmann
// linux version Fri Jul 14 09:15:27 CEST 2000
//
// Revision 1.1.1.1  2000/05/20 11:13:58  adelmann
// Initial working version, thick elements, without: field dump and collimators
//
// Revision 1.3  2000/01/28 07:22:40  adelmann
// Fixed some bugs with the Mad9pOutput
//
// Revision 1.2  2000/01/27 14:13:27  adelmann
// - Add  bunch->dataSink_m.saveStatDataGnuplotFormat( . )
//        DiscParticle write out
//        updateDotProgress
//
// Revision 1.1.1.1  2000/01/06 07:33:27  adelmann
// linux version works with gcc 991007 V2.96
//
// Revision 2.2  1999/10/29 05:02:05  adelmann
// *** empty log message ***
//
// Revision 2.1  1999/10/27 06:36:59  adelmann
// SGI-LINUX g++, with RF-Gap, REVSCATTER and read distribution from file
//
// Revision 1.1.1.1  1999/10/26 04:22:18  adelmann
// Classic 2.1 (with p)
//
// Revision 1.1.1.1  1999/10/26 04:14:36  adelmann
// Classic 2.1 (with p)
//
//
// ------------------------------------------------------------------------

#include "DataSink.h"
#include "config.h"
#include "AbstractObjects/OpalData.h"
#include "Physics/Physics.h"
#include "Utilities/Options.h"
#include "ValueDefinitions/RealVariable.h"
#include "Fields/Fieldmap.hh"
#include <iomanip>
#include <string>
#include <hdf5.h>
#include "H5Part.h"
extern "C" {
#include "H5Block.h"
#include "H5BlockTypes.h"
}

extern Inform *gmsg;

using namespace Options;
using namespace std;
//using Physics::m_e;

DataSink::DataSink()
{
    /// Constructor steps:

    /// Get timers from IPPL.
    H5PartTimer_m = IpplTimings::getTimer("H5PartTimer");
    StatMarkerTimer_m = IpplTimings::getTimer("StatMarkerTimer");

    /// Set file write flags to true. These will be set to false after first
    /// write operation.
    firstWriteToStat_m = true;
    firstWriteH5part_m = true;

    /// Define file names.
    string fn = OPAL.getInputFn();

    int pos=fn.find(string("."), 0);
    fn.erase(pos, fn.size() - pos);

    statFileName_m = fn + string(".stat");
    lBalFileName_m = fn + string(".lbal");

    fn += string(".h5");

    /// Open H5 file. Check that it opens correctly.
#ifdef PARALLEL_IO
    H5file_m=H5PartOpenFileParallel(fn.c_str(),H5PART_WRITE,MPI_COMM_WORLD);
#else
    H5file_m=H5PartOpenFile(fn.c_str(),H5PART_WRITE);
#endif

    if(!H5file_m) {
        ERRORMSG("h5 file open failed: exiting!" << endl);
        exit(0);
    }

    /// Write file attributes.
    writeH5FileAttributes();

    /// Set current record/time step to 0.
    H5call_m = 0;

    /// Determine the reference particle position shift.
    sshift_m  = 0.0;
    RealVariable *ar = dynamic_cast<RealVariable *>(OPAL.find("SSHIFT"));
    if (ar) {
        sshift_m = ar->getReal();
        *gmsg << "SSHIFT " << sshift_m << endl;
    }

}

DataSink::DataSink(int restartStep)
{
    /// Constructor steps:

    /// Get timers from IPPL.
    H5PartTimer_m = IpplTimings::getTimer("H5PartTimer");
    StatMarkerTimer_m = IpplTimings::getTimer("StatMarkerTimer");

    /// Set file write flags. Since this is a restart, assume H5 file is old.
    firstWriteToStat_m = true;
    firstWriteH5part_m = false;

    /// Get file name root.
    string fn = OPAL.getInputFn();

    int pos = fn.find(string(".in"), 0);
    fn.erase(pos, fn.size() - pos);

    statFileName_m = fn + string(".stat");
    lBalFileName_m = fn + string(".lbal");

    /// Test if .stat and .lbal files exist. If they do, new data will be appended.
    ofstream statFile(statFileName_m.c_str(), ios::in);
    ofstream lBalFile(lBalFileName_m.c_str(), ios::in);

    if (statFile.is_open()) {
        // File exists so we append data to end.
        firstWriteToStat_m = false;
        statFile.close();
        *gmsg << "Appending statistical data to existing data file: " << statFileName_m << endl;
    } else {
        statFile.clear();
        *gmsg << "Creating new file for statistical data: " << statFileName_m << endl;
    }

    if (lBalFile.is_open()) {
        // File exists so we append data to end.
        lBalFile.close();
        *gmsg << "Appending load balance data to existing data file: " << lBalFileName_m << endl;
    } else {
        lBalFile.clear();
        *gmsg << "Creating new file for load balance data: " << lBalFileName_m << endl;
    }

    /// Find if OPAL will append data to existing H5 file or create new one.
    if (OPAL.hasRestartFile()) {

        // Define H5 file name.
        fn=OPAL.getRestartFileName();
        if (fn.find(string("_Part"),0) == string::npos){
            int pos = fn.find(string(".h5"),0);
            fn.replace( pos, fn.length(), "_Part00002.h5");
        }
        else{

            int pos=fn.find(string("_Part"),0);
            string numstr=fn.substr(pos+5,5);
            int numint = atoi(numstr.c_str())+1;
            if (numint > 99999){
                ERRORMSG("Maximum \"h5\" file number is limited to 99999 !" << endl);
                exit(0);
            }

            char numCStr[5];
            sprintf(numCStr, "%05d", numint);
            string numstrnew(numCStr);

            fn.replace( pos+5, 5, numstrnew);
        }
        *gmsg<<"New phase space data will dump into the new file: " << fn << endl;

        // Open H5 file. Check that it opens correctly.
#ifdef PARALLEL_IO
        H5file_m=H5PartOpenFileParallel(fn.c_str(),H5PART_WRITE,MPI_COMM_WORLD);
#else
        H5file_m=H5PartOpenFile(fn.c_str(),H5PART_WRITE);
#endif

        if(!H5file_m) {
            ERRORMSG("h5 file open failed: exiting!" << endl);
            exit(0);
        }

        // Write attributes to new H5 file.
        writeH5FileAttributes();

        // This is a new file.
        firstWriteH5part_m = true;

        // Use dump frequency from old H5 file.
        int dumpfreq = 0;
        H5PartReadFileAttrib(H5file_m, "dump frequency", &dumpfreq);
        OPAL.setRestartDumpFreq(dumpfreq);

        // Set record/time step to 0.
        H5call_m = 0;

    } else {

        // Define file name.
        fn += string(".h5");

#ifdef PARALLEL_IO
        H5file_m=H5PartOpenFileParallel(fn.c_str(),H5PART_APPEND,MPI_COMM_WORLD);
#else
        H5file_m=H5PartOpenFile(fn.c_str(),H5PART_APPEND);
#endif

        INFOMSG(fn << endl);

        if(!H5file_m) {
            ERRORMSG("h5 file open failed: exiting!" << endl);
            exit(0);
        }

        int numStepsInFile = H5PartGetNumSteps(H5file_m);
        if (restartStep == -1) {
            restartStep = numStepsInFile;
        }

        if(numStepsInFile < restartStep || numStepsInFile > restartStep) {
            ERRORMSG("Must use last step in H5 restart file for restart" << endl
                     << "To solve this problem: 1) Use last step in H5 restart file" << endl
                     << "                       2) Use a value of -1 for restart step" << endl
                     << "                       3) Write new data to new H5 file" << endl
                     << endl << "Exiting!" << endl);
            exit(0);
        }

        // Use same dump frequency.
        int dumpfreq = 0;
        H5PartReadFileAttrib(H5file_m, "dump frequency", &dumpfreq);
        OPAL.setRestartDumpFreq(dumpfreq);

        // Set current record/time step to restart step.
        H5call_m = restartStep;
    }

    /// Determine reference particle position shift.
    sshift_m  = 0.0;
    RealVariable *ar = dynamic_cast<RealVariable *>(OPAL.find("SSHIFT"));
    if (ar) {
        sshift_m = ar->getReal();
        *gmsg << "SSHIFT " << sshift_m << endl;
    }
}

DataSink::~DataSink()
{
    H5PartCloseFile(H5file_m);
    Ippl::Comm->barrier();
}

void DataSink::writeH5FileAttributes() {

    /// Function steps:

    /// Write file attributes to describe phase space to H5 file.
    stringstream OPAL_version;
    OPAL_version << PACKAGE_STRING << " rev. " << SVN_VERSION;
    H5PartWriteFileAttribString(H5file_m,"OPAL_version",OPAL_version.str().c_str());

    H5PartWriteFileAttribString(H5file_m,"tUnit","s");
    H5PartWriteFileAttribString(H5file_m,"xUnit","m");
    H5PartWriteFileAttribString(H5file_m,"yUnit","m");
    H5PartWriteFileAttribString(H5file_m,"zUnit","m");
    H5PartWriteFileAttribString(H5file_m,"pxUnit","#beta#gamma");
    H5PartWriteFileAttribString(H5file_m,"pyUnit","#beta#gamma");
    H5PartWriteFileAttribString(H5file_m,"pzUnit","#beta#gamma");
    H5PartWriteFileAttribString(H5file_m,"idUnit","1");
    H5PartWriteFileAttribString(H5file_m,"SPOSUnit","m");
    H5PartWriteFileAttribString(H5file_m,"TIMEUnit","s");
    H5PartWriteFileAttribString(H5file_m,"#gammaUnit","1");
    H5PartWriteFileAttribString(H5file_m,"ENERGYUnit","MeV");
    H5PartWriteFileAttribString(H5file_m,"#varepsilonUnit","m rad");
    H5PartWriteFileAttribString(H5file_m,"#varepsilonrUnit","m rad");

    H5PartWriteFileAttribString(H5file_m,"#varepsilon-geomUnit","m rad");

    H5PartWriteFileAttribString(H5file_m,"#sigmaUnit"," ");
    H5PartWriteFileAttribString(H5file_m,"RMSXUnit","m");
    H5PartWriteFileAttribString(H5file_m,"RMSRUnit","m");
    H5PartWriteFileAttribString(H5file_m,"RMSPUnit","#beta#gamma");

    H5PartWriteFileAttribString(H5file_m,"maxdEUnit","MeV");
    H5PartWriteFileAttribString(H5file_m,"max#phiUnit","deg");

    H5PartWriteFileAttribString(H5file_m,"phizUnit","deg");
    H5PartWriteFileAttribString(H5file_m,"enezUnit","keV");

    H5PartWriteFileAttribString(H5file_m,"mpart","GeV");
    H5PartWriteFileAttribString(H5file_m,"qi","C");

    /// Write file attributes to describe fields of head/ref particle/tail.
    H5PartWriteFileAttribString(H5file_m,"spos-headUnit","m");
    H5PartWriteFileAttribString(H5file_m,"E-headUnit","MV/m");
    H5PartWriteFileAttribString(H5file_m,"B-headUnit","T");
    H5PartWriteFileAttribString(H5file_m,"spos-refUnit","m");
    H5PartWriteFileAttribString(H5file_m,"E-refUnit","MV/m");
    H5PartWriteFileAttribString(H5file_m,"B-refUnit","T");
    H5PartWriteFileAttribString(H5file_m,"spos-tailUnit","m");
    H5PartWriteFileAttribString(H5file_m,"E-tailUnit","MV/m");
    H5PartWriteFileAttribString(H5file_m,"B-tailUnit","T");
    H5PartWriteFileAttribString(H5file_m,"StepUnit"," ");
    H5PartWriteFileAttribString(H5file_m,"TrackStepUnit"," ");
    H5PartWriteFileAttribString(H5file_m,"NumBunchUnit"," ");
    H5PartWriteFileAttribString(H5file_m,"NumPartUnit"," ");
    H5PartWriteFileAttribString(H5file_m,"GammaBinUnit"," ");
    H5PartWriteFileAttribString(H5file_m,"RefPartRUnit","m");
    H5PartWriteFileAttribString(H5file_m,"RefPartPUnit","#beta#gamma");
    H5PartWriteFileAttribString(H5file_m,"SteptoLastInjUnit","");

    /// Write file dump frequency.
    h5part_int64_t dumpfreq = Options::psDumpFreq;
    H5PartWriteFileAttrib(H5file_m, "dump frequency", H5PART_INT64, &dumpfreq, 1);

    /// Reset first write flag.
    firstWriteH5part_m = false;


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

    if ( H5file_m->timegroup >= 0 ) {
        herr = H5Gclose ( H5file_m->timegroup );
        H5file_m->timegroup = -1;
    }

    if (Ippl::myNode() == 0) {
        struct stat st;
        off_t fsize;
        if (stat(OPAL.getInputFn().c_str(), &st) == 0) {
            fsize = st.st_size;
        }
        ContentLength = fsize / sizeof(char);
        FileContent = new char[ContentLength];

        filebuf inputFileBuffer;
        inputFileBuffer.open(OPAL.getInputFn().c_str(), ios::in);
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

    MPI_Bcast ( &ContentLength, 
                1, 
                MPI_LONG_LONG_INT, 
                0, 
                MPI_COMM_WORLD);
    
    herr = H5Gget_objinfo( H5file_m->file, group_name, 1, NULL );
    if (herr >= 0) { // there exists a group 'INPUT'
        delete[] FileContent;
        return;
    }
    
    group_id = H5Gcreate ( H5file_m->file, group_name, 0 );

    shape = H5Screate_simple(1, &ContentLength, &ContentLength);
    dataset_id = H5Dcreate( group_id,
                            dataset_name,
                            H5T_NATIVE_CHAR,
                            shape,
                            H5P_DEFAULT );
    H5Sclose ( shape );

    diskshape = H5Dget_space ( dataset_id );
    H5Sselect_hyperslab ( diskshape, 
                          H5S_SELECT_SET, 
                          &start, 
                          NULL, 
                          &write_length, 
                          NULL );

    memshape = H5Screate_simple(1, &write_length, &dmax);

    herr = H5Dwrite ( dataset_id,
                      H5T_NATIVE_CHAR,
                      memshape,
                      diskshape,
                      H5file_m->xfer_prop,
                      FileContent );

    H5Sclose ( memshape );
    H5Dclose ( dataset_id );
    H5Sclose ( diskshape );
    H5Gclose ( group_id );

    delete[] FileContent;
}

int DataSink::storeFieldmaps() {
    vector<string> fieldmap_list = Fieldmap::getListFieldmapNames();
    vector<file_size_name> fieldmap_list2;
    off_t fsize;
    off_t total_fsize = 0;
    struct stat st;
    int Nf = 0;
    int Np = Ippl::getNodes();
    if (Np > 20) Np = 20; // limit the number of processors that write data

    for (vector<string>::const_iterator it = fieldmap_list.begin(); it != fieldmap_list.end(); ++ it) {
        if (stat((*it).c_str(), &st) == 0) {
            fsize = st.st_size;
            ++ Nf;
        } else {
            continue;
        }
        fieldmap_list2.push_back(file_size_name((*it),fsize));
    }
    sort(fieldmap_list2.begin(), fieldmap_list2.end(), file_size_name::SortAsc);
//     for (int i = Nf - 1; i > -1; -- i) {
//         if (i % Np == Ippl::myNode()) {
//             cout << Ippl::myNode() << "\t" << fieldmap_list2[i].file_name_m << "\t" << fieldmap_list2[i].file_size_m << " bytes" << endl;
//             total_fsize += fieldmap_list2[i].file_size_m;
//         }
//     }

    char group_name[] = "INPUT/FIELDMAPS";
    int ContentLength;
    char* FileContent;
    herr_t herr;

//     herr = H5Gget_objinfo( H5file_m->file, group_name, 0 );
//     if (herr >= 0) {
//         return -1;
//     }

    if (Ippl::myNode() < Np) {
        for (int i = 0; i < Nf / Np; ++ i) {
            int lid = i * Np + Ippl::myNode();
            string filename = fieldmap_list2[lid].file_name_m;
            int ContentLength = fieldmap_list2[lid].file_size_m / sizeof(char);

            filebuf inputFileBuffer;
            inputFileBuffer.open(filename.c_str(), ios::in);
            istream inputFile(&inputFileBuffer);

            FileContent = new char[ContentLength];

            inputFile.get(FileContent, ContentLength, '\0');

            inputFileBuffer.close();
            delete[] FileContent;
        }
    }
    //group_id = H5Gcreate ( H5file_m->file, group_name, 0 );
    
    
}

void DataSink::writePhaseSpace(PartBunch &beam, Vector_t FDext[], double sposHead, double sposRef, double sposTail) {

    /// Function steps:

    /// Start timer.
    IpplTimings::startTimer(H5PartTimer_m);

    /// Calculate beam statistical parameters etc. Put them in the right format
    /// for H5 file write.
    beam.calcBeamParameters();

    double               actPos   = beam.get_sPos() +   sshift_m ;
    double               t        = beam.getT();
    Vektor< double, 3 >  rmin     = beam.get_origin();
    Vektor< double, 3 >  rmax     = beam.get_maxExtend();
    Vektor< double, 3 >  centroid = beam.get_centroid();

    size_t nTot                   = beam.getTotalNum();
    size_t nLoc                   = beam.getLocalNum();

    Vektor< double, 3 >  maxP(0.0);
    Vektor< double, 3 >  minP(0.0);

    Vektor<double, 3 > xsigma = beam.get_rrms();
    Vektor<double, 3 > psigma = beam.get_prms();
    Vektor<double, 3 > vareps = beam.get_norm_emit();
    Vektor<double, 3 > geomvareps = beam.get_emit();
    Vektor<double, 3 > RefPartR = beam.RefPart_R;
    Vektor<double, 3 > RefPartP = beam.RefPart_P;

    double meanEnergy = beam.get_meanEnergy();

    double sigma = ((xsigma[0]*xsigma[0])+(xsigma[1]*xsigma[1])) /
        (2.0*beam.get_gamma()*17.0e3*((geomvareps[0]*geomvareps[0]) + (geomvareps[1]*geomvareps[1])));

    beam.get_PBounds(minP,maxP);

    void *varray = malloc(nLoc*sizeof(double));
    double *farray = (double*)varray;
    h5part_int64_t *larray = (h5part_int64_t *)varray;

    ///Get the particle decomposition from all the compute nodes.
    size_t *locN = (size_t *) malloc(Ippl::getNodes()*sizeof(size_t));
    size_t  *globN = (size_t*) malloc(Ippl::getNodes()*sizeof(size_t));

    for(int i=0; i<Ippl::getNodes(); i++) {
        globN[i] = locN[i]=0;
    }
    locN[Ippl::myNode()] = nLoc;
    reduce(locN, locN + Ippl::getNodes(), globN, OpAddAssign());

    /// Set current record/time step.
    H5PartSetStep(H5file_m,H5call_m);
    H5PartSetNumParticles(H5file_m,nLoc);

    /// Write statistical data.
    H5PartWriteStepAttrib(H5file_m,"SPOS",     H5T_NATIVE_DOUBLE,&actPos,1);
    H5PartWriteStepAttrib(H5file_m,"#sigma",   H5T_NATIVE_DOUBLE,&sigma,1);
    H5PartWriteStepAttrib(H5file_m,"RMSX",     H5T_NATIVE_DOUBLE,&xsigma,3); //sigma
    H5PartWriteStepAttrib(H5file_m,"RMSP",     H5T_NATIVE_DOUBLE,&psigma,3); //sigma
    H5PartWriteStepAttrib(H5file_m,"maxX",     H5T_NATIVE_DOUBLE,&rmax,3);
    H5PartWriteStepAttrib(H5file_m,"minX",     H5T_NATIVE_DOUBLE,&rmin,3);
    H5PartWriteStepAttrib(H5file_m,"maxP",     H5T_NATIVE_DOUBLE,&maxP,3);
    H5PartWriteStepAttrib(H5file_m,"minP",     H5T_NATIVE_DOUBLE,&minP,3);
    H5PartWriteStepAttrib(H5file_m,"centroid", H5T_NATIVE_DOUBLE,&centroid,3);
    H5PartWriteStepAttrib(H5file_m,"TIME",     H5T_NATIVE_DOUBLE,&t,1);
    H5PartWriteStepAttrib(H5file_m,"ENERGY",   H5T_NATIVE_DOUBLE,&meanEnergy,1);
    H5PartWriteStepAttrib(H5file_m,"RefPartR", H5T_NATIVE_DOUBLE,&RefPartR,3);
    H5PartWriteStepAttrib(H5file_m,"RefPartP", H5T_NATIVE_DOUBLE,&RefPartP,3);

    /// Write head/reference particle/tail field information.
    FDext[1] *= 1e-6;
    FDext[3] *= 1e-6;
    FDext[5] *= 1e-6;
    H5PartWriteStepAttrib(H5file_m,"B-head", H5T_NATIVE_DOUBLE,&FDext[0],3);
    H5PartWriteStepAttrib(H5file_m,"E-head", H5T_NATIVE_DOUBLE,&FDext[1],3);
    H5PartWriteStepAttrib(H5file_m,"B-ref",  H5T_NATIVE_DOUBLE,&FDext[2],3);
    H5PartWriteStepAttrib(H5file_m,"E-ref",  H5T_NATIVE_DOUBLE,&FDext[3],3);
    H5PartWriteStepAttrib(H5file_m,"B-tail", H5T_NATIVE_DOUBLE,&FDext[4],3);
    H5PartWriteStepAttrib(H5file_m,"E-tail", H5T_NATIVE_DOUBLE,&FDext[5],3);

    sposHead +=   sshift_m ;
    sposRef  +=   sshift_m ;
    sposTail +=   sshift_m ;

    H5PartWriteStepAttrib(H5file_m,"spos-head", H5T_NATIVE_DOUBLE, &sposHead, 1);
    H5PartWriteStepAttrib(H5file_m,"spos-ref",  H5T_NATIVE_DOUBLE, &sposRef, 1);
    H5PartWriteStepAttrib(H5file_m,"spos-tail", H5T_NATIVE_DOUBLE, &sposTail, 1);

    /// Write number of compute nodes.
    H5PartWriteStepAttrib(H5file_m,"nloc",H5T_NATIVE_INT64, globN, Ippl::getNodes());

    ///Write particle mass and charge per particle. (Consider making these
    ///file attributes.)
    double mpart = 1.0e-9*beam.getM();
    double    qi = beam.getChargePerParticle();

    H5PartWriteStepAttrib(H5file_m,"mpart", H5T_NATIVE_DOUBLE, &mpart, 1);
    H5PartWriteStepAttrib(H5file_m,"qi", H5T_NATIVE_DOUBLE, &qi, 1);

    //Not provided by PartBunch at the moment
    //H5PartWriteStepAttrib(H5file_m,"#varepsilonr", H5T_NATIVE_DOUBLE, ,1);

    /// Write normalized emittance.
    H5PartWriteStepAttrib(H5file_m,"#varepsilon", H5T_NATIVE_DOUBLE, &vareps, 3);

    /// Write geometric emittance.
    H5PartWriteStepAttrib(H5file_m,"#varepsilon-geom", H5T_NATIVE_DOUBLE, &geomvareps, 3);

    /// Write beam phase space.
    for (size_t i=0; i<nLoc;i++)
        farray[i] =  beam.R[i](0);
    H5PartWriteDataFloat64(H5file_m,"x",farray);

    for (size_t i=0; i<nLoc;i++)
        farray[i] =  beam.R[i](1);
    H5PartWriteDataFloat64(H5file_m,"y",farray);

    for (size_t i=0; i<nLoc;i++)
        farray[i] =  beam.R[i](2);
    H5PartWriteDataFloat64(H5file_m,"z",farray);

    for (size_t i=0; i<nLoc;i++)
        farray[i] =  beam.P[i](0);
    H5PartWriteDataFloat64(H5file_m,"px",farray);

    for (size_t i=0; i<nLoc;i++)
        farray[i] =  beam.P[i](1);
    H5PartWriteDataFloat64(H5file_m,"py",farray);

    for (size_t i=0; i<nLoc;i++)
        farray[i] =  beam.P[i](2);
    H5PartWriteDataFloat64(H5file_m,"pz",farray);

    //for (size_t i=0; i<nLoc;i++)
    //  farray[i] =  phiz;
    //H5PartWriteDataFloat64(H5file_m,"phiz",farray);

    //  H5PartWriteDataFloat64(H5file_m,"enez",energy);

    /// Write particle id numbers.
    for (size_t i = 0; i < nLoc; i++)
        larray[i] =  beam.ID[i];
    H5PartWriteDataInt64(H5file_m,"id",larray);

    /// Write space charge field map if asked for.
    if (Options::rhoDump) {

        h5part_int64_t l[6];

        NDIndex<3> idx = beam.getFieldLayout().getLocalNDIndex();
        NDIndex<3> elem;
        h5part_int64_t herr = H5BlockDefine3DFieldLayout (H5file_m,
                                                          idx[0].min(), idx[0].max(),
                                                          idx[1].min(), idx[1].max(),
                                                          idx[2].min(), idx[2].max());

        if (herr < 0)
            *gmsg << "H5BlockDefine3DFieldLayout err " << herr << endl;

        h5part_float64_t *data = (h5part_float64_t *) malloc ( (idx[0].max()+1)  * (idx[1].max()+1) * (idx[2].max()+1) * sizeof ( *data ) );

        int ii = 0;
        // h5block uses the fortran convention of storing data:
        // INTEGER, DIMENSION(2,3) :: a
        // => {a(1,1), a(2,1), a(1,2), a(2,2), a(1,3), a(2,3)}
        for (int i=idx[2].min(); i<=idx[2].max(); ++i) {
            for (int j=idx[1].min(); j<=idx[1].max(); ++j) { 
                for (int k=idx[0].min(); k<=idx[0].max(); ++k) {
                    data[ii] = beam.getRho(k,j,i);
                    ii++;
                }
            }
        }

        herr = H5Block3dWriteScalarField (H5file_m, "rho", data );

        if (herr < 0)
            *gmsg << "H5Block3dWriteScalarField err " << herr << endl;

        /// Need this to align particles and fields when writing space charge map.
        herr = H5Block3dSetFieldOrigin (H5file_m,"rho",
                                        (h5part_float64_t)beam.get_origin()(0),
                                        (h5part_float64_t)beam.get_origin()(1),
                                        (h5part_float64_t)beam.get_origin()(2));

        herr = H5Block3dSetFieldSpacing (H5file_m,"rho",
                                         (h5part_float64_t)beam.get_hr()(0),
                                         (h5part_float64_t)beam.get_hr()(1),
                                         (h5part_float64_t)beam.get_hr()(2));
        if(data)
            free(data);
    }

    H5Fflush(H5file_m->file,H5F_SCOPE_GLOBAL);

    /// Step record/time step index.
    H5call_m++;

    if(varray)
        free(varray);

    /// %Stop timer.
    IpplTimings::stopTimer(H5PartTimer_m);
}



int DataSink::writePhaseSpace_cycl(PartBunch &beam, Vector_t FDext[]) {

    using Physics::m_p;
    IpplTimings::startTimer(H5PartTimer_m);

    beam.calcBeamParameters_cycl();

    // double               actPos   = beam.get_sPos();
    double               t        = beam.getT();

    Vektor< double, 3 >  rmin     = beam.get_origin();
    Vektor< double, 3 >  rmax     = beam.get_maxExtend();
    Vektor< double, 3 >  centroid = beam.get_centroid();

    size_t nTot                   = beam.getTotalNum();
    size_t nLoc                   = beam.getLocalNum();

    Vektor< double, 3 >  maxP(0.0);
    Vektor< double, 3 >  minP(0.0);


    
    Vektor<double, 3 > xsigma = beam.get_rrms();
    Vektor<double, 3 > psigma = beam.get_prms();
    Vektor<double, 3 > geomvareps = beam.get_emit();
    Vektor<double, 3 > vareps = beam.get_norm_emit();
    double meanEnergy = beam.get_meanEnergy();
    double meanGamma = meanEnergy/ (m_p*1000.0) + 1;

    beam.get_PBounds(minP,maxP);

    void *varray = malloc(nLoc*sizeof(double));
    double *farray = (double*)varray;
    h5part_int64_t *larray = (h5part_int64_t *)varray;

    double  pathLength = beam.getLPath();
    h5part_int64_t trackStep =(h5part_int64_t)beam.getTrackStep();
    h5part_int64_t numBunch =(h5part_int64_t)beam.getNumBunch();
    h5part_int64_t SteptoLastInj = (h5part_int64_t)beam.getSteptoLastInj();

    h5part_int64_t numPart[numBunch];

    // used for classify particles in restart multibunch run
    double gammaBin[numBunch];
    // calcualte gamma for each energy bin
    if ( beam.weHaveBins() )
        {

            beam.calcGammas_cycl();
            for (int ii=0; ii<beam.getLastemittedBin(); ii++ )
                {
                    numPart[ii] = beam.getNumPartInBin(ii);
                    gammaBin[ii]= beam.getBinGamma(ii);
                    *gmsg<<"Dump Bin: " <<ii<<", gamma: "<<gammaBin[ii]<<", Particles: "<<numPart[ii]<<endl;
                }

            // for multibunch, "meanenergy" "meanBin" represtants
            // the mean energy and gamma of the bin with BinID=0.
            meanEnergy = ( gammaBin[0] -1 ) * m_p * 1000.0;
            meanGamma  = gammaBin[0];
        }

    /* ------------------------------------------------------------------------
       Get the particle decomposition from all the nodes
    */
    size_t *locN = (size_t *) malloc(Ippl::getNodes()*sizeof(size_t));
    size_t  *globN = (size_t*) malloc(Ippl::getNodes()*sizeof(size_t));

    for(int i=0; i<Ippl::getNodes(); i++) {
        globN[i] = locN[i]=0;
    }
    locN[Ippl::myNode()] = nLoc;
    reduce(locN, locN + Ippl::getNodes(), globN, OpAddAssign());

    /* ------------------------------------------------------------------------ */
    H5PartSetStep(H5file_m,H5call_m);
    H5PartSetNumParticles(H5file_m,nLoc);

    /* write scalar data i.e the header */
    H5PartWriteStepAttrib(H5file_m,"SPOS",     H5T_NATIVE_DOUBLE,&pathLength,1);
    H5PartWriteStepAttrib(H5file_m,"Step",     H5T_NATIVE_INT64,&H5call_m,1);
    H5PartWriteStepAttrib(H5file_m,"RMSX",     H5T_NATIVE_DOUBLE,&xsigma,3); //sigma
    H5PartWriteStepAttrib(H5file_m,"RMSP",     H5T_NATIVE_DOUBLE,&psigma,3); //sigma
    H5PartWriteStepAttrib(H5file_m,"maxX",     H5T_NATIVE_DOUBLE,&rmax,3);
    H5PartWriteStepAttrib(H5file_m,"minX",     H5T_NATIVE_DOUBLE,&rmin,3);
    H5PartWriteStepAttrib(H5file_m,"maxP",     H5T_NATIVE_DOUBLE,&maxP,3);
    H5PartWriteStepAttrib(H5file_m,"minP",     H5T_NATIVE_DOUBLE,&minP,3);
    H5PartWriteStepAttrib(H5file_m,"centroid", H5T_NATIVE_DOUBLE,&centroid,3);
    H5PartWriteStepAttrib(H5file_m,"TIME",     H5T_NATIVE_DOUBLE,&t,1);
    H5PartWriteStepAttrib(H5file_m,"ENERGY",   H5T_NATIVE_DOUBLE,&meanEnergy,1);
    H5PartWriteStepAttrib(H5file_m,"B-ref",    H5T_NATIVE_DOUBLE,&FDext[0],3);
    H5PartWriteStepAttrib(H5file_m,"E-ref",    H5T_NATIVE_DOUBLE,&FDext[1],3);
    H5PartWriteStepAttrib(H5file_m,"nloc",     H5T_NATIVE_INT64, globN, Ippl::getNodes());
    H5PartWriteStepAttrib(H5file_m,"#varepsilon",      H5T_NATIVE_DOUBLE, &vareps, 3);  //unnormalized
    H5PartWriteStepAttrib(H5file_m,"#varepsilon-geom", H5T_NATIVE_DOUBLE, &geomvareps, 3);//normalized
    H5PartWriteStepAttrib(H5file_m,"TrackStep",        H5T_NATIVE_INT64, &trackStep,1);
    H5PartWriteStepAttrib(H5file_m,"NumBunch",         H5T_NATIVE_INT64, &numBunch,1);
    H5PartWriteStepAttrib(H5file_m,"GammaBin",         H5T_NATIVE_DOUBLE, &gammaBin, numBunch);
    H5PartWriteStepAttrib(H5file_m,"SteptoLastInj",        H5T_NATIVE_INT64, &SteptoLastInj,1);

    for (size_t i=0; i<nLoc;i++)
        farray[i] =  beam.R[i](0);
    H5PartWriteDataFloat64(H5file_m,"x",farray);

    for (size_t i=0; i<nLoc;i++)
        farray[i] =  beam.R[i](1);
    H5PartWriteDataFloat64(H5file_m,"y",farray);

    for (size_t i=0; i<nLoc;i++)
        farray[i] =  beam.R[i](2);
    H5PartWriteDataFloat64(H5file_m,"z",farray);

    for (size_t i=0; i<nLoc;i++)
        farray[i] =  beam.P[i](0);
    H5PartWriteDataFloat64(H5file_m,"px",farray);

    for (size_t i=0; i<nLoc;i++)
        farray[i] =  beam.P[i](1);
    H5PartWriteDataFloat64(H5file_m,"py",farray);

    for (size_t i=0; i<nLoc;i++)
        farray[i] =  beam.P[i](2);
    H5PartWriteDataFloat64(H5file_m,"pz",farray);

    for (size_t i=0; i<nLoc;i++)
        larray[i] =  beam.ID[i];
    H5PartWriteDataInt64(H5file_m,"id",larray);

#ifdef H5BLOCKSAVE
    h5part_int64_t l[6];

    NDIndex<3> idx = getFieldLayout().getLocalNDIndex();
    NDIndex<3> elem;
    h5part_int64_t herr = H5BlockDefine3DFieldLayout (
                                                      H5file_m,
                                                      idx[0].min(), idx[0].max(),
                                                      idx[1].min(), idx[1].max(),
                                                      idx[2].min(), idx[2].max());
    if (herr < 0)
        gmsg << "H5BlockDefine3DFieldLayout err " << herr << endl;

    h5part_float64_t *data = (h5part_float64_t *) malloc ( (idx[0].max()+1)  * (idx[1].max()+1) * (idx[2].max()+1) * sizeof ( *data ) );

    int ii = 0;
    // h5block uses the fortran convention of storing data:
    // INTEGER, DIMENSION(2,3) :: a
    // => {a(1,1), a(2,1), a(1,2), a(2,2), a(1,3), a(2,3)}
    for (int i=idx[2].min(); i<=idx[2].max(); ++i) {
        elem[2] = Index(i,i);
        for (int j=idx[1].min(); j<=idx[1].max(); ++j) {
            elem[1] = Index(j,j);
            for (int k=idx[0].min(); k<=idx[0].max(); ++k) {
                elem[0] = Index(k,k);
                data[ii] = 0.0; // EFDMag_m.localElement(elem);
                ii++;
            }
        }
    }

    herr = H5Block3dWriteScalarField (H5file_m, "EFmag", data );
    if (herr < 0)
        gmsg << "H5Block3dWriteScalarField err " << herr << endl;

    herr = H5Block3dSetFieldSpacing (H5file_m,"EFmag",
                                     (h5part_float64_t)hr_m[0],
                                     (h5part_float64_t)hr_m[1],
                                     (h5part_float64_t)hr_m[2]);

    if (Ippl::myNode() == 0) {
        for (h5part_int64_t p=0; p<Ippl::getNodes();p++) {
            herr = H5Block3dGetPartitionOfProc(H5file_m, p, &l[0], &l[1], &l[2], &l[3], &l[4], &l[5]);
            stringstream lstr;
            lstr << "layout" << p;
            H5BlockWriteFieldAttrib (H5file_m,"EFmag", lstr.str().c_str(), H5PART_INT64,l,6);
        }
    }
    if(data)
        free(data);
#endif

    H5Fflush(H5file_m->file,H5F_SCOPE_GLOBAL);

    H5call_m++;

    if(varray)
        free(varray);

    IpplTimings::stopTimer(H5PartTimer_m);

    return (int)( H5call_m-1 );
}

void DataSink::writeStatData(PartBunch &beam, Vector_t FDext[], double sposHead, double sposRef, double sposTail)
{
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

    if(Ippl::myNode() == 0) {
        if (firstWriteToStat_m) {
            os_statData.open(statFileName_m.c_str(),ios::out);
            os_statData.precision(15);
            os_statData.setf(ios::scientific, ios::floatfield);
            writeSDDSHeader(os_statData);

            os_lBalData.open(lBalFileName_m.c_str(),ios::out);
            os_lBalData.precision(15);
            os_lBalData.setf(ios::scientific, ios::floatfield);
            os_lBalData << "# " << Ippl::getNodes() << endl;

            firstWriteToStat_m = false;
        } else {
            os_statData.open(statFileName_m.c_str(),ios::app);
            os_statData.precision(15);
            os_statData.setf(ios::scientific, ios::floatfield);

            os_lBalData.open(lBalFileName_m.c_str(),ios::app);
            os_lBalData.precision(15);
            os_lBalData.setf(ios::scientific, ios::floatfield);
        }

	FDext[1] *= 1e-6;
	FDext[3] *= 1e-6;
	FDext[5] *= 1e-6;

        os_statData << beam.getT() << setw(pwi) << "\t"                                       // 1
                    << sposRef << setw(pwi) << "\t"                                           // 2

                    << beam.getTotalNum() << setw(pwi) << "\t"                                // 3
                    << beam.getTotalNum() * beam.getChargePerParticle() << setw(pwi) << "\t"  // 4

                    << beam.get_meanEnergy() << setw(pwi) << "\t"                             // 5

                    << beam.get_rrms()(0) << setw(pwi) << "\t"                                // 6
                    << beam.get_rrms()(1) << setw(pwi) << "\t"                                // 7
                    << beam.get_rrms()(2) << setw(pwi) << "\t"                                // 8

                    << beam.get_prms()(0) << setw(pwi) << "\t"                                // 9
                    << beam.get_prms()(1) << setw(pwi) << "\t"                                // 10
                    << beam.get_prms()(2) << setw(pwi) << "\t"                                // 11

                    << beam.get_norm_emit()(0) << setw(pwi) << "\t"                           // 12
                    << beam.get_norm_emit()(1) << setw(pwi) << "\t"                           // 13
                    << beam.get_norm_emit()(2) << setw(pwi) << "\t"                           // 14

                    << beam.get_rmean()(0)  << setw(pwi)<< "\t"                               // 15
                    << beam.get_rmean()(1)  << setw(pwi)<< "\t"                               // 16
                    << beam.get_rmean()(2)  << setw(pwi)<< "\t"                               // 17

                    << beam.get_maxExtend()(0) << setw(pwi) << "\t"                           // 18
                    << beam.get_maxExtend()(1) << setw(pwi) << "\t"                           // 19
                    << beam.get_maxExtend()(2) << setw(pwi) << "\t"                           // 20

            // Write out Courant Snyder parameters.
                    << 0.0  << setw(pwi) << "\t"                                              // 21
                    << 0.0  << setw(pwi) << "\t"                                              // 22

                    << 0.0 << setw(pwi) << "\t"                                               // 23
                    << 0.0 << setw(pwi) << "\t"                                               // 24

            // Write out dispersion.
                    << 0.0 << setw(pwi) << "\t"                                               // 25
                    << 0.0 << setw(pwi) << "\t"                                               // 26
                    << 0.0 << setw(pwi) << "\t"                                               // 27
                    << 0.0 << setw(pwi) << "\t"                                               // 28
               

	    // Write head/reference particle/tail field information.
		    << FDext[0](0) << setw(pwi) << "\t"                                         // 29 B-head x
		    << FDext[0](1) << setw(pwi) << "\t"                                         // 30 B-head y
		    << FDext[0](2) << setw(pwi) << "\t"                                         // 31 B-head z

		    << FDext[1](0) << setw(pwi) << "\t"                                         // 32 E-head x
		    << FDext[1](1) << setw(pwi) << "\t"                                         // 33 E-head y 
		    << FDext[1](2) << setw(pwi) << "\t"                                         // 34 E-head z

		    << FDext[2](0) << setw(pwi) << "\t"                                         // 35 B-ref x
		    << FDext[2](1) << setw(pwi) << "\t"                                         // 36 B-ref y 
		    << FDext[2](2) << setw(pwi) << "\t"                                         // 37 B-ref z

		    << FDext[3](0) << setw(pwi) << "\t"                                         // 38 E-ref x
		    << FDext[3](1) << setw(pwi) << "\t"                                         // 39 E-ref y
		    << FDext[3](2) << setw(pwi) << "\t"                                         // 40 E-ref z

		    << FDext[4](0) << setw(pwi) << "\t"                                         // 41 B-tail x
		    << FDext[4](1) << setw(pwi) << "\t"                                         // 42 B-tail y
		    << FDext[4](2) << setw(pwi) << "\t"                                         // 43 B-tail z

		    << FDext[5](0) << setw(pwi) << "\t"                                         // 44 E-tail x
		    << FDext[5](1) << setw(pwi) << "\t"                                         // 45 E-tail y 
		    << FDext[5](2) << setw(pwi) << "\t"                                         // 46 E-tail z
		    << endl;

        for (int p=0; p<Ippl::getNodes(); p++)
            os_lBalData << beam.getLoadBalance(p)  << setw(pwi) << "\t";
        os_lBalData << endl;
	
        os_statData.close();
        os_lBalData.close();
    }

    /// %Stop timer.
    IpplTimings::stopTimer(StatMarkerTimer_m);
}

void DataSink::writeSDDSHeader(ofstream &outputFile)
{
    outputFile << "SDDS1" << endl;
    outputFile << "&description text=\"Statistics data " << OPAL.getInputFn() << "\" " << endl;
    outputFile << ", contents=\"stat parameters\" &end" << endl;

    outputFile << "&parameter name=processors, type=long, ";
    outputFile << "description=\"Number of Processors\" &end" << endl;

    outputFile << "&column name=t, type=double, units=s, ";
    outputFile << "description=\"1 Time\" &end" << endl;
    outputFile << "&column name=s, type=double, units=m, ";
    outputFile << "description=\"2 Average Longitudinal Position\" &end" << endl;

    outputFile << "&column name=numParticles, type=long, units=1 ";
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

    outputFile << "&column name=beta_x, type=double, units=m , ";
    outputFile << "description=\"21 Beta function in x  \" &end" << endl;
    outputFile << "&column name=beta_y, type=double, units=m , ";
    outputFile << "description=\"22 Beta function in y  \" &end" << endl;

    outputFile << "&column name=alpha_x, type=double, units=1 , ";
    outputFile << "description=\"23 Alpha function in x  \" &end" << endl;
    outputFile << "&column name=alpha_y, type=double, units=1 , ";
    outputFile << "description=\"24 Alpha function in y  \" &end" << endl;

    outputFile << "&column name=Dx, type=double, units=m , ";
    outputFile << "description=\"25 Dispersion in x  \" &end" << endl;
    outputFile << "&column name=DDx, type=double, units=1 , ";
    outputFile << "description=\"26 Derivative of dispersion in x  \" &end" << endl;

    outputFile << "&column name=Dy, type=double, units=m , ";
    outputFile << "description=\"27 Dispersion in y  \" &end" << endl;
    outputFile << "&column name=DDy, type=double, units=1 , ";
    outputFile << "description=\"28 Derivative of dispersion in y  \" &end" << endl;

    outputFile << "&column name=Bx_head, type=double, units=T , ";
    outputFile << "description=\"29 Bx-Field component of head particle  \" &end" << endl;
    outputFile << "&column name=By_head, type=double, units=T , ";
    outputFile << "description=\"30 By-Field component of head particle  \" &end" << endl;
    outputFile << "&column name=Bz_head, type=double, units=T , ";
    outputFile << "description=\"31 Bz-Field component of head particle  \" &end" << endl;

    outputFile << "&column name=Ex_head, type=double, units=MV/m , ";
    outputFile << "description=\"32 Ex-Field component of head particle  \" &end" << endl;
    outputFile << "&column name=Ey_head, type=double, units=MV/m , ";
    outputFile << "description=\"33 Ey-Field component of head particle  \" &end" << endl;
    outputFile << "&column name=Ez_head, type=double, units=MV/m , ";
    outputFile << "description=\"34 Ez-Field component of head particle  \" &end" << endl;

    outputFile << "&column name=Bx_ref, type=double, units=T , ";
    outputFile << "description=\"35 Bx-Field component of ref particle  \" &end" << endl;
    outputFile << "&column name=By_ref, type=double, units=T , ";
    outputFile << "description=\"36 By-Field component of ref particle  \" &end" << endl;
    outputFile << "&column name=Bz_ref, type=double, units=T , ";
    outputFile << "description=\"37 Bz-Field component of ref particle  \" &end" << endl;

    outputFile << "&column name=Ex_ref, type=double, units=MV/m , ";
    outputFile << "description=\"38 Ex-Field component of ref particle  \" &end" << endl;
    outputFile << "&column name=Ey_ref, type=double, units=MV/m , ";
    outputFile << "description=\"39 Ey-Field component of ref particle  \" &end" << endl;
    outputFile << "&column name=Ez_ref, type=double, units=MV/m , ";
    outputFile << "description=\"40 Ez-Field component of ref particle  \" &end" << endl;

    outputFile << "&column name=Bx_tail, type=double, units=T , ";
    outputFile << "description=\"41 Bx-Field component of tail particle  \" &end" << endl;
    outputFile << "&column name=By_tail, type=double, units=T , ";
    outputFile << "description=\"42 By-Field component of tail particle  \" &end" << endl;
    outputFile << "&column name=Bz_tail, type=double, units=T , ";
    outputFile << "description=\"43 Bz-Field component of tail particle  \" &end" << endl;

    outputFile << "&column name=Ex_tail, type=double, units=MV/m , ";
    outputFile << "description=\"44 Ex-Field component of tail particle  \" &end" << endl;
    outputFile << "&column name=Ey_tail, type=double, units=MV/m , ";
    outputFile << "description=\"45 Ey-Field component of tail particle  \" &end" << endl;
    outputFile << "&column name=Ez_tail, type=double, units=MV/m , ";
    outputFile << "description=\"46 Ez-Field component of tail particle  \" &end" << endl;

    outputFile << "&data mode=ascii &end" << endl;

    outputFile << Ippl::getNodes() << endl;
}
/***************************************************************************
 * $RCSfile: DataSink.cpp,v $   $Author: adelmann $
 * $Revision: 1.3 $   $Date: 2004/06/02 19:38:54 $
 ***************************************************************************/

