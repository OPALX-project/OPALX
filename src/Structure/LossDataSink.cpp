#include <Structure/LossDataSink.h>

#include "Ippl.h"
#include "Algorithms/PBunchDefs.h"
#include "Utilities/Options.h"

using namespace Options;

LossDataSink::LossDataSink(std::string elem, bool hdf5Save):
    element_m(elem),
    h5hut_mode_m(hdf5Save)
{
    x_m.clear();
    y_m.clear();
    z_m.clear();
    px_m.clear();
    py_m.clear();
    pz_m.clear();
    id_m.clear();
    turn_m.clear();
    time_m.clear();
}

LossDataSink::LossDataSink(const LossDataSink &rsh):
    element_m(rsh.element_m),
    h5hut_mode_m(rsh.h5hut_mode_m)
{
    x_m.clear();
    y_m.clear();
    z_m.clear();
    px_m.clear();
    py_m.clear();
    pz_m.clear();
    id_m.clear();
    turn_m.clear();
    time_m.clear();
}

LossDataSink::LossDataSink() {
    LossDataSink(std::string("NULL"), false);
}

LossDataSink::~LossDataSink() {
}

void LossDataSink::openH5() {

    /// Open H5 file. Check that it opens correctly.
#ifdef PARALLEL_IO
    H5file_m = H5OpenFile(fn_m.c_str(), H5_O_WRONLY, Ippl::getComm());
#else
    H5file_m = H5OpenFile(fn_m.c_str(), H5_O_WRONLY, 0);
#endif

    if(H5file_m == (void*)H5_ERR) {
        ERRORMSG("h5 file open failed: exiting!" << endl);
        exit(0);
    }
}

void LossDataSink::writeHeaderH5() {

    h5_int64_t rc;
    // Write file attributes to describe phase space to H5 file.
    std::stringstream OPAL_version;
    OPAL_version << PACKAGE_NAME << " " << PACKAGE_VERSION << " svn rev. " << SVN_VERSION;
    rc = H5WriteFileAttribString(H5file_m, "OPAL_version", OPAL_version.str().c_str());
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteFileAttribString(H5file_m, "tUnit", "s");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteFileAttribString(H5file_m, "xUnit", "mm");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteFileAttribString(H5file_m, "yUnit", "mm");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteFileAttribString(H5file_m, "zUnit", "mm");
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
    rc = H5WriteFileAttribString(H5file_m, "idUnit", "1");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    /// in case of circular machines
    if (time_m.size() != 0) {
        rc = H5WriteFileAttribString(H5file_m, "turnUnit", "1");
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

        rc = H5WriteFileAttribString(H5file_m, "timeUnit", "s");
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    }

    rc = H5WriteFileAttribString(H5file_m, "SPOSUnit", "mm");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteFileAttribString(H5file_m, "TIMEUnit", "s");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5WriteFileAttribString(H5file_m, "mpart", "GeV");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5WriteFileAttribString(H5file_m, "qi", "C");
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    rc = H5AddAttachment(H5file_m, OpalData::getInstance()->getInputFn().c_str());
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
}

void LossDataSink::addParticle(const Vector_t x,const  Vector_t p, const size_t  id) {
    x_m.push_back(x(0));
    y_m.push_back(x(1));
    z_m.push_back(x(2));
    px_m.push_back(p(0));
    py_m.push_back(p(1));
    pz_m.push_back(p(2));
    id_m.push_back(id);
}

// For ring type simulation, dump the time and turn number
void LossDataSink::addParticle(const Vector_t x, const Vector_t p, const size_t id, const double time, const size_t turn) {
    addParticle(x, p, id);
    turn_m.push_back(turn);
    time_m.push_back(time);
}

void LossDataSink::save() {

    if (element_m == std::string("")) return;

    if (hasNoParticlesToDump()) return;

    if (h5hut_mode_m) {

        if (!Options::enableHDF5) return;

        fn_m = element_m + std::string(".h5");
	INFOMSG("Save " << fn_m << endl);
	openH5();
	writeHeaderH5();
	saveH5();
	H5CloseFile(H5file_m);
	Ippl::Comm->barrier();

    }
    else {

        fn_m = element_m + std::string(".loss");
        INFOMSG("Save " << fn_m << endl);
        if(OpalData::getInstance()->inRestartRun())
	    append();
        else
	    open();
        writeHeaderASCII();
        saveASCII();
        close();
        Ippl::Comm->barrier();

    }

    x_m.clear();
    y_m.clear();
    z_m.clear();
    px_m.clear();
    py_m.clear();
    pz_m.clear();
    id_m.clear();
    turn_m.clear();
    time_m.clear();
}

// Note: This was changed to calculate the global number of dumped particles
// because there are two cases to be considered:
// 1. ALL nodes have 0 lost particles -> nothing to be done.
// 2. Some nodes have 0 lost particles, some not -> H5 can handle that but all
// nodes HAVE to participate, otherwise H5 waits endlessly for a response from
// the nodes that didn't enter the saveH5 function. -DW
bool LossDataSink::hasNoParticlesToDump() {

    size_t nLoc = x_m.size();

    reduce(nLoc, nLoc, OpAddAssign());

    return nLoc == 0;
}

void LossDataSink::saveH5() {
    h5_int64_t rc;

    size_t nLoc = x_m.size();

    INFOMSG("saveH5 n = " << nLoc << endl);

    std::unique_ptr<char[]> varray(new char[(nLoc)*sizeof(double)]);
    double *farray = reinterpret_cast<double *>(varray.get());
    h5_int64_t *larray = reinterpret_cast<h5_int64_t *>(varray.get());

    std::unique_ptr<size_t[]> locN(new size_t[Ippl::getNodes()]);
    std::unique_ptr<size_t[]> globN(new size_t[Ippl::getNodes()]);

    for(int i = 0; i < Ippl::getNodes(); i++) {
        globN[i] = locN[i] = 0;
    }

    locN[Ippl::myNode()] = nLoc;
    reduce(locN.get(), locN.get() + Ippl::getNodes(), globN.get(), OpAddAssign());

    /// Set current record/time step.
    rc = H5SetStep(H5file_m, 0);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5PartSetNumParticles(H5file_m, nLoc);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    // Write all data
    for(size_t i = 0; i < nLoc; i++)
        farray[i] = x_m[i];
    rc = H5PartWriteDataFloat64(H5file_m, "x", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    for(size_t i = 0; i < nLoc; i++)
        farray[i] = y_m[i];
    rc = H5PartWriteDataFloat64(H5file_m, "y", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    for(size_t i = 0; i < nLoc; i++)
        farray[i] = z_m[i];
    rc = H5PartWriteDataFloat64(H5file_m, "z", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    for(size_t i = 0; i < nLoc; i++)
        farray[i] = px_m[i];
    rc = H5PartWriteDataFloat64(H5file_m, "px", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    for(size_t i = 0; i < nLoc; i++)
        farray[i] = py_m[i];
    rc = H5PartWriteDataFloat64(H5file_m, "py", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    for(size_t i = 0; i < nLoc; i++)
        farray[i] = pz_m[i];
    rc = H5PartWriteDataFloat64(H5file_m, "pz", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    /// Write particle id numbers.
    for(size_t i = 0; i < nLoc; i++)
        larray[i] =  id_m[i];
    rc = H5PartWriteDataInt64(H5file_m, "id", larray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    if (time_m.size() != 0) {
        for(size_t i = 0; i < nLoc; i++)
            farray[i] = time_m[i];
        rc = H5PartWriteDataFloat64(H5file_m, "time", farray);
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

        for(size_t i = 0; i < nLoc; i++)
            larray[i] =  turn_m[i];
        rc = H5PartWriteDataInt64(H5file_m, "turn", larray);
        if(rc != H5_SUCCESS)
            ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    }
    H5call_m++;
}

void LossDataSink::saveASCII() {

    /*
      ASCII output
    */
    int tag = Ippl::Comm->next_tag(IPPL_APP_TAG3, IPPL_APP_CYCLE);
    if(Ippl::Comm->myNode() == 0) {
        const unsigned partCount = x_m.size();

        if (time_m.size() != 0) {
            for(unsigned i = 0; i < partCount; i++) {
                os_m << element_m   << "   ";
                os_m << x_m[i] << "   ";
                os_m << y_m[i] << "   ";
                os_m << z_m[i] << "   ";
                os_m << px_m[i] << "   ";
                os_m << py_m[i] << "   ";
                os_m << pz_m[i] << "   ";
                os_m << id_m[i]   << "   ";
                os_m << turn_m[i] << "   ";
                os_m << time_m[i] << " " << std::endl;
            }
        }
        else {
            for(unsigned i = 0; i < partCount; i++) {
                os_m << element_m   << "   ";
                os_m << x_m[i] << "   ";
                os_m << y_m[i] << "   ";
                os_m << z_m[i] << "   ";
                os_m << px_m[i] << "   ";
                os_m << py_m[i] << "   ";
                os_m << pz_m[i] << "   ";
                os_m << id_m[i]   << "   " << std::endl;
            }
        }
        int notReceived =  Ippl::getNodes() - 1;
        while(notReceived > 0) {
            unsigned dataBlocks = 0;
            int node = COMM_ANY_NODE;
            Message *rmsg =  Ippl::Comm->receive_block(node, tag);
            if(rmsg == 0) {
                ERRORMSG("LossDataSink: Could not receive from client nodes output." << endl);
            }
            notReceived--;
            rmsg->get(&dataBlocks);
            if (time_m.size() != 0) {
                for(unsigned i = 0; i < dataBlocks; i++) {
                    long id, turn;
                    double rx, ry, rz, px, py, pz, time;
                    rmsg->get(&id);
                    rmsg->get(&rx);
                    rmsg->get(&ry);
                    rmsg->get(&rz);
                    rmsg->get(&px);
                    rmsg->get(&py);
                    rmsg->get(&pz);
                    rmsg->get(&turn);
                    rmsg->get(&time);
                    os_m << element_m << "   ";
                    os_m << rx << "   ";
                    os_m << ry << "   ";
                    os_m << rz << "   ";
                    os_m << px << "   ";
                    os_m << py << "   ";
                    os_m << pz << "   ";
                    os_m << id << "   ";
                    os_m << turn << "   ";
                    os_m << time << std::endl;
                }
            }
            else {
                for(unsigned i = 0; i < dataBlocks; i++) {
                    long id;
                    double rx, ry, rz, px, py, pz;
                    rmsg->get(&id);
                    rmsg->get(&rx);
                    rmsg->get(&ry);
                    rmsg->get(&rz);
                    rmsg->get(&px);
                    rmsg->get(&py);
                    rmsg->get(&pz);
                    os_m << element_m << "   ";
                    os_m << rx << "   ";
                    os_m << ry << "   ";
                    os_m << rz << "   ";
                    os_m << px << "   ";
                    os_m << py << "   ";
                    os_m << pz << "   ";
                    os_m << id << " " << std::endl;
                }
            }
            delete rmsg;
        }
    } else {
        Message *smsg = new Message();
        const unsigned msgsize = x_m.size();
        smsg->put(msgsize);
        if (time_m.size() != 0) {
            for(unsigned i = 0; i < msgsize; i++) {
                smsg->put(id_m[i]);
                smsg->put(x_m[i]);
                smsg->put(y_m[i]);
                smsg->put(z_m[i]);
                smsg->put(px_m[i]);
                smsg->put(py_m[i]);
                smsg->put(pz_m[i]);
                smsg->put(turn_m[i]);
                smsg->put(time_m[i]);
            }
        }
        else {
            for(unsigned i = 0; i < msgsize; i++) {
                smsg->put(id_m[i]);
                smsg->put(x_m[i]);
                smsg->put(y_m[i]);
                smsg->put(z_m[i]);
                smsg->put(px_m[i]);
                smsg->put(py_m[i]);
                smsg->put(pz_m[i]);
            }
        }
        bool res = Ippl::Comm->send(smsg, 0, tag);
        if(! res)
            ERRORMSG("LossDataSink Ippl::Comm->send(smsg, 0, tag) failed " << endl;);
    }
}
