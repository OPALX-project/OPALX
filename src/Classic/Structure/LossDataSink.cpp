#include "Structure/LossDataSink.h"

#include "Ippl.h"
#include "Utilities/Options.h"
#include "AbstractObjects/OpalData.h"
#include "Utilities/GeneralClassicException.h"
#include "Utilities/Util.h"

#include <boost/filesystem.hpp>
#include "OPALconfig.h"

#include <cassert>

#define REPORTONERROR(rc) reportOnError(rc, __FILE__, __LINE__)

extern Inform *gmsg;

std::map<double, std::string> LossDataSink::statFileEntries_s;

LossDataSink::LossDataSink(std::string elem, bool hdf5Save, ElementBase::ElementType type):
    element_m(elem),
    h5hut_mode_m(hdf5Save),
    H5file_m(0),
    H5call_m(0),
    type_m(type)
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
    h5hut_mode_m(rsh.h5hut_mode_m),
    H5file_m(rsh.H5file_m),
    H5call_m(rsh.H5call_m),
    RefPartR_m(rsh.RefPartR_m),
    RefPartP_m(rsh.RefPartP_m),
    globalTrackStep_m(rsh.globalTrackStep_m),
    refTime_m(rsh.refTime_m),
    spos_m(rsh.spos_m),
    type_m(rsh.type_m)
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
#define _unused(x) ((void)(x))

LossDataSink::~LossDataSink() {
    if (H5file_m) {
        h5_int64_t rc;
        rc = H5CloseFile(H5file_m);
        assert (rc != H5_ERR);
        _unused(rc);

        H5file_m = 0;
    }
    Ippl::Comm->barrier();
}

void LossDataSink::openH5(h5_int32_t mode) {
#if defined (USE_H5HUT2)
    h5_prop_t props = H5CreateFileProp ();
    MPI_Comm comm = Ippl::getComm();
    H5SetPropFileMPIOCollective (props, &comm);
    H5file_m = H5OpenFile (fn_m.c_str(), mode, props);

    if(H5file_m == (h5_file_t)H5_ERR) {
        throw GeneralClassicException("LossDataSink::openH5",
                                      "failed to open h5 file '" + fn_m + "'");
    }
#else

    H5file_m = H5OpenFile(fn_m.c_str(), mode, Ippl::getComm());
    if(H5file_m == (void*)H5_ERR) {
        throw GeneralClassicException("LossDataSink::openH5",
                                      "failed to open h5 file '" + fn_m + "'");
    }
#endif
}

void LossDataSink::writeHeaderH5() {

    // Write file attributes to describe phase space to H5 file.
    std::stringstream OPAL_version;
    OPAL_version << PACKAGE_NAME << " " << PACKAGE_VERSION_STR << " git rev. " << Util::getGitRevision();
    REPORTONERROR(H5WriteFileAttribString(H5file_m, "OPAL_version", OPAL_version.str().c_str()));
    REPORTONERROR(H5WriteFileAttribString(H5file_m, "tUnit", "s"));
    REPORTONERROR(H5WriteFileAttribString(H5file_m, "xUnit", "m"));
    REPORTONERROR(H5WriteFileAttribString(H5file_m, "yUnit", "m"));
    REPORTONERROR(H5WriteFileAttribString(H5file_m, "zUnit", "m"));
    REPORTONERROR(H5WriteFileAttribString(H5file_m, "pxUnit", "#beta#gamma"));
    REPORTONERROR(H5WriteFileAttribString(H5file_m, "pyUnit", "#beta#gamma"));
    REPORTONERROR(H5WriteFileAttribString(H5file_m, "pzUnit", "#beta#gamma"));
    REPORTONERROR(H5WriteFileAttribString(H5file_m, "idUnit", "1"));

    /// in case of circular machines
    if (hasTimeAttribute()) {
        REPORTONERROR(H5WriteFileAttribString(H5file_m, "turnUnit", "1"));

        REPORTONERROR(H5WriteFileAttribString(H5file_m, "timeUnit", "s"));
    }
    REPORTONERROR(H5WriteFileAttribString(H5file_m, "SPOSUnit", "mm"));
    REPORTONERROR(H5WriteFileAttribString(H5file_m, "TIMEUnit", "s"));

    REPORTONERROR(H5WriteFileAttribString(H5file_m, "mpart", "GeV"));
    REPORTONERROR(H5WriteFileAttribString(H5file_m, "qi", "C"));

    REPORTONERROR(H5AddAttachment(H5file_m, OpalData::getInstance()->getInputFn().c_str()));
}

void LossDataSink::addReferenceParticle(const Vector_t &x,
                                        const  Vector_t &p,
                                        double time,
                                        double spos,
                                        long long globalTrackStep
                                        ) {
    RefPartR_m.push_back(x);
    RefPartP_m.push_back(p);
    globalTrackStep_m.push_back((h5_int64_t)globalTrackStep);
    spos_m.push_back(spos);
    refTime_m.push_back(time);
}

void LossDataSink::addParticle(const Vector_t &x,const  Vector_t &p, const size_t  id) {
    x_m.push_back(x(0));
    y_m.push_back(x(1));
    z_m.push_back(x(2));
    px_m.push_back(p(0));
    py_m.push_back(p(1));
    pz_m.push_back(p(2));
    id_m.push_back(id);
}

// For ring type simulation, dump the time and turn number
void LossDataSink::addParticle(const Vector_t &x, const Vector_t &p, const size_t id, const double time, const size_t turn) {
    addParticle(x, p, id);
    turn_m.push_back(turn);
    time_m.push_back(time);
}

void LossDataSink::save(unsigned int numSets) {

    if (element_m == std::string("")) return;

    if (hasNoParticlesToDump()) return;

    if (h5hut_mode_m) {
        namespace fs = boost::filesystem;

        if (!Options::enableHDF5) return;

        fn_m = element_m + std::string(".h5");
	INFOMSG(level2 << "Save " << fn_m << endl);
        if (Options::openMode == Options::WRITE || !fs::exists(fn_m)) {
            openH5();
            writeHeaderH5();
        } else {
#ifdef H5HUT_API_VERSION
            openH5(H5_O_APPENDONLY);
#else
            openH5(H5_O_APPEND);
#endif
            H5call_m = H5GetNumSteps(H5file_m);
        }

        splitSets(numSets);
        saveStatistics(numSets);
        for (unsigned int i = 0; i < numSets; ++ i) {
            saveH5(i);
        }
	H5CloseFile(H5file_m);
	H5file_m = 0;
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
    spos_m.clear();
    refTime_m.clear();
    RefPartR_m.clear();
    RefPartP_m.clear();
    globalTrackStep_m.clear();
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

bool LossDataSink::hasTimeAttribute() {

    size_t tLoc = time_m.size();

    reduce(tLoc, tLoc, OpAddAssign());

    return tLoc > 0;
}





void LossDataSink::saveH5(unsigned int setIdx) {
    size_t startIdx = 0;
    size_t nLoc = x_m.size();
    if (setIdx + 1 < startSet_m.size()) {
        startIdx = startSet_m[setIdx];
        nLoc = startSet_m[setIdx + 1] - startSet_m[setIdx];
    }

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
    REPORTONERROR(H5SetStep(H5file_m, H5call_m));
    REPORTONERROR(H5PartSetNumParticles(H5file_m, nLoc));

    if (setIdx < spos_m.size()) {
        REPORTONERROR(H5WriteStepAttribFloat64(H5file_m, "SPOS", &(spos_m[setIdx]), 1));
        REPORTONERROR(H5WriteStepAttribFloat64(H5file_m, "TIME", &(refTime_m[setIdx]), 1));
        REPORTONERROR(H5WriteStepAttribFloat64(H5file_m, "RefPartR", (h5_float64_t *)&(RefPartR_m[setIdx]), 3));
        REPORTONERROR(H5WriteStepAttribFloat64(H5file_m, "RefPartP", (h5_float64_t *)&(RefPartP_m[setIdx]), 3));
        REPORTONERROR(H5WriteStepAttribInt64(H5file_m, "GlobalTrackStep", &(globalTrackStep_m[setIdx]), 1));
    }
    // Write all data
    size_t j = startIdx;
    for(size_t i = 0; i < nLoc; ++ i, ++ j)
        farray[i] = x_m[j];
    REPORTONERROR(H5PartWriteDataFloat64(H5file_m, "x", farray));

    j = startIdx;
    for(size_t i = 0; i < nLoc; ++ i, ++ j)
        farray[i] = y_m[j];
    REPORTONERROR(H5PartWriteDataFloat64(H5file_m, "y", farray));

    j = startIdx;
    for(size_t i = 0; i < nLoc; ++ i, ++ j)
        farray[i] = z_m[j];
    REPORTONERROR(H5PartWriteDataFloat64(H5file_m, "z", farray));

    j = startIdx;
    for(size_t i = 0; i < nLoc; ++ i, ++ j)
        farray[i] = px_m[j];
    REPORTONERROR(H5PartWriteDataFloat64(H5file_m, "px", farray));

    j = startIdx;
    for(size_t i = 0; i < nLoc; ++ i, ++ j)
        farray[i] = py_m[j];
    REPORTONERROR(H5PartWriteDataFloat64(H5file_m, "py", farray));

    j = startIdx;
    for(size_t i = 0; i < nLoc; ++ i, ++ j)
        farray[i] = pz_m[j];
    REPORTONERROR(H5PartWriteDataFloat64(H5file_m, "pz", farray));

    /// Write particle id numbers.
    j = startIdx;
    for(size_t i = 0; i < nLoc; ++ i, ++ j)
        larray[i] =  id_m[j];
    REPORTONERROR(H5PartWriteDataInt64(H5file_m, "id", larray));

    if (hasTimeAttribute()) {
        j = startIdx;
        for(size_t i = 0; i < nLoc; ++ i, ++ j)
            farray[i] = time_m[j];
        REPORTONERROR(H5PartWriteDataFloat64(H5file_m, "time", farray));

        j = startIdx;
        for(size_t i = 0; i < nLoc; ++ i, ++ j)
            larray[i] =  turn_m[j];
        REPORTONERROR(H5PartWriteDataInt64(H5file_m, "turn", larray));
    }

    ++ H5call_m;
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

void LossDataSink::splitSets(unsigned int numSets) {
    if (numSets <= 1) return;

    const size_t nLoc = x_m.size();
    size_t avgNumPerSet = nLoc / numSets;
    std::vector<size_t> numPartsInSet(numSets, avgNumPerSet);
    size_t test = numSets * avgNumPerSet;
    for (unsigned int j = 0; j < (nLoc - numSets * avgNumPerSet); ++ j) {
        ++ numPartsInSet[j];
        ++ test;
    }

    if (/*nLoc > 0 && */time_m.size() == nLoc) {
        std::vector<double> meanT(2 * numSets, 0.0);
        std::vector<double> timeRange(2 * numSets, 0.0);
        double minmaxT[] = {time_m[0], -time_m[0]};
        size_t i = 0;
        for (unsigned int j = 0; j < numSets; ++ j) {
            const size_t &numThisSet = numPartsInSet[j];
            for (size_t k = 0; k < numThisSet; ++ k, ++ i) {
                meanT[2 * j] += time_m[i];
                minmaxT[0] = std::min(minmaxT[0], time_m[i]);
                minmaxT[1] = std::min(minmaxT[1], -time_m[i]);
            }
            meanT[2 * j + 1] = numThisSet;
        }

        reduce(&(meanT[0]), &(meanT[0]) + 2 * numSets, &(meanT[0]), OpAddAssign());
        reduce(minmaxT, minmaxT + 2, minmaxT, OpMinAssign());

        for (unsigned int j = 0; j < numSets; ++ j) {
            meanT[2 * j] /= meanT[2 * j + 1];
        }

        timeRange[0] = minmaxT[0];
        for (unsigned int j = 1; j < numSets; ++ j) {
            timeRange[2 * j - 1] = 0.5 * (meanT[2 * (j - 1)] + meanT[2 * j]);
            timeRange[2 * j] = timeRange[2 * j - 1];
        }
        timeRange[2 * numSets - 1] = -minmaxT[1];

        size_t j = 0;
        size_t idxPrior = 0;
        startSet_m.push_back(0);

        for (size_t idx = 0; idx < nLoc; ++ idx) {
            if (time_m[idx] > timeRange[2 * j + 1]) {
                startSet_m.push_back(idx);
                numPartsInSet[j] = idx - idxPrior;
                idxPrior = idx;
                ++ j;
            }
        }
        numPartsInSet[numSets - 1] = nLoc - idxPrior;
        startSet_m.push_back(nLoc);

        i = 0;
        for (unsigned int j = 0; j < numSets; ++ j) {
            const size_t &numThisSet = numPartsInSet[j];
            for (size_t k = 0; k < numThisSet; ++ k, ++ i) {
                meanT[2 * j] += time_m[i];
                minmaxT[0] = std::min(minmaxT[0], time_m[i]);
                minmaxT[1] = std::min(minmaxT[1], -time_m[i]);
            }
            meanT[2 * j + 1] = numThisSet;
        }

        reduce(&(meanT[0]), &(meanT[0]) + 2 * numSets, &(meanT[0]), OpAddAssign());

        for (unsigned int j = 0; j < numSets; ++ j) {
            meanT[2 * j] /= meanT[2 * j + 1];
        }

    }
}

void LossDataSink::saveStatistics(unsigned int numSets) {
    if (type_m != ElementBase::MONITOR || !hasTimeAttribute()) return;

    for (unsigned int setIdx = 0; setIdx < numSets; ++ setIdx) {

        size_t startIdx = 0;
        size_t nLoc = x_m.size();
        if (setIdx + 1 < startSet_m.size()) {
            startIdx = startSet_m[setIdx];
            nLoc = startSet_m[setIdx + 1] - startSet_m[setIdx];
        }

        double tmean(0.0), trms(0.0);
        Vector_t rmean(0.0), pmean(0.0), rrms(0.0), prms(0.0), rprms(0.0), normEmit(0.0);
        Vector_t rsqsum(0.0), psqsum(0.0), rpsum(0.0), eps2(0.0), eps_norm(0.0), fac(0.0);

        double part[6];

        const unsigned int totalSize = 45;
        double plainData[totalSize];
        double rmax[] = {0.0, 0.0, 0.0};

        {
            Util::KahanAccumulation data[totalSize];
            Util::KahanAccumulation *centroid = data + 1;
            Util::KahanAccumulation *moments = data + 7;
            Util::KahanAccumulation *others = data + 43;

            data[0].sum = nLoc;

            unsigned int idx = startIdx;
            for(unsigned long k = 0; k < nLoc; ++ k, ++ idx) {
                part[1] = px_m[idx];
                part[3] = py_m[idx];
                part[5] = pz_m[idx];
                part[0] = x_m[idx];
                part[2] = y_m[idx];
                part[4] = z_m[idx];

                for(int i = 0; i < 6; i++) {
                    centroid[i] += part[i];
                    for(int j = 0; j <= i; j++) {
                        moments[i * 6 + j] += part[i] * part[j];
                    }
                }
                others[0] += time_m[idx];
                others[1] += std::pow(time_m[idx], 2);

                rmax[0] = std::max(rmax[0], std::abs(x_m[idx]));
                rmax[1] = std::max(rmax[1], std::abs(y_m[idx]));
                rmax[2] = std::max(rmax[2], std::abs(z_m[idx]));
            }

            for(int i = 0; i < 6; i++) {
                for(int j = 0; j < i; j++) {
                    moments[j * 6 + i] = moments[i * 6 + j];
                }
            }

            for (unsigned int i = 0; i < totalSize; ++ i) {
                plainData[i] = data[i].sum;
            }
        }

        reduce(plainData, plainData + totalSize, plainData, OpAddAssign());
        reduce(rmax, rmax + 3, rmax, OpMaxAssign());

        if (Ippl::myNode() != 0) continue;

        double *centroid = plainData + 1;
        double *moments = plainData + 7;
        double *others = plainData + 43;
        double nTotal = plainData[0];
        for(unsigned int i = 0 ; i < 3u; i++) {
            rmean(i) = centroid[2 * i] / nTotal;
            pmean(i) = centroid[(2 * i) + 1] / nTotal;
            rsqsum(i) = moments[2 * i * 6 + 2 * i] - nTotal * rmean(i) * rmean(i);
            psqsum(i) = std::max(0.0, moments[(2 * i + 1) * 6 + (2 * i) + 1] - nTotal * pmean(i) * pmean(i));
            rpsum(i) = moments[(2 * i) * 6 + (2 * i) + 1] - nTotal * rmean(i) * pmean(i);
        }
        tmean = others[0] / nTotal;
        trms = sqrt(std::max(0.0, (others[1] / nTotal - std::pow(tmean, 2))));

        eps2 = (rsqsum * psqsum - rpsum * rpsum) / (1.0 * nTotal * nTotal);

        rpsum /= nTotal;

        for(unsigned int i = 0 ; i < 3u; i++) {
            rrms(i) = sqrt(std::max(0.0, rsqsum(i)) / nTotal);
            prms(i) = sqrt(std::max(0.0, psqsum(i)) / nTotal);
            eps_norm(i)  =  std::sqrt(std::max(0.0, eps2(i)));
            double tmp = rrms(i) * prms(i);
            fac(i) = (tmp == 0) ? 0.0 : 1.0 / tmp;
        }

        rprms = rpsum * fac;

        std::stringstream statOut;
        statOut.precision(8);
        statOut << element_m << "\t"
                << spos_m[setIdx] << "\t"
                << refTime_m[setIdx] * 1e9 << "\t"
                << (unsigned int)floor(nTotal + 0.5) << "\t"
                << rrms(0) << "\t"
                << rrms(1) << "\t"
                << rrms(2) << "\t"
                << trms * 1e9 << "\t"
                << prms(0) << "\t"
                << prms(1) << "\t"
                << prms(2) << "\t"
                << eps_norm(0) << "\t"
                << eps_norm(1) << "\t"
                << eps_norm(2) << "\t"
                << rmean(0) << "\t"
                << rmean(1) << "\t"
                << rmean(2) << "\t"
                << tmean * 1e9 << "\t"
                << RefPartR_m[setIdx](0) << "\t"
                << RefPartR_m[setIdx](1) << "\t"
                << RefPartR_m[setIdx](2) << "\t"
                << RefPartP_m[setIdx](0) << "\t"
                << RefPartP_m[setIdx](1) << "\t"
                << RefPartP_m[setIdx](2) << "\t"
                << rmax[0] << "\t"
                << rmax[1] << "\t"
                << rmax[2] << "\t"
                << rprms(0) << "\t"
                << rprms(1) << "\t"
                << rprms(2) << "\t"
                << std::endl;

        statFileEntries_s.insert(std::make_pair(spos_m[setIdx], statOut.str()));
    }
}

void LossDataSink::writeStatistics() {
    if (Ippl::myNode() != 0 || statFileEntries_s.size() == 0) return;

    namespace fs = boost::filesystem;

    std::string fileName = OpalData::getInstance()->getInputBasename() + std::string("_Monitors.stat");
    std::ofstream statOut;

    if (Options::openMode == Options::APPEND && fs::exists(fileName)) {
        Util::rewindLinesSDDS(fileName, statFileEntries_s.begin()->first);

        statOut.open(fileName, std::ios::app);
    } else {
        statOut.open(fileName);
        std::string indent("        ");

        statOut << "SDDS1\n";
        statOut << "&description \n"
                << indent << "text=\"Statistics data of monitors\", \n"
                << indent << "contents=\"stat parameters\"\n"
                << "&end\n";
        statOut << "&parameter \n"
                << indent << "name=revision, \n"
                << indent << "type=string, \n"
                << indent << "description=\"git revision of opal\"\n"
                << "&end\n";
        statOut << "&column \n"
                << indent << "name=name, \n"
                << indent << "type=string \n"
                << indent << "description=\"1 Monitor name\"\n"
                << "&end\n";
        statOut << "&column \n"
                << indent << "name=s, \n"
                << indent << "type=double, \n"
                << indent << "units=m, \n"
                << indent << "description=\"2 Longitudinal Position\"\n"
                << "&end\n";
        statOut << "&column \n"
                << indent << "name=t, \n"
                << indent << "type=double, \n"
                << indent << "units=ns, \n"
                << indent << "description=\"3 Passage Time Reference Particle\"\n"
                << "&end\n";
        statOut << "&column \n"
                << indent << "name=numParticles, \n"
                << indent << "type=long, \n"
                << indent << "units=1, \n"
                << indent << "description=\"4 Number of Macro Particles\"\n"
                << "&end\n";
        statOut << "&column \n"
                << indent << "name=rms_x, \n"
                << indent << "type=double, \n"
                << indent << "units=m, \n"
                << indent << "description=\"5 RMS Beamsize in x\"\n"
                << "&end\n";
        statOut << "&column \n"
                << indent << "name=rms_y, \n"
                << indent << "type=double, \n"
                << indent << "units=m, \n"
                << indent << "description=\"6 RMS Beamsize in y\"\n"
                << "&end\n";
        statOut << "&column \n"
                << indent << "name=rms_s, \n"
                << indent << "type=double, \n"
                << indent << "units=m, \n"
                << indent << "description=\"7 RMS Beamsize in s\"\n"
                << "&end\n";
        statOut << "&column \n"
                << indent << "name=rms_t, \n"
                << indent << "type=double, \n"
                << indent << "units=ns, \n"
                << indent << "description=\"8 RMS Passage Time\"\n"
                << "&end\n";
        statOut << "&column \n"
                << indent << "name=rms_px, \n"
                << indent << "type=double, \n"
                << indent << "units=1, \n"
                << indent << "description=\"9 RMS Momenta in x\"\n"
                << "&end\n";
        statOut << "&column \n"
                << indent << "name=rms_py, \n"
                << indent << "type=double, \n"
                << indent << "units=1, \n"
                << indent << "description=\"10 RMS Momenta in y\"\n"
                << "&end\n";
        statOut << "&column \n"
                << indent << "name=rms_ps, \n"
                << indent << "type=double, \n"
                << indent << "units=1, \n"
                << indent << "description=\"11 RMS Momenta in s\"\n"
                << "&end\n";
        statOut << "&column \n"
                << indent << "name=emit_x, \n"
                << indent << "type=double, \n"
                << indent << "units=m, \n"
                << indent << "description=\"12 Normalized Emittance x\"\n"
                << "&end\n";
        statOut << "&column \n"
                << indent << "name=emit_y, \n"
                << indent << "type=double, \n"
                << indent << "units=m, \n"
                << indent << "description=\"13 Normalized Emittance y\"\n"
                << "&end\n";
        statOut << "&column \n"
                << indent << "name=emit_s, \n"
                << indent << "type=double, \n"
                << indent << "units=m, \n"
                << indent << "description=\"14 Normalized Emittance s\"\n"
                << "&end\n";
        statOut << "&column \n"
                << indent << "name=mean_x, \n"
                << indent << "type=double, \n"
                << indent << "units=m, \n"
                << indent << "description=\"15 Mean Beam Position in x\"\n"
                << "&end\n";
        statOut << "&column \n"
                << indent << "name=mean_y, \n"
                << indent << "type=double, \n"
                << indent << "units=m, \n"
                << indent << "description=\"16 Mean Beam Position in y\"\n"
                << "&end\n";
        statOut << "&column \n"
                << indent << "name=mean_s, \n"
                << indent << "type=double, \n"
                << indent << "units=m, \n"
                << indent << "description=\"17 Mean Beam Position in s\"\n"
                << "&end\n";
        statOut << "&column \n"
                << indent << "name=mean_t, \n"
                << indent << "type=double, \n"
                << indent << "units=ns, \n"
                << indent << "description=\"18 Mean Passage Time\"\n"
                << "&end\n";
        statOut << "&column \n"
                << indent << "name=ref_x, \n"
                << indent << "type=double, \n"
                << indent << "units=m, \n"
                << indent << "description=\"19 x coordinate of reference particle in lab cs\"\n"
                << "&end\n";
        statOut << "&column \n"
                << indent << "name=ref_y, \n"
                << indent << "type=double, \n"
                << indent << "units=m, \n"
                << indent << "description=\"20 y coordinate of reference particle in lab cs\"\n"
                << "&end\n";
        statOut << "&column \n"
                << indent << "name=ref_z, \n"
                << indent << "type=double, \n"
                << indent << "units=m, \n"
                << indent << "description=\"21 z coordinate of reference particle in lab cs\"\n"
                << "&end\n";
        statOut << "&column \n"
                << indent << "name=ref_px, \n"
                << indent << "type=double, \n"
                << indent << "units=1, \n"
                << indent << "description=\"22 x momentum of reference particle in lab cs\"\n"
                << "&end\n";
        statOut << "&column \n"
                << indent << "name=ref_py, \n"
                << indent << "type=double, \n"
                << indent << "units=1, \n"
                << indent << "description=\"23 y momentum of reference particle in lab cs\"\n"
                << "&end\n";
        statOut << "&column \n"
                << indent << "name=ref_pz, \n"
                << indent << "type=double, \n"
                << indent << "units=1, \n"
                << indent << "description=\"24 z momentum of reference particle in lab cs\"\n"
                << "&end\n";
        statOut << "&column \n"
                << indent << "name=max_x, \n"
                << indent << "type=double, \n"
                << indent << "units=m, \n"
                << indent << "description=\"25 Max Beamsize in x\"\n"
                << "&end\n";
        statOut << "&column \n"
                << indent << "name=max_y, \n"
                << indent << "type=double, \n"
                << indent << "units=m, \n"
                << indent << "description=\"26 Max Beamsize in y\"\n"
                << "&end\n";
        statOut << "&column \n"
                << indent << "name=max_s, \n"
                << indent << "type=double, \n"
                << indent << "units=m, \n"
                << indent << "description=\"27 Max Beamsize in s\"\n"
                << "&end\n";
        statOut << "&column \n"
                << indent << "name=xpx, \n"
                << indent << "type=double, \n"
                << indent << "units=1, \n"
                << indent << "description=\"28 Correlation xpx\"\n"
                << "&end\n";
        statOut << "&column \n"
                << indent << "name=ypy, \n"
                << indent << "type=double, \n"
                << indent << "units=1, \n"
                << indent << "description=\"29 Correlation ypyy\"\n"
                << "&end\n";
        statOut << "&column \n"
                << indent << "name=zpz, \n"
                << indent << "type=double, \n"
                << indent << "units=1, \n"
                << indent << "description=\"30 Correlation zpz\"\n"
                << "&end\n";
        statOut << "&data \n"
                << indent << "mode=ascii,\n"
                << indent << "no_row_counts=1\n"
                << "&end\n";

        statOut << PACKAGE_NAME << " " << PACKAGE_VERSION << " git rev. " << Util::getGitRevision() << std::endl;
    }

    for (auto entry: statFileEntries_s) {
        statOut << entry.second;
    }

    statOut.close();
    statFileEntries_s.clear();
}