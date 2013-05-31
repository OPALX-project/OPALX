#include <fstream>
#include <ios>

#include "Fields/FM3DH5Block_nonescale.hh"
#include "Fields/Fieldmap.icc"
#include "H5hut.h"
#include "Physics/Physics.h"

extern Inform *gmsg;

using namespace std;
using Physics::mu_0;

FM3DH5Block_nonescale::FM3DH5Block_nonescale(string aFilename):
    Fieldmap(aFilename) {
    Inform msg("FM3DH5 ");
    h5_err_t h5err;
    h5_size_t grid_rank;
    h5_size_t grid_dims[3];
    h5_size_t field_dims;
    char name[20];
    h5_size_t len_name = 20;

    Type = T3DDynamicH5Block;
    h5_file_t *file = H5OpenFile(aFilename.c_str(), H5_O_RDONLY, Ippl::getComm());


    if(file) {
        h5_int64_t last_step = H5GetNumSteps(file) - 1;
        h5err = H5SetStep(file, last_step);
        if(h5err != H5_SUCCESS)
            ERRORMSG("H5 rc= " << h5err << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        h5_int64_t num_fields = H5BlockGetNumFields(file);


        for(h5_ssize_t i = 0; i < num_fields; ++ i) {
            /*
              Work around API changes in H5Block
            */
            h5_int64_t ftype;
            h5err = H5BlockGetFieldInfo(file, (h5_size_t)i, name, len_name, &grid_rank, grid_dims, &field_dims, &ftype);
            if(h5err != H5_SUCCESS)
                ERRORMSG("H5 rc= " << h5err << " in " << __FILE__ << " @ line " << __LINE__ << endl);
            if(strcmp(name, "Efield") == 0) {
                num_gridpx_m = grid_dims[0];
                num_gridpy_m = grid_dims[1];
                num_gridpz_m = grid_dims[2];
            }
        }
        h5err = H5Block3dGetFieldSpacing(file, "Efield", &hx_m, &hy_m, &hz_m);
        if(h5err != H5_SUCCESS)
            ERRORMSG("H5 rc= " << h5err << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        h5err = H5Block3dGetFieldOrigin(file, "Efield", &xbegin_m, &ybegin_m, &zbegin_m);
        if(h5err != H5_SUCCESS)
            ERRORMSG("H5 rc= " << h5err << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        xend_m = xbegin_m + (num_gridpx_m - 1) * hx_m;
        yend_m = ybegin_m + (num_gridpy_m - 1) * hy_m;
        zend_m = zbegin_m + (num_gridpz_m - 1) * hz_m;

        //         xcentral_idx_m = static_cast<int>(fabs(xbegin_m) / hx_m);
        //         ycentral_idx_m = static_cast<int>(fabs(ybegin_m) / hy_m);


        h5err = H5ReadFileAttribFloat64(file, "Resonance Frequency(Hz)", &frequency_m);
        if(h5err != H5_SUCCESS)
            ERRORMSG("H5 rc= " << h5err << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        frequency_m *= Physics::two_pi;

        h5err = H5CloseFile(file);
        if(h5err != H5_SUCCESS)
            ERRORMSG("H5 rc= " << h5err << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    }
}


FM3DH5Block_nonescale::~FM3DH5Block_nonescale() {

}

void FM3DH5Block_nonescale::readMap() {
    if(FieldstrengthEz_m.empty()) {
        Inform msg("FM3DH5_NS ");
        h5_file_t *file = H5OpenFile(Filename_m.c_str(), H5_O_RDONLY, Ippl::getComm());

        if(file) {
            h5_int64_t h5err;

            long field_size = 0;
            int Nnodes = Ippl::getNodes();//min(20, Ippl::getNodes());
            int Nz_avrg = static_cast<int>(floor(0.5 + num_gridpz_m / Nnodes));
            int Nz_diff = Nz_avrg * Nnodes - num_gridpz_m;
            int signNz = Nz_diff > 0 ? 1 : -1;
            int *Nz_read_start = new int[Ippl::getNodes() + 1];
            int *Nz_read_length = new int[Ippl::getNodes()];
            int N_read_start;
            int start = 0;
            // int rbuf_size;

            h5_int64_t last_step = H5GetNumSteps(file) - 1;
            h5err = H5SetStep(file, last_step);
            if(h5err != H5_SUCCESS)
                ERRORMSG("H5 rc= " << h5err << " in " << __FILE__ << " @ line " << __LINE__ << endl);

            for(int i = 0; i < abs(Nz_diff); ++ i) {
                Nz_read_length[i] = Nz_avrg - signNz;
                Nz_read_start[i] = start;
                start += Nz_read_length[i];
            }
            for(int i = abs(Nz_diff); i < Nnodes; ++ i) {
                Nz_read_length[i] = Nz_avrg;
                Nz_read_start[i] = start;
                start += Nz_read_length[i];
            }
            for(int i = Nnodes; i < Ippl::getNodes(); ++ i) {
                Nz_read_length[i] = 0;
                Nz_read_start[i] = start;
            }
            Nz_read_start[Ippl::getNodes()] = start;

            N_read_start = Nz_read_start[Ippl::myNode()] * num_gridpx_m * num_gridpy_m;

            // rbuf_size = max(Nz_avrg, Nz_avrg - signNz);
            // std::unique_ptr<double> rbuf(new double[Ippl::getNodes() * rbuf_size]);

            h5err = H5Block3dSetView(file,
                                     0, num_gridpx_m - 1,
                                     0, num_gridpy_m - 1,
                                     Nz_read_start[Ippl::myNode()], Nz_read_start[Ippl::myNode() + 1] - 1);
            if(h5err != H5_SUCCESS)
                ERRORMSG("H5 rc= " << h5err << " in " << __FILE__ << " @ line " << __LINE__ << endl);

            field_size = (num_gridpx_m * num_gridpy_m * num_gridpz_m);
            FieldstrengthEx_m.resize(field_size);
            FieldstrengthEy_m.resize(field_size);
            FieldstrengthEz_m.resize(field_size);
            FieldstrengthHx_m.resize(field_size);
            FieldstrengthHy_m.resize(field_size);
            FieldstrengthHz_m.resize(field_size);

            h5err = H5Block3dReadVector3dFieldFloat64(file, "Efield", &(FieldstrengthEx_m[N_read_start]), &(FieldstrengthEy_m[N_read_start]), &(FieldstrengthEz_m[N_read_start]));
            if(h5err != H5_SUCCESS)
                ERRORMSG("H5 rc= " << h5err << " in " << __FILE__ << " @ line " << __LINE__ << endl);
            h5err = H5Block3dReadVector3dFieldFloat64(file, "Hfield", &(FieldstrengthHx_m[N_read_start]), &(FieldstrengthHy_m[N_read_start]), &(FieldstrengthHz_m[N_read_start]));
            if(h5err != H5_SUCCESS)
                ERRORMSG("H5 rc= " << h5err << " in " << __FILE__ << " @ line " << __LINE__ << endl);

            for(int i = 0; i < Nnodes; ++ i) {
                int N_read_start = Nz_read_start[i] * num_gridpx_m * num_gridpy_m;
                int N_read_length = Nz_read_length[i] * num_gridpx_m * num_gridpy_m;
                MPI_Bcast(&(FieldstrengthEx_m[N_read_start]), N_read_length, MPI_DOUBLE, i, Ippl::getComm());
                MPI_Bcast(&(FieldstrengthEy_m[N_read_start]), N_read_length, MPI_DOUBLE, i, Ippl::getComm());
                MPI_Bcast(&(FieldstrengthEz_m[N_read_start]), N_read_length, MPI_DOUBLE, i, Ippl::getComm());
                MPI_Bcast(&(FieldstrengthHx_m[N_read_start]), N_read_length, MPI_DOUBLE, i, Ippl::getComm());
                MPI_Bcast(&(FieldstrengthHy_m[N_read_start]), N_read_length, MPI_DOUBLE, i, Ippl::getComm());
                MPI_Bcast(&(FieldstrengthHz_m[N_read_start]), N_read_length, MPI_DOUBLE, i, Ippl::getComm());
            }

            h5err = H5CloseFile(file);
            if(h5err != H5_SUCCESS)
                ERRORMSG("H5 rc= " << h5err << " in " << __FILE__ << " @ line " << __LINE__ << endl);
            delete[] Nz_read_start;
            delete[] Nz_read_length;


            for(long i = 0; i < num_gridpx_m * num_gridpy_m * num_gridpz_m; i++) {
                FieldstrengthEz_m[i] *= 1.0e6 ;
                FieldstrengthEx_m[i] *= 1.0e6 ;
                FieldstrengthEy_m[i] *= 1.0e6 ;
                FieldstrengthHx_m[i] *= 1.0e6 * mu_0 ;
                FieldstrengthHy_m[i] *= 1.0e6 * mu_0 ;
                FieldstrengthHz_m[i] *= 1.0e6 * mu_0 ;
            }

            msg << typeset_msg("read in fieldmap '" + Filename_m  + "'", "info") << "\n"
                << endl;

        } else {
            WARNMSG("could not read file '" << Filename_m << "'")
        }
    }
}

void FM3DH5Block_nonescale::freeMap() {
    if(!FieldstrengthEz_m.empty()) {
        Inform msg("FM3DH5_NS ");
        FieldstrengthEx_m.clear();
        FieldstrengthEy_m.clear();
        FieldstrengthEz_m.clear();
        FieldstrengthHx_m.clear();
        FieldstrengthHy_m.clear();
        FieldstrengthHz_m.clear();

        msg << typeset_msg("freed fieldmap '" + Filename_m + "'", "info") << "\n"
            << endl;
    }
}

bool FM3DH5Block_nonescale::getFieldstrength(const Vector_t &R, Vector_t &E, Vector_t &B) const {
    const int index_x = static_cast<int>(floor((R(0) - xbegin_m) / hx_m));
    const double lever_x = (R(0) - xbegin_m) / hx_m - index_x;

    const int index_y = static_cast<int>(floor((R(1) - ybegin_m) / hy_m));
    const double lever_y = (R(1) - ybegin_m) / hy_m - index_y;

    const int index_z = (int)floor((R(2) - zbegin_m) / hz_m);
    const double lever_z = (R(2) - zbegin_m) / hz_m - index_z;

    if((index_z < 0) || (index_z + 2 > num_gridpz_m)) {
        return false;
    }

    if(index_x + 2 > num_gridpx_m || index_y + 2 > num_gridpy_m) {
        return false;// do nothing just let E and B unchange in the element with this field map, we check if a particle outside domain in BoundaryGeometry;
    }

    if(index_x < 0 || index_y < 0) {
        //*gmsg<<"Field strength at R "<<R<<" out of field maps"<<endl;
        return false;// if x or y < xmin or ymin, do nothing just let E and B unchange in the element with this field map, we check if a particle outside domain in BoundaryGeometry, not here. That's resonable if we have field mapped without covering entire boundary of an element.
    }
    const long index1 = index_x + (index_y + index_z * num_gridpy_m) * num_gridpx_m;

    E(0) += (1.0 - lever_x) * (1.0 - lever_y) * (1.0 - lever_z) * FieldstrengthEx_m[index1]
            + lever_x           * (1.0 - lever_y) * (1.0 - lever_z) * FieldstrengthEx_m[index1 + 1]
            + (1.0 - lever_x)   * lever_y         * (1.0 - lever_z) * FieldstrengthEx_m[index1 + num_gridpx_m]
            + lever_x           * lever_y         * (1.0 - lever_z) * FieldstrengthEx_m[index1 + num_gridpx_m + 1]
            + (1.0 - lever_x)   * (1.0 - lever_y) * lever_z         * FieldstrengthEx_m[index1 + num_gridpx_m * num_gridpy_m]
            + lever_x           * (1.0 - lever_y) * lever_z         * FieldstrengthEx_m[index1 + num_gridpx_m * num_gridpy_m + 1]
            + (1.0 - lever_x)   * lever_y         * lever_z         * FieldstrengthEx_m[index1 + num_gridpx_m * num_gridpy_m + num_gridpx_m]
            + lever_x           * lever_y         * lever_z         * FieldstrengthEx_m[index1 + num_gridpx_m * num_gridpy_m + num_gridpx_m + 1];

    E(1) += (1.0 - lever_x) * (1.0 - lever_y) * (1.0 - lever_z) * FieldstrengthEy_m[index1]
            + lever_x           * (1.0 - lever_y) * (1.0 - lever_z) * FieldstrengthEy_m[index1 + 1]
            + (1.0 - lever_x)   * lever_y         * (1.0 - lever_z) * FieldstrengthEy_m[index1 + num_gridpx_m]
            + lever_x           * lever_y         * (1.0 - lever_z) * FieldstrengthEy_m[index1 + num_gridpx_m + 1]
            + (1.0 - lever_x)   * (1.0 - lever_y) * lever_z         * FieldstrengthEy_m[index1 + num_gridpx_m * num_gridpy_m]
            + lever_x           * (1.0 - lever_y) * lever_z         * FieldstrengthEy_m[index1 + num_gridpx_m * num_gridpy_m + 1]
            + (1.0 - lever_x)   * lever_y         * lever_z         * FieldstrengthEy_m[index1 + num_gridpx_m * num_gridpy_m + num_gridpx_m]
            + lever_x           * lever_y         * lever_z         * FieldstrengthEy_m[index1 + num_gridpx_m * num_gridpy_m + num_gridpx_m + 1];

    E(2) += (1.0 - lever_x) * (1.0 - lever_y) * (1.0 - lever_z) * FieldstrengthEz_m[index1]
            + lever_x           * (1.0 - lever_y) * (1.0 - lever_z) * FieldstrengthEz_m[index1 + 1]
            + (1.0 - lever_x)   * lever_y         * (1.0 - lever_z) * FieldstrengthEz_m[index1 + num_gridpx_m]
            + lever_x           * lever_y         * (1.0 - lever_z) * FieldstrengthEz_m[index1 + num_gridpx_m + 1]
            + (1.0 - lever_x)   * (1.0 - lever_y) * lever_z         * FieldstrengthEz_m[index1 + num_gridpx_m * num_gridpy_m]
            + lever_x           * (1.0 - lever_y) * lever_z         * FieldstrengthEz_m[index1 + num_gridpx_m * num_gridpy_m + 1]
            + (1.0 - lever_x)   * lever_y         * lever_z         * FieldstrengthEz_m[index1 + num_gridpx_m * num_gridpy_m + num_gridpx_m]
            + lever_x           * lever_y         * lever_z         * FieldstrengthEz_m[index1 + num_gridpx_m * num_gridpy_m + num_gridpx_m + 1];

    B(0) += ((1.0 - lever_x) * (1.0 - lever_y) * (1.0 - lever_z) * FieldstrengthHx_m[index1]
             + lever_x           * (1.0 - lever_y) * (1.0 - lever_z) * FieldstrengthHx_m[index1 + 1]
             + (1.0 - lever_x)   * lever_y         * (1.0 - lever_z) * FieldstrengthHx_m[index1 + num_gridpx_m]
             + lever_x           * lever_y         * (1.0 - lever_z) * FieldstrengthHx_m[index1 + num_gridpx_m + 1]
             + (1.0 - lever_x)   * (1.0 - lever_y) * lever_z         * FieldstrengthHx_m[index1 + num_gridpx_m * num_gridpy_m]
             + lever_x           * (1.0 - lever_y) * lever_z         * FieldstrengthHx_m[index1 + num_gridpx_m * num_gridpy_m + 1]
             + (1.0 - lever_x)   * lever_y         * lever_z         * FieldstrengthHx_m[index1 + num_gridpx_m * num_gridpy_m + num_gridpx_m]
             + lever_x           * lever_y         * lever_z         * FieldstrengthHx_m[index1 + num_gridpx_m * num_gridpy_m + num_gridpx_m + 1]);

    B(1) += ((1.0 - lever_x) * (1.0 - lever_y) * (1.0 - lever_z) * FieldstrengthHy_m[index1]
             + lever_x           * (1.0 - lever_y) * (1.0 - lever_z) * FieldstrengthHy_m[index1 + 1]
             + (1.0 - lever_x)   * lever_y         * (1.0 - lever_z) * FieldstrengthHy_m[index1 + num_gridpx_m]
             + lever_x           * lever_y         * (1.0 - lever_z) * FieldstrengthHy_m[index1 + num_gridpx_m + 1]
             + (1.0 - lever_x)   * (1.0 - lever_y) * lever_z         * FieldstrengthHy_m[index1 + num_gridpx_m * num_gridpy_m]
             + lever_x           * (1.0 - lever_y) * lever_z         * FieldstrengthHy_m[index1 + num_gridpx_m * num_gridpy_m + 1]
             + (1.0 - lever_x)   * lever_y         * lever_z         * FieldstrengthHy_m[index1 + num_gridpx_m * num_gridpy_m + num_gridpx_m]
             + lever_x           * lever_y         * lever_z         * FieldstrengthHy_m[index1 + num_gridpx_m * num_gridpy_m + num_gridpx_m + 1]);

    B(2) += ((1.0 - lever_x) * (1.0 - lever_y) * (1.0 - lever_z) * FieldstrengthHz_m[index1]
             + lever_x           * (1.0 - lever_y) * (1.0 - lever_z) * FieldstrengthHz_m[index1 + 1]
             + (1.0 - lever_x)   * lever_y         * (1.0 - lever_z) * FieldstrengthHz_m[index1 + num_gridpx_m]
             + lever_x           * lever_y         * (1.0 - lever_z) * FieldstrengthHz_m[index1 + num_gridpx_m + 1]
             + (1.0 - lever_x)   * (1.0 - lever_y) * lever_z         * FieldstrengthHz_m[index1 + num_gridpx_m * num_gridpy_m]
             + lever_x           * (1.0 - lever_y) * lever_z         * FieldstrengthHz_m[index1 + num_gridpx_m * num_gridpy_m + 1]
             + (1.0 - lever_x)   * lever_y         * lever_z         * FieldstrengthHz_m[index1 + num_gridpx_m * num_gridpy_m + num_gridpx_m]
             + lever_x           * lever_y         * lever_z         * FieldstrengthHz_m[index1 + num_gridpx_m * num_gridpy_m + num_gridpx_m + 1]);


    return false;
}

bool FM3DH5Block_nonescale::getFieldDerivative(const Vector_t &R, Vector_t &E, Vector_t &B, const DiffDirection &dir) const {
    return false;
}

void FM3DH5Block_nonescale::getFieldDimensions(double &zBegin, double &zEnd, double &rBegin, double &rEnd) const {
    zBegin = zbegin_m;
    zEnd = zend_m;
    rBegin = xbegin_m;
    rEnd = xend_m;
}
void FM3DH5Block_nonescale::getFieldDimensions(double &xIni, double &xFinal, double &yIni, double &yFinal, double &zIni, double &zFinal) const {}

void FM3DH5Block_nonescale::swap() { }

void FM3DH5Block_nonescale::getInfo(Inform *msg) {
    (*msg) << Filename_m << " (3D dynamic); zini= " << zbegin_m << " m; zfinal= " << zend_m << " m;" << endl;
}

double FM3DH5Block_nonescale::getFrequency() const {
    return frequency_m;
}

void FM3DH5Block_nonescale::setFrequency(double freq) {
    frequency_m = freq;
}

void FM3DH5Block_nonescale::getOnaxisEz(vector<pair<double, double> > & F) {
    double Ez_max = 0.0, dz = (zend_m - zbegin_m) / (num_gridpz_m - 1);
    const int index_x = -static_cast<int>(floor(xbegin_m / hx_m));
    const double lever_x = -xbegin_m / hx_m - index_x;

    const int index_y = -static_cast<int>(floor(ybegin_m / hy_m));
    const double lever_y = -ybegin_m / hy_m - index_y;

    long index1 = index_x + index_y * num_gridpx_m;

    F.resize(num_gridpz_m);

    for(int i = 0; i < num_gridpz_m; ++ i) {
        F[i].first = dz * i;
        F[i].second = (1.0 - lever_x)   * (1.0 - lever_y) * FieldstrengthEz_m[index1 + num_gridpx_m * num_gridpy_m]
                      + lever_x           * (1.0 - lever_y) * FieldstrengthEz_m[index1 + num_gridpx_m * num_gridpy_m + 1]
                      + (1.0 - lever_x)   * lever_y         * FieldstrengthEz_m[index1 + num_gridpx_m * num_gridpy_m + num_gridpx_m]
                      + lever_x           * lever_y         * FieldstrengthEz_m[index1 + num_gridpx_m * num_gridpy_m + num_gridpx_m + 1];


        if(fabs(F[i].second) > Ez_max) {
            Ez_max = fabs(F[i].second);
        }
        index1 += num_gridpy_m * num_gridpx_m;
    }

    for(int i = 0; i < num_gridpz_m; ++ i) {
        F[i].second /= Ez_max;
    }
}
