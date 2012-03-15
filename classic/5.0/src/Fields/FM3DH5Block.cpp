#include <fstream>
#include <ios>

#include "Fields/FM3DH5Block.hh"
#include "Fields/Fieldmap.icc"
#include "H5Part.h"
#include "H5Block.h"
#include "Physics/Physics.h"

extern Inform *gmsg;

using namespace std;
using Physics::mu_0;

FM3DH5Block::FM3DH5Block(string aFilename):
    Fieldmap(aFilename),
    FieldstrengthEx_m(NULL),
    FieldstrengthEy_m(NULL),
    FieldstrengthEz_m(NULL),
    FieldstrengthHx_m(NULL),
    FieldstrengthHy_m(NULL),
    FieldstrengthHz_m(NULL) {
    Inform msg("FM3DH5 ");
    h5part_int64_t h5err;
    h5part_int64_t grid_rank;
    h5part_int64_t grid_dims[3];
    h5part_int64_t field_dims;
    char name[20];
    h5part_int64_t len_name = 20;
    long field_size = 0;

    Type = T3DDynamicH5Block;
    H5PartFile *file = H5PartOpenFileParallel(aFilename.c_str(), H5PART_READ, MPI_COMM_WORLD);


    if(file) {
        h5part_int64_t last_step = H5PartGetNumSteps(file) - 1;
        H5PartSetStep(file, last_step);
        h5part_int64_t num_fields = H5BlockGetNumFields(file);


        for(h5part_int64_t i = 0; i < num_fields; ++ i) {
            /*
              Work around API changes in H5Block
            */
#if H5PART_VER_MINOR >= 6
            h5part_int64_t ftype;
            h5err = H5BlockGetFieldInfo(file, i, name, len_name, &grid_rank, grid_dims, &field_dims, &ftype);
#else
            h5err = H5BlockGetFieldInfo(file, i, name, len_name, &grid_rank, grid_dims, &field_dims);
#endif
            if(strcmp(name, "Efield") == 0) {
                num_gridpx_m = grid_dims[0];
                num_gridpy_m = grid_dims[1];
                num_gridpz_m = grid_dims[2];
            }
        }
        h5err = H5Block3dGetFieldSpacing(file, "Efield", &hx_m, &hy_m, &hz_m);
        h5err = H5Block3dGetFieldOrigin(file, "Efield", &xbegin_m, &ybegin_m, &zbegin_m);
        xend_m = xbegin_m + (num_gridpx_m - 1) * hx_m;
        yend_m = ybegin_m + (num_gridpy_m - 1) * hy_m;
        zend_m = zbegin_m + (num_gridpz_m - 1) * hz_m;

        //         xcentral_idx_m = static_cast<int>(fabs(xbegin_m) / hx_m);
        //         ycentral_idx_m = static_cast<int>(fabs(ybegin_m) / hy_m);


        h5err = H5PartReadFileAttrib(file, "Resonance Frequency(Hz)", &frequency_m);
        frequency_m *= Physics::two_pi;

        H5PartCloseFile(file);
    }
}


FM3DH5Block::~FM3DH5Block() {
    if(FieldstrengthEz_m != NULL) {
        delete[] FieldstrengthEz_m;
        delete[] FieldstrengthEx_m;
        delete[] FieldstrengthEy_m;
        delete[] FieldstrengthHz_m;
        delete[] FieldstrengthHx_m;
        delete[] FieldstrengthHy_m;
    }
}

void FM3DH5Block::readMap() {
    if(FieldstrengthEz_m == NULL) {
        Inform msg("FM3DH5 ");
        H5PartFile *file = H5PartOpenFileParallel(Filename_m.c_str(), H5PART_READ, MPI_COMM_WORLD);

        if(file) {
            h5part_int64_t h5err;
            h5part_int64_t grid_rank;
            h5part_int64_t grid_dims[3];
            h5part_int64_t field_dims;
            h5part_int64_t len_name = 20;
            char name[len_name];

            long field_size = 0;
            int Nnodes = Ippl::getNodes();//min(20, Ippl::getNodes());
            int Nz_avrg = static_cast<int>(floor(0.5 + num_gridpz_m / Nnodes));
            int Nz_diff = Nz_avrg * Nnodes - num_gridpz_m;
            int signNz = Nz_diff > 0 ? 1 : -1;
            int *Nz_read_start = new int[Ippl::getNodes() + 1];
            int *Nz_read_length = new int[Ippl::getNodes()];
            int N_read_start;
            int start = 0;
            int rbuf_size;
            double *rbuf;

            long ii;
            int index_x;
            int index_y;
            double lever_x;
            double lever_y;
            double Ezmax = 0.0;

            h5part_int64_t last_step = H5PartGetNumSteps(file) - 1;
            H5PartSetStep(file, last_step);

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

            rbuf_size = max(Nz_avrg, Nz_avrg - signNz);
            rbuf = (double *) malloc(Ippl::getNodes() * rbuf_size * sizeof(double));

            h5err = H5BlockDefine3DFieldLayout(file,
                                               0, num_gridpx_m - 1,
                                               0, num_gridpy_m - 1,
                                               Nz_read_start[Ippl::myNode()], Nz_read_start[Ippl::myNode() + 1] - 1);

            field_size = (num_gridpx_m * num_gridpy_m * num_gridpz_m);
            FieldstrengthEx_m = new h5part_float64_t[field_size];
            FieldstrengthEy_m = new h5part_float64_t[field_size];
            FieldstrengthEz_m = new h5part_float64_t[field_size];
            FieldstrengthHx_m = new h5part_float64_t[field_size];
            FieldstrengthHy_m = new h5part_float64_t[field_size];
            FieldstrengthHz_m = new h5part_float64_t[field_size];
#if H5PART_VER_MINOR >= 6
            h5err = H5Block3dRead3dVectorFieldFloat64(file, "Efield", &(FieldstrengthEx_m[N_read_start]), &(FieldstrengthEy_m[N_read_start]), &(FieldstrengthEz_m[N_read_start]));
            h5err = H5Block3dRead3dVectorFieldFloat64(file, "Hfield", &(FieldstrengthHx_m[N_read_start]), &(FieldstrengthHy_m[N_read_start]), &(FieldstrengthHz_m[N_read_start]));
#else
            h5err = H5Block3dRead3dVectorField(file, "Efield", &(FieldstrengthEx_m[N_read_start]), &(FieldstrengthEy_m[N_read_start]), &(FieldstrengthEz_m[N_read_start]));
            h5err = H5Block3dRead3dVectorField(file, "Hfield", &(FieldstrengthHx_m[N_read_start]), &(FieldstrengthHy_m[N_read_start]), &(FieldstrengthHz_m[N_read_start]));
#endif

            for(int i = 0; i < Nnodes; ++ i) {
                int N_read_start = Nz_read_start[i] * num_gridpx_m * num_gridpy_m;
                int N_read_length = Nz_read_length[i] * num_gridpx_m * num_gridpy_m;
                MPI_Bcast(&(FieldstrengthEx_m[N_read_start]), N_read_length, MPI_DOUBLE, i, MPI_COMM_WORLD);
                MPI_Bcast(&(FieldstrengthEy_m[N_read_start]), N_read_length, MPI_DOUBLE, i, MPI_COMM_WORLD);
                MPI_Bcast(&(FieldstrengthEz_m[N_read_start]), N_read_length, MPI_DOUBLE, i, MPI_COMM_WORLD);
                MPI_Bcast(&(FieldstrengthHx_m[N_read_start]), N_read_length, MPI_DOUBLE, i, MPI_COMM_WORLD);
                MPI_Bcast(&(FieldstrengthHy_m[N_read_start]), N_read_length, MPI_DOUBLE, i, MPI_COMM_WORLD);
                MPI_Bcast(&(FieldstrengthHz_m[N_read_start]), N_read_length, MPI_DOUBLE, i, MPI_COMM_WORLD);
            }

            H5PartCloseFile(file);

            delete[] Nz_read_start;
            delete[] Nz_read_length;

            index_x = static_cast<int>(floor(-xbegin_m / hx_m));
            index_y = static_cast<int>(floor(-ybegin_m / hy_m));
            lever_x = -xbegin_m / hx_m - index_x;
            lever_y = -ybegin_m / hy_m - index_y;


            ii = index_x + index_y * num_gridpx_m;
            for(int i = 0; i < num_gridpz_m; i++) {
                double E = fabs((1. - lever_x) * (1. - lever_y) * FieldstrengthEz_m[ii] +
                                lever_x * (1. - lever_y) * FieldstrengthEz_m[ii+1] +
                                (1. - lever_x) * lever_y * FieldstrengthEz_m[ii+num_gridpx_m] +
                                lever_x * lever_y * FieldstrengthEz_m[ii+num_gridpx_m+1]);

                if(E > Ezmax) {
                    Ezmax = E;
                }
                ii += num_gridpx_m * num_gridpy_m;
            }

            for(long i = 0; i < num_gridpx_m * num_gridpy_m * num_gridpz_m; i++) {
                FieldstrengthEz_m[i] *= 1.0e6 / Ezmax;
                FieldstrengthEx_m[i] *= 1.0e6 / Ezmax;
                FieldstrengthEy_m[i] *= 1.0e6 / Ezmax;
                FieldstrengthHx_m[i] *= 1.0e6 * mu_0 / Ezmax;
                FieldstrengthHy_m[i] *= 1.0e6 * mu_0 / Ezmax;
                FieldstrengthHz_m[i] *= 1.0e6 * mu_0 / Ezmax;
            }

            msg << typeset_msg("read in fieldmap '" + Filename_m  + "'", "info") << "\n"
                << endl;
        } else {
            WARNMSG("could not read file '" << Filename_m << "'")
        }
    }
}

void FM3DH5Block::freeMap() {
    if(FieldstrengthEz_m != NULL) {
        Inform msg("FM3DD ");
        delete[] FieldstrengthEx_m;
        delete[] FieldstrengthEy_m;
        delete[] FieldstrengthEz_m;
        delete[] FieldstrengthHx_m;
        delete[] FieldstrengthHy_m;
        delete[] FieldstrengthHz_m;

        msg << typeset_msg("freed fieldmap '" + Filename_m + "'", "info") << "\n"
            << endl;
    }
}

bool FM3DH5Block::getFieldstrength(const Vector_t &R, Vector_t &E, Vector_t &B) const {
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
        return true;
    }

    const long index1 = index_x + (index_y + index_z * num_gridpy_m) * num_gridpx_m;
    const long index2 = index1 + num_gridpz_m;

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

bool FM3DH5Block::getFieldstrength_fdiff(const Vector_t &R, Vector_t &E, Vector_t &B, const DiffDirection &dir) const {

}

void FM3DH5Block::getFieldDimensions(double &zBegin, double &zEnd, double &rBegin, double &rEnd) const {
    zBegin = zbegin_m;
    zEnd = zend_m;
    rBegin = xbegin_m;
    rEnd = xend_m;
}

void FM3DH5Block::swap() {
    //     if (swap_m) swap_m = false;
    //     else swap_m = true;
}

void FM3DH5Block::getInfo(Inform *msg) {
    (*msg) << Filename_m << " (3D dynamic); zini= " << zbegin_m << " m; zfinal= " << zend_m << " m;" << endl;
}

double FM3DH5Block::getFrequency() const {
    return frequency_m;
}

void FM3DH5Block::setFrequency(double freq) {
    frequency_m = freq;
}

void FM3DH5Block::getOnaxisEz(vector<pair<double, double> > & F) {
    double Ez_max = 0.0, dz = (zend_m - zbegin_m) / (num_gridpz_m - 1);
    double tmpDouble;
    int tmpInt;
    string tmpString;
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
