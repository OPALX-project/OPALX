#include "IPPLInitializer.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

#ifdef USENAMESPACE
namespace initializer {
#endif

void IPPLInitializer::indens(){
    long i,j,k,pkindex;
    int proc;
    double rn, rn1, rn2, scal;

    // Generate Gaussian white noise field:
    proc = Parallel.GetMyPE();
    init_random(DataBase.seed, proc);
    for (i=nx1; i<=nx2; ++i){
        for (j=ny1; j<=ny2; ++j){
            for (k=nz1; k<=nz2; ++k){
                do {
                    rn1 = -1 + 2*genrand_real();
                    rn2 = -1 + 2*genrand_real();
                    rn = rn1*rn1 + rn2*rn2;
                } while (rn > 1.0 || rn == 0);
                NDIndex<3> l(Index(i, i), Index(j, j), Index(k, k));
                rho_m.localElement(l) = dcomplex(rn2 * sqrt(-2.0*log(rn)/rn), 0.0);
            }
        }   
    }

    // Go to the k-space:
    fft_m->transform("tok", rho_m);

    // Multiply by the power spectrum:
    scal = pow(DataBase.ngrid, 1.5);
    pkindex = 0;
    for (i=nx1; i<=nx2; ++i){
        for (j=ny1; j<=ny2; ++j){
            for (k=nz1; k<=nz2; ++k){
                NDIndex<3> l(Index(i, i), Index(j, j), Index(k, k));
                rho_m.localElement(l) = dcomplex(std::real(rho_m.localElement(l)) * sqrt(Pk[pkindex])/scal, std::imag(rho_m.localElement(l)) * sqrt(Pk[pkindex])/scal);
                pkindex++;
            }
        }
    }

    // Fix k=0 point:
    if (nx1 == 0 && ny1 == 0 && nz1 == 0) {
        NDIndex<3> l(Index(0, 0), Index(0, 0), Index(0, 0));
        rho_m.localElement(l) = dcomplex(0.0, 0.0);
    }

    return;
}

void IPPLInitializer::test_reality() {
    const long int ngrid=DataBase.ngrid;
    int i,j,k;
    float test_zero, max_rho, min_rho, ave_rho, scal;
    Inform msg ("mc4 ");
    
    min_rho = 1.0e30;
    max_rho = -1.0e30;
    ave_rho = 0.0;
    for (i=0; i<ngx; ++i){
        for (j=0; j<ngy; ++j){
            for (k=0; k<ngz; ++k){

                NDIndex<3> l(Index(nx1+i, nx1+i), Index(ny1+j, ny1+j), Index(nz1+k, nz1+k));
                ave_rho += (real)std::real(rho_m.localElement(l));
                if ((real)std::real(rho_m.localElement(l)) > max_rho) max_rho = (real)std::real(rho_m.localElement(l));
                if ((real)std::real(rho_m.localElement(l)) < min_rho) min_rho = (real)std::real(rho_m.localElement(l));

            }
        }
    }
    ave_rho   = ave_rho/My_Ng;

    MPI_Allreduce(MPI_IN_PLACE, &max_rho, 1, MPI_FLOAT, MPI_MAX, Parallel.GetMpiComm());
    MPI_Allreduce(MPI_IN_PLACE, &min_rho, 1, MPI_FLOAT, MPI_MIN, Parallel.GetMpiComm());
    MPI_Allreduce(MPI_IN_PLACE, &ave_rho, 1, MPI_FLOAT, MPI_SUM, Parallel.GetMpiComm());

    msg << endl << "Min and max value of density in k space: " << min_rho << " " << max_rho << endl;
    ave_rho = ave_rho/NumPEs;
    msg << "Average value of density in k space: " << ave_rho << endl;

    // k to real space!
    fft_m->transform("toreal", rho_m);
    msg << "FFT done ... " << endl;

    scal = 1.0 / pow((float)DataBase.box_size, (float)1.5); // inverse FFT box scaling
    scal *= (ngrid*ngrid*ngrid); //inverse FFT scaling
    for (i=0; i<ngx; ++i){
        for (j=0; j<ngy; ++j){
            for (k=0; k<ngz; ++k){
                NDIndex<3> l(Index(nx1+i, nx1+i), Index(ny1+j, ny1+j), Index(nz1+k, nz1+k));
                rho_m.localElement(l) = dcomplex(std::real(rho_m.localElement(l))*scal, std::imag(rho_m.localElement(l)));
            }
        }
    }

    min_rho = 1.0e30;
    max_rho = -1.0e30;
    test_zero = 0.0;
    ave_rho = 0.0;
    for (i=nx1; i<=nx2; ++i){
        for (j=ny1; j<=ny2; ++j){
            for (k=nz1; k<=nz2; ++k){
                NDIndex<3> l(Index(i, i), Index(j, j), Index(k, k));
                if (test_zero < fabs((real)std::imag(rho_m.localElement(l))))
                    test_zero = fabs((real)std::imag(rho_m.localElement(l)));
                ave_rho += (real)std::real(rho_m.localElement(l));
                if ((real)std::real(rho_m.localElement(l)) > max_rho) max_rho = (real)std::real(rho_m.localElement(l));
                if ((real)std::real(rho_m.localElement(l)) < min_rho) min_rho = (real)std::real(rho_m.localElement(l));
            }
        }
    }
    ave_rho = ave_rho/My_Ng;
    
    MPI_Allreduce(MPI_IN_PLACE, &test_zero, 1, MPI_FLOAT, MPI_MAX, Parallel.GetMpiComm());	
    MPI_Allreduce(MPI_IN_PLACE, &max_rho, 1, MPI_FLOAT, MPI_MAX, Parallel.GetMpiComm());
    MPI_Allreduce(MPI_IN_PLACE, &min_rho, 1, MPI_FLOAT, MPI_MIN, Parallel.GetMpiComm());
    MPI_Allreduce(MPI_IN_PLACE, &ave_rho, 1, MPI_FLOAT, MPI_SUM, Parallel.GetMpiComm());

    msg << endl << "Max value of the imaginary part of density is " << test_zero << endl;
    ave_rho = ave_rho/NumPEs;
    msg << "Average value of density in code units: " << ave_rho << endl;
    return;
}

void IPPLInitializer::gravity_force(integer axis) { // 0=x, 1=y, 2=z
    long i, j, k, k_i, k_j, k_k, index;
    const real tpi=2.0*pi;
    const integer ngrid=DataBase.ngrid;
    const integer nq=ngrid/2;
    real kk, Grad, Green, scal;
    real xx, yy, zz;
    float min_F, max_F, ave_F, test_zero;
    dcomplex temp;
    Inform msg ("mc4 ");

    /* Multiply by the Green's function (-1/k^2), 
       and take the gradient in k-space (-ik_i) */
    switch (axis) {
    case 0:
        xx = 1.0;
        yy = 0.0;
        zz = 0.0;
        break;
    case 1:
        xx = 0.0;
        yy = 1.0;
        zz = 0.0;
        break;
    case 2:
        xx = 0.0;
        yy = 0.0;
        zz = 1.0;
        break;
    }
    for (i=0; i<ngx; ++i){
        k_i = i+nx1;
        if (k_i >= nq) {k_i = -MOD(ngrid-k_i,ngrid);}
        for (j=0; j<ngy; ++j){
            k_j = j+ny1;
            if (k_j >= nq) {k_j = -MOD(ngrid-k_j,ngrid);}
            for (k=0; k<ngz; ++k){
                k_k = k+nz1;
                if (k_k >= nq) {k_k = -MOD(ngrid-k_k,ngrid);}
                Grad = -tpi/ngrid*(k_i*xx+k_j*yy+k_k*zz);
                kk = tpi/ngrid*sqrt(k_i*k_i+k_j*k_j+k_k*k_k);
                Green = -1.0/(kk*kk);
                if (kk == 0.0) Green = 0.0;
                NDIndex<3> l(Index(nx1+i, nx1+i), Index(ny1+j, ny1+j), Index(nz1+k, nz1+k));
                temp = rho_m.localElement(l);
                rho_m.localElement(l) = dcomplex(-Grad * Green * std::imag(temp), Grad * Green * std::real(temp));
            }
        }
    }
    // Fix k=0 point:
    if (nx1 == 0 && ny1 == 0 && nz1 == 0) {
        NDIndex<3> l(Index(0, 0), Index(0, 0), Index(0, 0));
        rho_m.localElement(l) = dcomplex(0.0, std::imag(rho_m.localElement(l)));
    }

    fft_m->transform("toreal", rho_m);
    msg << "FFT done " << endl;

    /* compensation for inverse FFT scaling */
    scal = 1.0 / pow((float)DataBase.box_size, (float)1.5);
    scal *= ngrid*ngrid*ngrid;
    for (i=0; i<ngx; ++i){
        for (j=0; j<ngy; ++j){
            for (k=0; k<ngz; ++k){
                NDIndex<3> l(Index(nx1+i, nx1+i), Index(ny1+j, ny1+j), Index(nz1+k, nz1+k));
                rho_m.localElement(l) = dcomplex(std::real(rho_m.localElement(l))*scal, std::imag(rho_m.localElement(l)));
            }
        }
    }
}

void IPPLInitializer::gravity_potential() {
    long i, j, k, k_i, k_j, k_k, index;
    const real tpi=2.0*pi;
    const integer ngrid=DataBase.ngrid;
    const integer nq=ngrid/2;
    real kk, Green, scal;
    float min_phi, max_phi, ave_phi, test_zero;
    Inform msg ("mc4 ");

    /* Multiply by the Green's function (-1/k^2) */
    for (i=0; i<ngx; ++i){
        k_i = i+nx1;
        if (k_i >= nq) {k_i = -MOD(ngrid-k_i,ngrid);}
        for (j=0; j<ngy; ++j){
            k_j = j+ny1;
            if (k_j >= nq) {k_j = -MOD(ngrid-k_j,ngrid);}
            for (k=0; k<ngz; ++k){
                k_k = k+nz1;
                if (k_k >= nq) {k_k = -MOD(ngrid-k_k,ngrid);}
                index = (i*ngy+j)*ngz+k;
                kk = tpi/ngrid*sqrt(k_i*k_i+k_j*k_j+k_k*k_k);
                Green = -1.0/(kk*kk);
                if (kk == 0.0) Green = 0.0;
                NDIndex<3> l(Index(nx1+i, nx1+i), Index(ny1+j, ny1+j), Index(nz1+k, nz1+k));
                rho_m.localElement(l) *= dcomplex(Green, Green);
            }
        }
    }
    // Fix k=0 point:
    if (nx1 == 0 && ny1 == 0 && nz1 == 0) {
        NDIndex<3> l(Index(0, 0), Index(0, 0), Index(0, 0));
        rho_m.localElement(l) = dcomplex(0.0, std::imag(rho_m.localElement(l)));
    }

    fft_m->transform("toreal", rho_m);
    msg  << "FFT done " << endl;

    /* compensation for inverse FFT scaling */
    scal = 1.0 / pow((float)DataBase.box_size, (float)1.5);
    scal *= ngrid*ngrid*ngrid;
    for (i=0; i<ngx; ++i){
        for (j=0; j<ngy; ++j){
            for (k=0; k<ngz; ++k){
                NDIndex<3> l(Index(nx1+i, nx1+i), Index(ny1+j, ny1+j), Index(nz1+k, nz1+k));
                rho_m.localElement(l) = dcomplex(std::real(rho_m.localElement(l))*scal, std::imag(rho_m.localElement(l)));
            }   
        }
    }

    return;
}


void IPPLInitializer::non_gaussian(integer axis) { // 0=x, 1=y, 2=z 
    const real f_NL=DataBase.f_NL;
    const integer ngrid=DataBase.ngrid;
    const integer nq=ngrid/2;
    const real tpi=2.0*pi;
    const real tpiL=tpi/DataBase.box_size; // k0, physical units
    long i, j, k, k_i, k_j, k_k, index;
    float ave_phi2;
    float xx, yy, zz;
    real Grad, kk, trans_f, fff;
    dcomplex temp;
    CxField_t phi;
    Inform msg ("mc4 ");

    phi = rho_m; /* After gravity solve this array holds potential */
    // Find average value for Phi^2:
    ave_phi2 = 0.0;
    for (i=nx1; i<=nx2; ++i){
        for (j=ny1; j<=ny2; ++j){
            for (k=nz1; k<=nz2; ++k){
                NDIndex<3> l(Index(i, i), Index(j, j), Index(k, k));
                ave_phi2 += (float)(std::real(phi.localElement(l))*std::real(phi.localElement(l)));
            }
        }
    }
    ave_phi2 = ave_phi2/My_Ng;
    MPI_Allreduce(MPI_IN_PLACE, &ave_phi2, 1, MPI_FLOAT, MPI_SUM, Parallel.GetMpiComm());
    ave_phi2 = ave_phi2/NumPEs;

    // Add f_NL:
    fff = 1.0;
    for (i=nx1; i<=nx2; ++i){
        for (j=ny1; j<=ny2; ++j){
            for (k=nz1; k<=nz2; ++k){
                NDIndex<3> l(Index(i, i), Index(j, j), Index(k, k));
                phi.localElement(l) += dcomplex(f_NL*(std::real(phi.localElement(l)) * std::real(phi.localElement(l)) - ave_phi2)*fff, 0.0);
            }
        }
    }

    // Go back to k-space (no forward scaling):
    fft_m->transform("tok", phi);
    msg << "FFT done " << endl;

    // Add transfer function and gradient in k-space:
    switch (axis) {
    case 0:
        xx = 1.0;
        yy = 0.0;
        zz = 0.0;
        break;
    case 1:
        xx = 0.0;
        yy = 1.0;
        zz = 0.0;
        break;
    case 2:
        xx = 0.0;
        yy = 0.0;
        zz = 1.0;
        break;
    }
    for (i=0; i<ngx; ++i){
        k_i = i+nx1;
        if (k_i >= nq) {k_i = -MOD(ngrid-k_i,ngrid);}
        for (j=0; j<ngy; ++j){
            k_j = j+ny1;
            if (k_j >= nq) {k_j = -MOD(ngrid-k_j,ngrid);}
            for (k=0; k<ngz; ++k){
                k_k = k+nz1;
                if (k_k >= nq) {k_k = -MOD(ngrid-k_k,ngrid);}
                kk = tpiL*sqrt(k_i*k_i+k_j*k_j+k_k*k_k);
                trans_f = Cosmology.TransferFunction(kk);
                Grad = -tpi/ngrid*(k_i*xx+k_j*yy+k_k*zz);
                NDIndex<3> l(Index(nx1+i, nx1+i), Index(ny1+j, ny1+j), Index(nz1+k, nz1+k));
                temp = phi.localElement(l);
                phi.localElement(l) = dcomplex(-Grad * trans_f * std::imag(temp), Grad * trans_f * std::real(temp));
            }
        }
    }
    if (nx1 == 0 && ny1 == 0 && nz1 == 0) {
        NDIndex<3> l(Index(0, 0), Index(0, 0), Index(0, 0));
        phi.localElement(l) = dcomplex(0.0, std::imag(phi.localElement(l))); // fix k=0
    }

    // FFT to real space to get the force:
    fft_m->transform("toreal", phi);
    msg << "FFT done " << endl;

    //FIXME:!!
    real scal = pow(ngrid, 3.0);  /* inverse FFT normalization */
    for (i=nx1; i<=nx2; ++i){
        for (j=ny1; j<=ny2; ++j){
            for (k=nz1; k<=nz2; ++k){
                NDIndex<3> l(Index(i, i), Index(j, j), Index(k, k));
                phi.localElement(l) = dcomplex(std::real(phi.localElement(l))/scal, std::imag(phi.localElement(l)));
            }
        }
    }

    return;
}

void IPPLInitializer::set_particles(real z_in, real d_z, real ddot, integer axis){ // 0=x, 1=y, 2=z
   const integer ngrid=DataBase.ngrid;
   long i, j, k, index;
   real pos_0, xx, yy, zz;
   float  move, max_move, ave_move;
      
   switch (axis) {
      case 0:
         xx = 1.0;
         yy = 0.0;
         zz = 0.0;
         break;
      case 1:
         xx = 0.0;
         yy = 1.0;
         zz = 0.0;
         break;
      case 2:
         xx = 0.0;
         yy = 0.0;
         zz = 1.0;
         break;
   }

   /* Move particles according to Zeldovich approximation 
   particles will be put in rho array; real array will 
   hold positions, imaginary will hold velocities */

   for (i=0; i<ngx; ++i){
      for (j=0; j<ngy; ++j){
         for (k=0; k<ngz; ++k){
            index = (i*ngy+j)*ngz+k;
            pos_0 = (i+nx1)*xx + (j+ny1)*yy + (k+nz1)*zz;
            NDIndex<3> l(Index(nx1+i, nx1+i), Index(ny1+j, ny1+j), Index(nz1+k, nz1+k));
            dcomplex tmp = rho_m.localElement(l);
            rho_m.localElement(l) = dcomplex(pos_0 - d_z*std::real(tmp), -ddot*std::real(tmp));
         }
      }
   }
      
   return;
}

void IPPLInitializer::output(integer axis, real* pos, real* vel){ // 0=x, 1=y, 2=z
   long i, j, k, index;
      
   index = 0;
   for (i=0; i<ngx; ++i){
       for (j=0; j<ngy; ++j){
           for (k=0; k<ngz; ++k){
               NDIndex<3> l(Index(nx1+i, nx1+i), Index(ny1+j, ny1+j), Index(nz1+k, nz1+k));
               dcomplex val = rho_m.localElement(l);
               pos[index] = (real)std::real(val);
               vel[index] = (real)std::imag(val);
               index++;
           }
        }
    }
      
   return;
}


//template <class T, unsigned int Dim>
void IPPLInitializer::init_particles(ParticleAttrib< Vektor<T, 3> > &pos, 
            ParticleAttrib< Vektor<T, 3> > &vel, FieldLayout_t *layout, Mesh_t *mesh,
            InputParser& par, const char *tfName, MPI_Comm comm) 
{
    integer i;
    real d_z, ddot, f_NL;
    double t1, t2;
    Inform msg ("mc4 ");

    // Preliminaries:
    StartMonitor();
    OnClock("Initialization");
    // XXX: we use our own communicator
    Parallel.InitMPI(&t1, comm);
    MyPE = Parallel.GetMyPE();
    MasterPE = Parallel.GetMasterPE();
    NumPEs = Parallel.GetNumPEs();
    //FIXME
    //DataBase.GetParameters(bdata, Parallel);
    DataBase.GetParameters(par, Parallel);
    Cosmology.SetParameters(DataBase, tfName);

    this->flayout = layout;
    this->mesh = mesh;
    // create additional objects for the use in the FFT's
    rho_m.initialize(*mesh, *flayout, GuardCellSizes<3>(1));

    // Find what part of the global grid i got:
    NDIndex<3> localidx =  flayout->getLocalNDIndex();
    nx1 = localidx[0].first();
    nx2 = localidx[0].last();
    ny1 = localidx[1].first();
    ny2 = localidx[1].last();
    nz1 = localidx[2].first();
    nz2 = localidx[2].last();

    ngx = localidx[0].length();
    ngy = localidx[1].length();
    ngz = localidx[2].length();

    My_Ng = ngx*ngy*ngz;

    // Allocate arrays:
    Pk   = (real *)malloc(My_Ng*sizeof(T));

    // create the FFT object
    bool compressTemps = true;
    fft_m = new FFT_t(flayout->getDomain(), compressTemps);
    fft_m->setDirectionName(+1, "toreal");
    fft_m->setDirectionName(-1, "tok");

    OffClock("Initialization");

    // Set P(k) array:
    OnClock("Creating Pk(x,y,z)");
    initspec();
    OffClock("Creating Pk(x,y,z)");

    // Growth factor for the initial redshift:
    Cosmology.GrowthFactor(DataBase.z_in, &d_z, &ddot);

    msg << "redshift: " << DataBase.z_in << "; growth factor = " << d_z << "; ";
    msg << "derivative = " << ddot << endl;


    f_NL=DataBase.f_NL;
    for (i=0; i<Ndim; ++i){
        OnClock("Creating rho(x,y,z)");
        indens();          // Set density in k-space:
        OffClock("Creating rho(x,y,z)");
#ifdef TESTREALITY
        // Check if density field is real:
        test_reality();
#endif

        if (f_NL == 0.0){ // Solve for the force immediately 
            OnClock("Poisson solve");
            gravity_force(i);  // Calculate i-th component of the force
            OffClock("Poisson solve");
        }
        else { // Solve for Phi, apply f_NL, apply T(k), solve for the force
            OnClock("Poisson solve");
            gravity_potential();  // Calculate gravitational potential
            OffClock("Poisson solve");
            OnClock("Non-Gaussianity");
            non_gaussian(i);  // Create non-gaussian force
            OffClock("Non-Gaussianity");
        }
   
        OnClock("Particle move");
        set_particles(DataBase.z_in, d_z, ddot, i); // Zeldovich move
        OffClock("Particle move");

        real *tpos = (real*)malloc(My_Ng*sizeof(T));
        real *tvel = (real*)malloc(My_Ng*sizeof(T));

        if (i == 0) {
            OnClock("Output");
	    output(i, tpos, tvel);
            OffClock("Output");
            if (DataBase.Omega_nu > 0.0){
                OnClock("Neutrinos");
                buf_pos = (T *)malloc(ngy*ngz*sizeof(T));
                buf_vel = (T *)malloc(ngy*ngz*sizeof(T));
                exchange_buffers(tpos, tvel);
                inject_neutrinos(tpos, tvel, i);
                OffClock("Neutrinos");
            }
            apply_periodic(tpos);
        }   // Put particle data
        else if (i == 1) {
            OnClock("Output");
            output(i, tpos, tvel);
            OffClock("Output");
            if (DataBase.Omega_nu > 0.0){
                OnClock("Neutrinos");
                exchange_buffers(tpos, tvel);
                inject_neutrinos(tpos, tvel, i);
                OffClock("Neutrinos");
            }
            apply_periodic(tpos);
        }
        else {
            OnClock("Output");
            output(i, tpos, tvel);
            OffClock("Output");
            if (DataBase.Omega_nu > 0.0){
                OnClock("Neutrinos");
                exchange_buffers(tpos, tvel);
                inject_neutrinos(tpos, tvel, i);
                free(buf_pos);
                free(buf_vel);
                OffClock("Neutrinos");
            }
            apply_periodic(tpos);
        }
	
        //FIXME: store particles to IPPL container + NEUTRINOS
        //FIXME: use pos_i, vel_i

	for(size_t j=0; j<My_Ng; ++j){
	    pos[j](i) = tpos[j];
	    vel[j](i) = tvel[j];
	}
	free(tpos);
	free(tvel);
    }
    //FIXME
    // Add thermal velocity to neutrino particles:
    //if (DataBase.Omega_nu > 0.0) thermal_vel(vel_x, vel_y, vel_z);


    // Clean all the allocations:
    free(Pk);

    MPI_Barrier(comm);
    if (MyPE == MasterPE) PrintClockSummary(stdout);
    Parallel.FinalMPI(&t2);	

}

#ifdef USENAMESPACE
}
#endif

