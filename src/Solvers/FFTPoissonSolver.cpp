// -*- C++ -*-
/***************************************************************************
 *
 *
 * FFTPoissonSolver.cc
 *
 *
 *
 *
 *
 *
 *
 ***************************************************************************/

// include files
#include "Solvers/FFTPoissonSolver.h"
#include "Algorithms/PartBunch.h"
#include "Physics/Physics.h"
#include "Utility/IpplTimings.h"
#include "BasicActions/Option.h"
#include "Utilities/Options.h"
#include <fstream>
//////////////////////////////////////////////////////////////////////////////
// a little helper class to specialize the action of the Green's function
// calculation.  This should be specialized for each dimension
// to the proper action for computing the Green's function.  The first
// template parameter is the full type of the Field to compute, and the second
// is the dimension of the data, which should be specialized.

#ifdef OPAL_NOCPLUSPLUS11_NULLPTR
#define nullptr NULL
#endif

template<unsigned int Dim>
struct SpecializedGreensFunction { };

template<>
struct SpecializedGreensFunction<3> {
    template<class T, class FT, class FT2>
    static void calculate(Vektor<T, 3> &hrsq, FT &grn, FT2 *grnI) {
        grn = grnI[0] * hrsq[0] + grnI[1] * hrsq[1] + grnI[2] * hrsq[2];
        grn = 1.0 / sqrt(grn);
        grn[0][0][0] = grn[0][0][1];
    }
};

////////////////////////////////////////////////////////////////////////////

// constructor


FFTPoissonSolver::FFTPoissonSolver(Mesh_t *mesh, FieldLayout_t *fl, std::string greensFunction, std::string bcz):
    mesh_m(mesh),
    layout_m(fl),
    mesh2_m(nullptr),
    layout2_m(nullptr),
    mesh3_m(nullptr),
    layout3_m(nullptr),
    mesh4_m(nullptr),
    layout4_m(nullptr),
    greensFunction_m(greensFunction)
{
    int i;

    bcz_m = (bcz==std::string("PERIODIC"));   // for DC beams, the z direction has periodic boundary conditions

    domain_m = layout_m->getDomain();

    // For efficiency in the FFT's, we can use a parallel decomposition
    // which can be serial in the first dimension.
    e_dim_tag decomp[3];
    e_dim_tag decomp2[3];
    for(int d = 0; d < 3; ++d) {
        decomp[d] = layout_m->getRequestedDistribution(d);
        decomp2[d] = layout_m->getRequestedDistribution(d);
    }

    if (bcz_m) {
        // The FFT's require double-sized field sizes in order to
        // simulate an isolated system.  The FFT of the charge density field, rho,
        // would otherwise mimic periodic boundary conditions, i.e. as if there were
        // several beams set a periodic distance apart.  The double-sized fields in x and
        // alleviate this problem, in z we have periodic BC's
        for(i = 0; i < 2; i++) {
            hr_m[i] = mesh_m->get_meshSpacing(i);
            nr_m[i] = domain_m[i].length();
            domain2_m[i] = Index(2 * nr_m[i] + 1);
        }

        hr_m[2] = mesh_m->get_meshSpacing(2);
        nr_m[2] = domain_m[2].length();
        domain2_m[2] = Index(2 * nr_m[2] + 1);

        for(i = 0; i < 2 * 3; ++i) {
            bc_m[i] = new ZeroFace<double, 3, Mesh_t, Center_t>(i);
            vbc_m[i] = new ZeroFace<Vector_t, 3, Mesh_t, Center_t>(i);
        }
        // z-direction
        bc_m[4] = new ParallelPeriodicFace<double,3,Mesh_t,Center_t>(4);
        bc_m[5] = new ParallelPeriodicFace<double,3,Mesh_t,Center_t>(5);
        vbc_m[4] = new ZeroFace<Vector_t, 3, Mesh_t, Center_t>(4);
        vbc_m[5] = new ZeroFace<Vector_t, 3, Mesh_t, Center_t>(5);
    }
    else {
        // The FFT's require double-sized field sizes in order to
        // simulate an isolated system.  The FFT of the charge density field, rho,
        // would otherwise mimic periodic boundary conditions, i.e. as if there were
        // several beams set a periodic distance apart.  The double-sized fields
        // alleviate this problem.
        for(i = 0; i < 3; i++) {
            hr_m[i] = mesh_m->get_meshSpacing(i);
            nr_m[i] = domain_m[i].length();
            domain2_m[i] = Index(2 * nr_m[i] + 1);
        }

        for(i = 0; i < 2 * 3; ++i) {
            bc_m[i] = new ZeroFace<double, 3, Mesh_t, Center_t>(i);
            vbc_m[i] = new ZeroFace<Vector_t, 3, Mesh_t, Center_t>(i);
        }
    }

    // create double sized mesh and layout objects for the use in the FFT's
    mesh2_m = std::unique_ptr<Mesh_t>(new Mesh_t(domain2_m));
    layout2_m = std::unique_ptr<FieldLayout_t>(new FieldLayout_t(*mesh2_m, decomp));
#ifdef OPAL_DKS
    rho2_m.initialize(*mesh2_m, *layout2_m, false);
#else
    rho2_m.initialize(*mesh2_m, *layout2_m);
#endif


    NDIndex<3> tmpdomain;
    // Create the domain for the transformed (complex) fields.  Do this by
    // taking the domain from the doubled mesh, permuting it to the right, and
    // setting the 2nd dimension to have n/2 + 1 elements.
    domain3_m[0] = Index(2 * nr_m[3-1] + 1);
    domain3_m[1] = Index(nr_m[0] + 2);

    for(i = 2; i < 3; ++i)
        domain3_m[i] = Index(2 * nr_m[i-1] + 1);

    // create mesh and layout for the new real-to-complex FFT's, for the
    // complex transformed fields
    mesh3_m = std::unique_ptr<Mesh_t>(new Mesh_t(domain3_m));
    layout3_m = std::unique_ptr<FieldLayout_t>(new FieldLayout_t(*mesh3_m, decomp2));

    rho2tr_m.initialize(*mesh3_m, *layout3_m);
    imgrho2tr_m.initialize(*mesh3_m, *layout3_m);
    grntr_m.initialize(*mesh3_m, *layout3_m);

    // helper field for sin
    greentr_m.initialize(*mesh3_m, *layout3_m);

    for(i = 0; i < 3; i++) {
        domain4_m[i] = Index(nr_m[i] + 2);
    }
    mesh4_m = std::unique_ptr<Mesh_t>(new Mesh_t(domain4_m));
    layout4_m = std::unique_ptr<FieldLayout_t>(new FieldLayout_t(*mesh4_m, decomp));

    tmpgreen.initialize(*mesh4_m, *layout4_m);

    // create a domain used to indicate to the FFT's how to construct it's
    // temporary fields.  This is the same as the complex field's domain,
    // but permuted back to the left.
    tmpdomain = layout3_m->getDomain();
    for(i = 0; i < 3; ++i)
        domainFFTConstruct_m[i] = tmpdomain[(i+1) % 3];

    // create the FFT object
    fft_m = std::unique_ptr<FFT_t>(new FFT_t(layout2_m->getDomain(), domainFFTConstruct_m));

    // these are fields that are used for calculating the Green's function.
    // they eliminate some calculation at each time-step.
    for(i = 0; i < 3; ++i) {
        grnIField_m[i].initialize(*mesh2_m, *layout2_m);
        grnIField_m[i][domain2_m] = where(lt(domain2_m[i], nr_m[i]),
                                          domain2_m[i] * domain2_m[i],
                                          (2 * nr_m[i] - domain2_m[i]) *
                                          (2 * nr_m[i] - domain2_m[i]));
    }

    GreensFunctionTimer_m = IpplTimings::getTimer("SF: GreensFTotal");

    if(greensFunction_m == std::string("INTEGRATED")) {
        IntGreensFunctionTimer1_m = IpplTimings::getTimer("SF: IntGreenF1");
        IntGreensFunctionTimer2_m = IpplTimings::getTimer("SF: IntGreenF2");
        IntGreensFunctionTimer3_m = IpplTimings::getTimer("SF: IntGreenF3");
        IntGreensMirrorTimer1_m = IpplTimings::getTimer("SF: MirrorRho1");

        ShIntGreensFunctionTimer1_m = IpplTimings::getTimer("SF: ShIntGreenF1");
        ShIntGreensFunctionTimer2_m = IpplTimings::getTimer("SF: ShIntGreenF2");
        ShIntGreensFunctionTimer3_m = IpplTimings::getTimer("SF: ShIntGreenF3");
        ShIntGreensFunctionTimer4_m = IpplTimings::getTimer("SF: ShIntGreenF4");
        IntGreensMirrorTimer2_m = IpplTimings::getTimer("SF: MirrorRho2");
    } else {
        GreensFunctionTimer1_m = IpplTimings::getTimer("SF: GreenF1");
        GreensFunctionTimer4_m = IpplTimings::getTimer("SF: GreenF4");
    }

    ComputePotential_m = IpplTimings::getTimer("ComputePotential");

#ifdef OPAL_DKS

    if (IpplInfo::DKSEnabled) {
      int dkserr;

      dksbase.setAPI("Cuda", 4);
      dksbase.setDevice("-gpu", 4);
      dksbase.initDevice();

      if (Ippl::myNode() == 0) {

	//create stream for greens function
	dksbase.createStream(streamGreens);
	dksbase.createStream(streamFFT);

	//create fft plans for multiple reuse
	int dimsize[3] = {2*nr_m[0], 2*nr_m[1], 2*nr_m[2]};

	dksbase.setupFFT(3, dimsize);

	//allocate memory
	int sizegreen = tmpgreen.getLayout().getDomain().size();
	int sizerho2_m = rho2_m.getLayout().getDomain().size();
	int sizecomp = grntr_m.getLayout().getDomain().size();

	tmpgreen_ptr = dksbase.allocateMemory<double>(sizegreen, dkserr);
	rho2_m_ptr = dksbase.allocateMemory<double>(sizerho2_m, dkserr);
	rho2real_m_ptr = dksbase.allocateMemory<double>(sizerho2_m, dkserr);

	grntr_m_ptr = dksbase.allocateMemory< std::complex<double>  >(sizecomp, dkserr);
	rho2tr_m_ptr = dksbase.allocateMemory< std::complex<double> > (sizecomp, dkserr);

	//send rho2real_m_ptr to other mpi processes
	//send streamFFT to other processes
	for (int p = 1; p < Ippl::getNodes(); p++) {
	  dksbase.sendPointer( rho2real_m_ptr, p, Ippl::getComm() );
	}
      } else {
	//create stream for FFT data transfer
	dksbase.createStream(streamFFT);
	//receive pointer
	rho2real_m_ptr = dksbase.receivePointer(0, Ippl::getComm(), dkserr);
      }
    }

#endif

}


FFTPoissonSolver::FFTPoissonSolver(PartBunch &beam, std::string greensFunction):
    mesh_m(&beam.getMesh()),
    layout_m(&beam.getFieldLayout()),
    mesh2_m(nullptr),
    layout2_m(nullptr),
    mesh3_m(nullptr),
    layout3_m(nullptr),
    mesh4_m(nullptr),
    layout4_m(nullptr),
    greensFunction_m(greensFunction)
{
    int i;
    domain_m = layout_m->getDomain();

    INFOMSG("FFTPoissonSolver(PartBunch &beam ..." << endl;);



    // For efficiency in the FFT's, we can use a parallel decomposition
    // which can be serial in the first dimension.
    e_dim_tag decomp[3];
    e_dim_tag decomp2[3];
    for(int d = 0; d < 3; ++d) {
        decomp[d] = layout_m->getRequestedDistribution(d);
        decomp2[d] = layout_m->getRequestedDistribution(d);
    }

    // The FFT's require double-sized field sizes in order to (more closely
    // do not understand this ...)
    // simulate an isolated system.  The FFT of the charge density field, rho,
    // would otherwise mimic periodic boundary conditions, i.e. as if there were
    // several beams set a periodic distance apart.  The double-sized fields
    // alleviate this problem.
    for(i = 0; i < 3; i++) {
        hr_m[i] = mesh_m->get_meshSpacing(i);
        nr_m[i] = domain_m[i].length();
        domain2_m[i] = Index(2 * nr_m[i] + 1);
    }
    // create double sized mesh and layout objects for the use in the FFT's
    mesh2_m = std::unique_ptr<Mesh_t>(new Mesh_t(domain2_m));
    layout2_m = std::unique_ptr<FieldLayout_t>(new FieldLayout_t(*mesh2_m, decomp));
#ifdef OPAL_DKS
    rho2_m.initialize(*mesh2_m, *layout2_m, false);
#else
    rho2_m.initialize(*mesh2_m, *layout2_m);
#endif

    NDIndex<3> tmpdomain;
    // Create the domain for the transformed (complex) fields.  Do this by
    // taking the domain from the doubled mesh, permuting it to the right, and
    // setting the 2nd dimension to have n/2 + 1 elements.
    domain3_m[0] = Index(2 * nr_m[3-1] + 1);
    domain3_m[1] = Index(nr_m[0] + 2);

    for(i = 2; i < 3; ++i)
        domain3_m[i] = Index(2 * nr_m[i-1] + 1);

    // create mesh and layout for the new real-to-complex FFT's, for the
    // complex transformed fields
    mesh3_m = std::unique_ptr<Mesh_t>(new Mesh_t(domain3_m));
    layout3_m = std::unique_ptr<FieldLayout_t>(new FieldLayout_t(*mesh3_m, decomp2));
    rho2tr_m.initialize(*mesh3_m, *layout3_m);
    imgrho2tr_m.initialize(*mesh3_m, *layout3_m);
    grntr_m.initialize(*mesh3_m, *layout3_m);

    // helper field for sin
    greentr_m.initialize(*mesh3_m, *layout3_m);

    for(i = 0; i < 3; i++) {
        domain4_m[i] = Index(nr_m[i] + 2);
    }
    mesh4_m = std::unique_ptr<Mesh_t>(new Mesh_t(domain4_m));
    layout4_m = std::unique_ptr<FieldLayout_t>(new FieldLayout_t(*mesh4_m, decomp));

    tmpgreen.initialize(*mesh4_m, *layout4_m);

    // create a domain used to indicate to the FFT's how to construct it's
    // temporary fields.  This is the same as the complex field's domain,
    // but permuted back to the left.
    tmpdomain = layout3_m->getDomain();
    for(i = 0; i < 3; ++i)
        domainFFTConstruct_m[i] = tmpdomain[(i+1) % 3];

    // create the FFT object
    fft_m = std::unique_ptr<FFT_t>(new FFT_t(layout2_m->getDomain(), domainFFTConstruct_m));

    // these are fields that are used for calculating the Green's function.
    // they eliminate some calculation at each time-step.
    for(i = 0; i < 3; ++i) {
        grnIField_m[i].initialize(*mesh2_m, *layout2_m);
        grnIField_m[i][domain2_m] = where(lt(domain2_m[i], nr_m[i]),
                                          domain2_m[i] * domain2_m[i],
                                          (2 * nr_m[i] - domain2_m[i]) *
                                          (2 * nr_m[i] - domain2_m[i]));
    }

#ifdef OPAL_DKS
    if (IpplInfo::DKSEnabled) {
      int dkserr;

      dksbase.setAPI("Cuda", 4);
      dksbase.setDevice("-gpu", 4);
      dksbase.initDevice();

      if (Ippl::myNode() == 0) {

	//create stream for greens function
	dksbase.createStream(streamGreens);
	dksbase.createStream(streamFFT);

	//create fft plans for multiple reuse
	int dimsize[3] = {2*nr_m[0], 2*nr_m[1], 2*nr_m[2]};

	dksbase.setupFFT(3, dimsize);

	//allocate memory
	int sizegreen = tmpgreen.getLayout().getDomain().size();
	int sizerho2_m = rho2_m.getLayout().getDomain().size();
	int sizecomp = grntr_m.getLayout().getDomain().size();

	tmpgreen_ptr = dksbase.allocateMemory<double>(sizegreen, dkserr);
	rho2_m_ptr = dksbase.allocateMemory<double>(sizerho2_m, dkserr);
	rho2real_m_ptr = dksbase.allocateMemory<double>(sizerho2_m, dkserr);

	grntr_m_ptr = dksbase.allocateMemory< std::complex<double>  >(sizecomp, dkserr);
	rho2tr_m_ptr = dksbase.allocateMemory< std::complex<double> > (sizecomp, dkserr);

	//send rho2real_m_ptr to other mpi processes
	//send streamFFT to other processes
	for (int p = 1; p < Ippl::getNodes(); p++) {
	  dksbase.sendPointer( rho2real_m_ptr, p, Ippl::getComm() );
	}
      } else {
	//create stream for FFT data transfer
	dksbase.createStream(streamFFT);
	//receive pointer
	rho2real_m_ptr = dksbase.receivePointer(0, Ippl::getComm(), dkserr);
      }
    }
#endif

}

////////////////////////////////////////////////////////////////////////////
// destructor
FFTPoissonSolver::~FFTPoissonSolver() {
    // delete the FFT object
    //~ delete fft_m;

    // delete the mesh and layout objects
    //~ if(mesh2_m != 0)     delete mesh2_m;
    //~ if(layout2_m != 0)   delete layout2_m;
    //~ if(mesh3_m != 0)     delete mesh3_m;
    //~ if(layout3_m != 0)   delete layout3_m;

#ifdef OPAL_DKS
  //free all the allocated memory
  if (IpplInfo::DKSEnabled) {
    if (Ippl::myNode() == 0) {
      //get number of elements
      int sizegreen = tmpgreen.getLayout().getDomain().size();
      int sizerho2_m = rho2_m.getLayout().getDomain().size();
      int sizecomp = grntr_m.getLayout().getDomain().size();

      //free memory
      dksbase.freeMemory<double>(tmpgreen_ptr, sizegreen);
      dksbase.freeMemory<double>(rho2_m_ptr, sizerho2_m);
      dksbase.freeMemory< std::complex<double> >(grntr_m_ptr, sizecomp);

      //wait for other processes to close handle to rho2real_m_ptr before freeing memory
      MPI_Barrier(Ippl::getComm());
      dksbase.freeMemory<double>(rho2real_m_ptr, sizerho2_m);
      dksbase.freeMemory< std::complex<double> >(rho2tr_m_ptr, sizecomp);
    } else {
      dksbase.closeHandle(rho2real_m_ptr);
      MPI_Barrier(Ippl::getComm());
    }
  }
#endif
}

////////////////////////////////////////////////////////////////////////////
// given a charge-density field rho and a set of mesh spacings hr,
// compute the electric potential from the image charge by solving
// the Poisson's equation

void FFTPoissonSolver::computePotential(Field_t &rho, Vector_t hr, double zshift) {

    // use grid of complex doubled in both dimensions
    // and store rho in lower left quadrant of doubled grid
    rho2_m = 0.0;

    rho2_m[domain_m] = rho[domain_m];

    // needed in greens function
    hr_m = hr;
    // FFT double-sized charge density
    // we do a backward transformation so that we dont have to account for the normalization factor
    // that is used in the forward transformation of the IPPL FFT
    fft_m->transform(-1, rho2_m, rho2tr_m);

    // must be called if the mesh size has changed
    // have to check if we can do G with h = (1,1,1)
    // and rescale later

    // Do image charge.
    // The minus sign is due to image charge.
    // Convolute transformed charge density with shifted green's function.
    IpplTimings::startTimer(GreensFunctionTimer_m);
    shiftedIntGreensFunction(zshift);
    IpplTimings::stopTimer(GreensFunctionTimer_m);

    // Multiply transformed charge density and
    // transformed Green's function. Don't divide
    // by (2*nx_m)*(2*ny_m), as Ryne does; this
    // normalization is done in POOMA's fft routine.
    imgrho2tr_m = - rho2tr_m * grntr_m;

    // Inverse FFT to find image charge potential, rho2_m equals the electrostatic potential.
    fft_m->transform(+1, imgrho2tr_m, rho2_m);

    // Re-use rho to store image potential. Flip z coordinate since this is a mirror image.
    Index I = nr_m[0];
    Index J = nr_m[1];
    Index K = nr_m[2];
    rho[I][J][K] = rho2_m[I][J][nr_m[2] - K - 1];

}


////////////////////////////////////////////////////////////////////////////
// given a charge-density field rho and a set of mesh spacings hr,
// compute the electric field and put in eg by solving the Poisson's equation

void FFTPoissonSolver::computePotential(Field_t &rho, Vector_t hr) {

    IpplTimings::startTimer(ComputePotential_m);

    // use grid of complex doubled in both dimensions
    // and store rho in lower left quadrant of doubled grid
    rho2_m = 0.0;

    rho2_m[domain_m] = rho[domain_m];

    // needed in greens function
    hr_m = hr;

    if (!IpplInfo::DKSEnabled) {
      // FFT double-sized charge density
      // we do a backward transformation so that we dont have to account for the normalization factor
      // that is used in the forward transformation of the IPPL FFT
      fft_m->transform(-1, rho2_m, rho2tr_m);

      // must be called if the mesh size has changed
      // have to check if we can do G with h = (1,1,1)
      // and rescale later
      IpplTimings::startTimer(GreensFunctionTimer_m);
      if(greensFunction_m == std::string("INTEGRATED"))
          integratedGreensFunction();
      else
          greensFunction();
      IpplTimings::stopTimer(GreensFunctionTimer_m);
      // multiply transformed charge density
      // and transformed Green function
      // Don't divide by (2*nx_m)*(2*ny_m), as Ryne does;
      // this normalization is done in POOMA's fft routine.
      rho2tr_m *= grntr_m;

      // inverse FFT, rho2_m equals to the electrostatic potential
      fft_m->transform(+1, rho2tr_m, rho2_m);
      // end convolution
    } else {
#ifdef OPAL_DKS
        dksbase.syncDevice();
        MPI_Barrier(Ippl::getComm());

        if (Ippl::myNode() == 0) {
            IpplTimings::startTimer(GreensFunctionTimer_m);
            integratedGreensFunction();
            IpplTimings::stopTimer(GreensFunctionTimer_m);
            //transform the greens function
            int dimsize[3] = {2*nr_m[0], 2*nr_m[1], 2*nr_m[2]};
            dksbase.callR2CFFT(rho2_m_ptr, grntr_m_ptr, 3, dimsize, streamGreens);
        }
        MPI_Barrier(Ippl::getComm());

        //transform rho2_m keep pointer to GPU memory where results are stored in rho2tr_m_ptr
        fft_m->transformDKSRC(-1, rho2_m, rho2real_m_ptr, rho2tr_m_ptr, dksbase, streamFFT, false);

        if (Ippl::myNode() == 0) {
            //transform the greens function
            //int dimsize[3] = {2*nr_m[0], 2*nr_m[1], 2*nr_m[2]};
            //dksbase.callR2CFFT(rho2_m_ptr, grntr_m_ptr, 3, dimsize, streamGreens);

            //multiply fields and free unneeded memory
            int sizecomp = grntr_m.getLayout().getDomain().size();
            dksbase.syncDevice();
            dksbase.callMultiplyComplexFields(rho2tr_m_ptr, grntr_m_ptr, sizecomp);
        }

        MPI_Barrier(Ippl::getComm());

        //inverse FFT and transfer result back to rho2_m
        fft_m->transformDKSCR(+1, rho2_m, rho2real_m_ptr, rho2tr_m_ptr, dksbase);

        MPI_Barrier(Ippl::getComm());
#endif
    }

    // back to physical grid
    // reuse the charge density field to store the electrostatic potential
    rho[domain_m] = rho2_m[domain_m];

    IpplTimings::stopTimer(ComputePotential_m);
}
///////////////////////////////////////////////////////////////////////////
// calculate the FFT of the Green's function for the given field
void FFTPoissonSolver::greensFunction() {

    //hr_m[0]=hr_m[1]=hr_m[2]=1;

    Vector_t hrsq(hr_m * hr_m);
    IpplTimings::startTimer(GreensFunctionTimer1_m);
    SpecializedGreensFunction<3>::calculate(hrsq, rho2_m, grnIField_m);
    IpplTimings::stopTimer(GreensFunctionTimer1_m);
    // Green's function calculation complete at this point.
    // The next step is to FFT it.
    // FFT of Green's function

    //  ofstream fstr;
    //    fstr.precision(9);
    //    fstr.open("green",ios::out);

    //    for (int i=0;i<domain2_m[0].length();i++)
    //      fstr << i << " " << rho2_m[0][0][i] << endl;
    //    fstr.close();

    //IFF: +1 is forward
    // we do a backward transformation so that we dont have to account for the normalization factor
    // that is used in the forward transformation of the IPPL FFT
    IpplTimings::startTimer(GreensFunctionTimer4_m);
    fft_m->transform(-1, rho2_m, grntr_m);
    IpplTimings::stopTimer(GreensFunctionTimer4_m);

    //    ofstream fstr2;
    //    fstr2.precision(9);
    //    fstr2.open("green_fft",ios::out);
    //    for (int i=0;i<domain2_m[0].length();i++)
    //      fstr2 << i << " " << grntr_m[0][0][i] << endl;
    //    fstr2.close();
}

/** If the beam has a longitudinal size >> transverse size the
 * direct Green function at each mesh point is not efficient
 * (needs a lot of mesh points along the transverse size to
 * get a good resolution)
 *
 * If the charge density function is uniform within each cell
 * the following Green's function can be defined:
 *
 * \f[ \overline{G}(x_i - x_{i'}, y_j - y_{j'}, z_k - z_{k'}  cout << I << endl;
 cout << J << endl;
 cout << K << endl;
 cout << IE << endl;
 cout << JE << endl;
 cout << KE << endl;

 ) = \int_{x_{i'} - h_x/2}^{x_{i'} + h_x/2} dx' \int_{y_{j'} - h_y/2}^{y_{j'} + h_y/2} dy' \int_{z_{k'} - h_z/2}^{z_{k'} + h_z/2} dz' G(x_i - x_{i'}, y_j - y_{j'}, z_k - z_{k'}).
 * \f]
 */
void FFTPoissonSolver::integratedGreensFunction() {



    IpplTimings::startTimer(IntGreensFunctionTimer1_m);

    /**
     * This integral can be calculated analytically in a closed from:
     */


    NDIndex<3> idx =  layout4_m->getLocalNDIndex();
    double cellVolume = hr_m[0] * hr_m[1] * hr_m[2];
    tmpgreen = 0.0;

    for(int k = idx[2].first(); k <= idx[2].last() + 1; k++) {
        for(int j = idx[1].first(); j <=  idx[1].last() + 1; j++) {
            for(int i = idx[0].first(); i <= idx[0].last() + 1; i++) {

                Vector_t vv = Vector_t(0.0);
                vv(0) = i * hr_m[0] - hr_m[0] / 2;
                vv(1) = j * hr_m[1] - hr_m[1] / 2;
                vv(2) = k * hr_m[2] - hr_m[2] / 2;

                double r = sqrt(vv(0) * vv(0) + vv(1) * vv(1) + vv(2) * vv(2));
                double tmpgrn  = -vv(2) * vv(2) * atan(vv(0) * vv(1) / (vv(2) * r)) / 2;
                tmpgrn += -vv(1) * vv(1) * atan(vv(0) * vv(2) / (vv(1) * r)) / 2;
                tmpgrn += -vv(0) * vv(0) * atan(vv(1) * vv(2) / (vv(0) * r)) / 2;
                tmpgrn += vv(1) * vv(2) * log(vv(0) + r);
                tmpgrn += vv(0) * vv(2) * log(vv(1) + r);
                tmpgrn += vv(0) * vv(1) * log(vv(2) + r);

                tmpgreen[i][j][k] = tmpgrn / cellVolume;

            }
        }
    }

    IpplTimings::stopTimer(IntGreensFunctionTimer1_m);

    IpplTimings::startTimer(IntGreensFunctionTimer2_m);

    //assign seems to have problems when we need values that are on another CPU, i.e. [I+1]
    /*assign(rho2_m[I][J][K] ,
      tmpgreen[I+1][J+1][K+1] - tmpgreen[I][J+1][K+1] -
      tmpgreen[I+1][J][K+1] + tmpgreen[I][J][K+1] -
      tmpgreen[I+1][J+1][K] + tmpgreen[I][J+1][K] +
      tmpgreen[I+1][J][K] - tmpgreen[I][J][K]);*/

    Index I = nr_m[0] + 1;
    Index J = nr_m[1] + 1;
    Index K = nr_m[2] + 1;

    // the actual integration
    rho2_m = 0.0;
    rho2_m[I][J][K]  = tmpgreen[I+1][J+1][K+1];
    rho2_m[I][J][K] += tmpgreen[I+1][J][K];
    rho2_m[I][J][K] += tmpgreen[I][J+1][K];
    rho2_m[I][J][K] += tmpgreen[I][J][K+1];
    rho2_m[I][J][K] -= tmpgreen[I+1][J+1][K];
    rho2_m[I][J][K] -= tmpgreen[I+1][J][K+1];
    rho2_m[I][J][K] -= tmpgreen[I][J+1][K+1];
    rho2_m[I][J][K] -= tmpgreen[I][J][K];

    IpplTimings::stopTimer(IntGreensFunctionTimer2_m);

    mirrorRhoField();

    IpplTimings::startTimer(IntGreensFunctionTimer3_m);

    fft_m->transform(-1, rho2_m, grntr_m);

    IpplTimings::stopTimer(IntGreensFunctionTimer3_m);
}

void FFTPoissonSolver::integratedGreensFunctionDKS() {

#ifdef OPAL_DKS
  IpplTimings::startTimer(IntGreensFunctionTimer1_m);
  /**
   * This integral can be calculated analytically in a closed from:
   */
  NDIndex<3> idx =  layout4_m->getDomain();
  dksbase.callGreensIntegral(tmpgreen_ptr, idx[0].length(), idx[1].length(), idx[2].length(),
			     nr_m[0]+1, nr_m[1]+1, hr_m[0], hr_m[1], hr_m[2], streamGreens);

  IpplTimings::stopTimer(IntGreensFunctionTimer1_m);

  IpplTimings::startTimer(IntGreensFunctionTimer2_m);

  Index I = nr_m[0] + 1;
  Index J = nr_m[1] + 1;
  Index K = nr_m[2] + 1;

  dksbase.callGreensIntegration(rho2_m_ptr, tmpgreen_ptr, nr_m[0]+1, nr_m[1]+1, nr_m[2]+1,
				streamGreens);

  IpplTimings::stopTimer(IntGreensFunctionTimer2_m);

  dksbase.callMirrorRhoField(rho2_m_ptr, nr_m[0], nr_m[1], nr_m[2], streamGreens);
#endif

}

void FFTPoissonSolver::shiftedIntGreensFunction(double zshift) {

    tmpgreen = 0.0;
    Field_t grn2(*mesh4_m, *layout4_m);
    grn2 = 0.0;
    NDIndex<3> idx =  layout4_m->getLocalNDIndex();
    double cellVolume = hr_m[0] * hr_m[1] * hr_m[2];
    IpplTimings::startTimer(ShIntGreensFunctionTimer1_m);
    for(int k = idx[2].first(); k <= idx[2].last(); k++) {
        for(int j = idx[1].first(); j <= idx[1].last(); j++) {
            for(int i = idx[0].first(); i <= idx[0].last(); i++) {

                Vector_t vv = Vector_t(0.0);
                vv(0) = i * hr_m[0] - hr_m[0] / 2;
                vv(1) = j * hr_m[1] - hr_m[1] / 2;
                vv(2) = k * hr_m[2] - hr_m[2] / 2 + zshift;

                double r = sqrt(vv(0) * vv(0) + vv(1) * vv(1) + vv(2) * vv(2));
                double tmpgrn  = -vv(2) * vv(2) * atan(vv(0) * vv(1) / (vv(2) * r)) / 2;
                tmpgrn += -vv(1) * vv(1) * atan(vv(0) * vv(2) / (vv(1) * r)) / 2;
                tmpgrn += -vv(0) * vv(0) * atan(vv(1) * vv(2) / (vv(0) * r)) / 2;
                tmpgrn += vv(1) * vv(2) * log(vv(0) + r);
                tmpgrn += vv(0) * vv(2) * log(vv(1) + r);
                tmpgrn += vv(0) * vv(1) * log(vv(2) + r);

                tmpgreen[i][j][k] = tmpgrn / cellVolume;

            }
        }
    }
    IpplTimings::stopTimer(ShIntGreensFunctionTimer1_m);

    IpplTimings::startTimer(ShIntGreensFunctionTimer2_m);
    for(int k = idx[2].first(); k <= idx[2].last(); k++) {
        for(int j = idx[1].first(); j <= idx[1].last(); j++) {
            for(int i = idx[0].first(); i <= idx[0].last(); i++) {

                Vector_t vv = Vector_t(0.0);
                vv(0) = i * hr_m[0] - hr_m[0] / 2;
                vv(1) = j * hr_m[1] - hr_m[1] / 2;
                vv(2) = k * hr_m[2] - hr_m[2] / 2 + zshift - nr_m[2] * hr_m[2];

                double r = sqrt(vv(0) * vv(0) + vv(1) * vv(1) + vv(2) * vv(2));
                double tmpgrn  = -vv(2) * vv(2) * atan(vv(0) * vv(1) / (vv(2) * r)) / 2;
                tmpgrn += -vv(1) * vv(1) * atan(vv(0) * vv(2) / (vv(1) * r)) / 2;
                tmpgrn += -vv(0) * vv(0) * atan(vv(1) * vv(2) / (vv(0) * r)) / 2;
                tmpgrn += vv(1) * vv(2) * log(vv(0) + r);
                tmpgrn += vv(0) * vv(2) * log(vv(1) + r);
                tmpgrn += vv(0) * vv(1) * log(vv(2) + r);

                grn2[i][j][k] = tmpgrn / cellVolume;

            }
        }
    }
    IpplTimings::stopTimer(ShIntGreensFunctionTimer2_m);

    /**
     ** (x[0:nr_m[0]-1]^2 + y[0:nr_m[1]-1]^2 + (z_c + z[0:nr_m[2]-1])^2)^{-0.5}
     ** (x[nr_m[0]:1]^2   + y[0:nr_m[1]-1]^2 + (z_c + z[0:nr_m[2]-1])^2)^{-0.5}
     ** (x[0:nr_m[0]-1]^2 + y[nr_m[1]:1]^2   + (z_c + z[0:nr_m[2]-1])^2)^{-0.5}
     ** (x[nr_m[0]:1]^2   + y[nr_m[1]:1]^2   + (z_c + z[0:nr_m[2]-1])^2)^{-0.5}
     **
     ** (x[0:nr_m[0]-1]^2 + y[0:nr_m[1]-1]^2 + (z_c - z[nr_m[2]:1])^2)^{-0.5}
     ** (x[nr_m[0]:1]^2   + y[0:nr_m[1]-1]^2 + (z_c - z[nr_m[2]:1])^2)^{-0.5}
     ** (x[0:nr_m[0]-1]^2 + y[nr_m[1]:1]^2   + (z_c - z[nr_m[2]:1])^2)^{-0.5}
     ** (x[nr_m[0]:1]^2   + y[nr_m[1]:1]^2   + (z_c - z[nr_m[2]:1])^2)^{-0.5}
     */

    IpplTimings::startTimer(ShIntGreensFunctionTimer3_m);

    Index I = nr_m[0] + 1;
    Index J = nr_m[1] + 1;
    Index K = nr_m[2] + 1;

    // the actual integration
    rho2_m = 0.0;
    rho2_m[I][J][K]  = tmpgreen[I+1][J+1][K+1];
    rho2_m[I][J][K] += tmpgreen[I+1][J][K];
    rho2_m[I][J][K] += tmpgreen[I][J+1][K];
    rho2_m[I][J][K] += tmpgreen[I][J][K+1];
    rho2_m[I][J][K] -= tmpgreen[I+1][J+1][K];
    rho2_m[I][J][K] -= tmpgreen[I+1][J][K+1];
    rho2_m[I][J][K] -= tmpgreen[I][J+1][K+1];
    rho2_m[I][J][K] -= tmpgreen[I][J][K];

    tmpgreen = 0.0;
    tmpgreen[I][J][K]  = grn2[I+1][J+1][K+1];
    tmpgreen[I][J][K] += grn2[I+1][J][K];
    tmpgreen[I][J][K] += grn2[I][J+1][K];
    tmpgreen[I][J][K] += grn2[I][J][K+1];
    tmpgreen[I][J][K] -= grn2[I+1][J+1][K];
    tmpgreen[I][J][K] -= grn2[I+1][J][K+1];
    tmpgreen[I][J][K] -= grn2[I][J+1][K+1];
    tmpgreen[I][J][K] -= grn2[I][J][K];
    IpplTimings::stopTimer(ShIntGreensFunctionTimer3_m);

    mirrorRhoField(tmpgreen);

    IpplTimings::startTimer(ShIntGreensFunctionTimer4_m);
    fft_m->transform(-1, rho2_m, grntr_m);
    IpplTimings::stopTimer(ShIntGreensFunctionTimer4_m);

}

void FFTPoissonSolver::mirrorRhoField() {

    IpplTimings::startTimer(IntGreensMirrorTimer1_m);
    Index aI(0, 2 * nr_m[0]);
    Index aJ(0, 2 * nr_m[1]);

    Index J(0, nr_m[1]);
    Index K(0, nr_m[2]);

    Index IE(nr_m[0] + 1, 2 * nr_m[0]);
    Index JE(nr_m[1] + 1, 2 * nr_m[1]);
    Index KE(nr_m[2] + 1, 2 * nr_m[2]);

    Index mirroredIE = 2 * nr_m[0] - IE;
    Index mirroredJE = 2 * nr_m[1] - JE;
    Index mirroredKE = 2 * nr_m[2] - KE;

    rho2_m[0][0][0] = rho2_m[0][0][1];

    rho2_m[IE][J ][K ] = rho2_m[mirroredIE][J         ][K         ];
    rho2_m[aI][JE][K ] = rho2_m[aI        ][mirroredJE][K         ];
    rho2_m[aI][aJ][KE] = rho2_m[aI        ][aJ        ][mirroredKE];

    IpplTimings::stopTimer(IntGreensMirrorTimer1_m);
}

void FFTPoissonSolver::mirrorRhoField(Field_t & ggrn2) {

    IpplTimings::startTimer(IntGreensMirrorTimer2_m);
    Index aI(0, 2 * nr_m[0]);
    Index aK(0, 2 * nr_m[2]);

    Index I(0, nr_m[0]);
    Index J(0, nr_m[1]);

    Index IE(nr_m[0] + 1, 2 * nr_m[0]);
    Index JE(nr_m[1] + 1, 2 * nr_m[1]);
    Index KE(nr_m[2] + 1, 2 * nr_m[2]);

    Index mirroredIE = 2*nr_m[0] - IE;
    Index mirroredJE = 2*nr_m[1] - JE;
    Index shiftedKE  = KE - nr_m[2];

    rho2_m[I ][J ][KE] = ggrn2[I          ][J         ][shiftedKE];
    rho2_m[IE][J ][aK] = rho2_m[mirroredIE][J         ][aK       ];
    rho2_m[aI][JE][aK] = rho2_m[aI        ][mirroredJE][aK       ];

    IpplTimings::stopTimer(IntGreensMirrorTimer2_m);
}

Inform &FFTPoissonSolver::print(Inform &os) const {
    os << "* ************* F F T P o i s s o n S o l v e r ************************************ " << endl;
    os << "* h " << hr_m << '\n';
    os << "* ********************************************************************************** " << endl;
    return os;
}

/***************************************************************************
 * $RCSfile: FFTPoissonSolver.cc,v $   $Author: adelmann $
 * $Revision: 1.6 $   $Date: 2001/08/16 09:36:08 $
 ***************************************************************************/