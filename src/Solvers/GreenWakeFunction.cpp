#include "Solvers/GreenWakeFunction.hh"
#include "Physics/Physics.h"
//#ifdef ENABLE_WAKE_TESTS
#include "TestLambda.h"	// used for tests
//#endif

#include <fstream>
#include <string>
#include <vector>
#include <istream>
#include <iostream>  // Needed for stream I/O
#include <fstream>   // Needed for file I/O
#include <iomanip>   // Needed for I/O manipulators

#include "Utilities/OpalException.h"

//IFF: TEST
//#define ENABLE_WAKE_DEBUG
//#define ENABLE_WAKE_DUMP
//#define ENABLE_WAKE_TESTS

//#define readWakeFromFile
//#define WakeFile "test.sdds"


/**
 * @brief	just a Testfunction!  Calculate the energy of the Wakefunction with the lambda
 *
 *
 * @todo		In this code one can only apply either the longitudinal wakefield or the transversal wakefield. One should implement that both wakefields can be applied to the particle beam 
 * @todo		NBins must be set equal to MT of the fieldsolver. This should be changed. (the length of lineDensity_m must be NBins and not to MT)
 *
 * @param[in] 	ref
 * @param[in]	NBIN number of Bins
 * @param[in]	Z0 impedance of the tube
 * @param[in]	radius radius of the tube
 * @param[in]	sigma material constant 
 * @param[in]	c speed of light
 * @param[in]	acMode 1 for AC and 2 for DC
 * @param[in]	tau material constant
 * @param[in]	direction 0 for transversal and 1 for Longitudinal
 * @param[in]	constLength  true if the length of the particle bunch is considered as constant
 * @param[in]   fname read wake from file
 */
GreenWakeFunction::GreenWakeFunction(PartData &ref,  vector<OpalFilter*> filters, int NBIN, double Z0, double radius, double sigma, double c, int acMode, double tau, int direction, bool constLength, string fname):
    WakeFunction(ref),
    filters_m(filters.begin(),filters.end()),
    lineDensity_m(),
    Z0(Z0),
    radius(radius),
    sigma(sigma),
    c(c),
    acMode(acMode),
    tau(tau),
    direction(direction),
    constLength(constLength),
    filename_m(fname),
    NBin(NBIN)
{
#ifdef ENABLE_WAKE_DEBUG
    *gmsg << "* ************* W A K E ************************************************************ " << endl;
    *gmsg << "* Entered GreenWakeFunction::GreenWakeFunction " << '\n';
    *gmsg << "* ********************************************************************************** " << endl;
#endif
    FftWField = 0;
}

/**
 * @brief	given a vector of length N, distribute the indexes among the available processors
 *
 *
 * @todo		make this function general available 
 *
 *
 * @param[in] 	length of vector
 * @param[out]	first: lowIndex, second: hiIndex
 */
pair<int,int> GreenWakeFunction::distrIndices(int vectLen) {

    pair<int,int> dist;

    //IFF: properly distribute remainder
    int rem = vectLen - (vectLen / Ippl::getNodes()) * Ippl::getNodes();
    int tmp = (rem > Ippl::myNode()) ? 1 : 0;
    int locBunchRange = vectLen / Ippl::getNodes() + tmp;
    
    dist.first = locBunchRange * Ippl::myNode() + (1-tmp)*rem;
    dist.second = dist.first + locBunchRange - 1;
    
    return dist;
}

void GreenWakeFunction::apply(PartBunch &bunch)
{
#ifdef ENABLE_WAKE_TESTS
    // overwrite the line density
    testApply(bunch);
#else

    double spacing, mindist;
    Vector_t rmin, rmax;

    bunch.calcBeamParameters();
    bunch.get_bounds(rmin,rmax);
    //FIXME IFF: do we have unitless r's here? is that what we want?

    mindist = rmin(2);
    switch(direction) {
        case LONGITUDINAL:  spacing = abs(rmax(2) - rmin(2)); break; //FIXME: Kann mann das Spacing immer ändern?
        case TRANSVERSAL:   spacing = rmax(0)*rmax(0)+rmax(1)*rmax(1); break;
        default:            throw OpalException("GreenWakeFunction","invalid direction specified");
    }
    assert (NBin > 0);
    spacing /= NBin;
    //spacing /= (NBin-1); //FIXME: why -1?

    // Calculate the Wakefield if needed
    if(FftWField == 0) {
        FftWField = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (2*NBin-1));
        if(filename_m != "")
            setWakeFromFile(NBin, spacing);
        else
            CalcWakeFFT(Z0, radius, sigma, c, acMode, tau, direction, spacing, NBin, FftWField);
    } else if(!constLength)
        CalcWakeFFT(Z0, radius, sigma, c, acMode, tau, direction, spacing, NBin, FftWField);

    // Calculate the line density of the particle bunch
    bunch.boundp(); //IFF: added boundp (need if no FS available)
    bunch.calcLineDensity();
    bunch.getLineDensity(lineDensity_m);

#ifdef ENABLE_WAKE_DEBUG
    *gmsg << "* ************* W A K E ************************************************************ " << endl;
    *gmsg << "* GreenWakeFunction::apply  lineDensity_m.size() = "<< lineDensity_m.size() << '\n';
    *gmsg << "* ********************************************************************************** " << endl;
#endif

    // smooth the line density of the particle bunch
    for (vector<OpalFilter*>::const_iterator fit = filters_m.begin(); fit != filters_m.end(); ++fit)
        (*fit)->apply(lineDensity_m);

    // set some constants
    double charge;
    double K = 0;	// constant to normalize the lineDensity_m to 1
    for(int i=0; i<lineDensity_m.size(); i++) 
        K += lineDensity_m[i];
    K = 1/K;
    double OutEnergy[NBin];
    charge = bunch.Q[1]; //FIXME: why 1? bunch,getTotalCharge()?

    // compute the kick due to the wakefield
    compEnergy(K, charge, NBin, FftWField, lineDensity_m, OutEnergy);

    // Add the right OutEnergy[i] to all the particles
    //FIXME: can we specify LONG AND TRANS?
    switch(direction) {
        case LONGITUDINAL: 
            for(int i=0; i<bunch.getLocalNum(); i++) {

                //FIXME: Stimmt das????????? (von den einheiten) 
                // calculate bin containing particle
                int idx = (int)(floor((bunch.R[i](2)-mindist)/spacing));
                //IFF: should be ok
                if(idx==NBin) idx--;
                assert (idx >= 0 && idx < NBin);
                double dE = OutEnergy[idx]; 
                bunch.Ef[i](2) += dE;

            }
            break;

        case TRANSVERSAL:
            for(int i=0; i<bunch.getLocalNum(); i++) {

                // calculate bin containing particle
                int idx = (int)(floor((bunch.R[i](2)-mindist)/spacing)); 
                //IFF: should be ok
                if(idx==NBin) idx--;
                assert (idx >= 0 && idx < NBin);
                double dE = OutEnergy[idx]; 
                
                // ACHTUNG spacing auch in transversal richtung
                double dist = sqrt( bunch.R[i](0)*bunch.R[i](0) + bunch.R[i](1)*bunch.R[i](1) );
                assert (dist > 0);
                bunch.Ef[i](0) += dE * bunch.R[i](0)/dist;
                bunch.Ef[i](1) += dE * bunch.R[i](1)/dist;

            }
            break;

        default: 
            throw OpalException("GreenWakeFunction","invalid direction specified");
    }

#ifdef ENABLE_WAKE_DUMP
    ofstream  f2("OutEnergy.dat"); 
    f2 << "# Energy of the Wake calculated in Opal" << endl;
    f2 << "# Z0 = " << Z0 << endl;
    f2 << "# radius = " << radius << endl;
    f2 << "# sigma = " << sigma << endl;
    f2 << "# c = " << c << endl;
    f2 << "# acMode = " << acMode << endl;
    f2 << "# tau = " << tau << endl;
    f2 << "# direction = " << direction << endl;
    f2 << "# spacing = " << spacing << endl;
    f2 << "# Lbunch = " << NBin << endl;
    for (int i=0; i<NBin; i++){
        f2 << i+1 << " " << OutEnergy[i] << endl; 
    }
    //for (int i=0; i<lineDensity_m.size(); i++){
    //    f2 << i << " " << lineDensity_m[i] << endl; 
    //}
    (f2).close(); 
#endif

#endif //ENABLE_WAKE_TESTS
}

/**
 * @brief	Just a test function
 */
void GreenWakeFunction::testApply(PartBunch &bunch)
{
    double spacing, mindist;
    mindist = 0;
    spacing = 1e-6; //IFF: charge in testLambda.h in 1um spacings
    NBin = 294;

    if (FftWField==0) {
        FftWField = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (2*NBin-1));
        CalcWakeFFT(Z0, radius, sigma, c, acMode, tau, direction, spacing, NBin, FftWField);
    } else if(constLength == false) {
        CalcWakeFFT(Z0, radius, sigma, c, acMode, tau, direction, spacing, NBin, FftWField);
    }
    // determine K and charge
    double charge;
    double K = 0.20536314319923724e-9; //K normalizes nC data in lambda.h?
    //double K = 0.20536314319923724; 
    double OutEnergy[NBin];
    //charge = 0.8e-10; // pC //IFF: WHY?
    charge = 0.8e-9; // nC

    compEnergy(K, charge, NBin, FftWField, testLambda, OutEnergy);

#ifdef ENABLE_WAKE_TESTS
    ofstream  f2("OutEnergy.dat"); 
    f2 << "# Energy of the Wake calculated in Opal" << endl;
    f2 << "# Z0 = " << Z0<<endl;
    f2 << "# radius = " << radius<<endl;
    f2 << "# sigma = " << sigma<<endl;
    f2 << "# c = " << c<<endl;
    f2 << "# acMode = " << acMode<<endl;
    f2 << "# tau = " << tau<<endl;
    f2 << "# direction = " << direction<<endl;
    f2 << "# spacing = " << spacing<<endl;
    f2 << "# Lbunch = " << NBin<<endl;
    for (int i=0;i<NBin;i++){
        f2 << i+1 <<" "<< OutEnergy[i] << endl; 
    }
    (f2).close(); 
#endif

}

GreenWakeFunction::~GreenWakeFunction()
{
    fftw_free(FftWField);
}

/**
 * @brief	just a Testfunction!  Calculate the energy of the Wakefunction with the lambda
 *
 *
 * @param[in] 	K a constant
 * @param[in] 	charge a constant
 * @param[in] 	LengthWake length of the Wake function
 * @param[in] 	FftWField FFT from the zero padded wake function
 * @param[in] 	lambda the distribution of the Particles
 * @param[out]	OutEnergy this is the Output
 */
void GreenWakeFunction::compEnergy(const double K, const double charge, const int LengthWake, const fftw_complex* FftWField, const double* lambda, double* OutEnergy)
{
    fftw_plan p;
    int N = 2*LengthWake -1;
    // Allocate Space for the Result
    double *energy = (double*)fftw_malloc(N*sizeof(double));

    // Allocate Space for the zero padded lambda and its Fourier Transformed
    double *pLambda = (double*)fftw_malloc(N*sizeof(double));
    fftw_complex* FftLambda = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    fftw_complex* final = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)* N);

    // fill the arrays with data
    for (int i=0;i<LengthWake;i++){
        pLambda[i] = lambda[i];
    }
    // make the Zero padding
    for (int i=LengthWake;i<N;i++){
        pLambda[i] = 0;
    }

    //FFT of the lambda
    p =  fftw_plan_dft_r2c_1d(N, pLambda, FftLambda, FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);	

    // convolution -> multiplication in Fourier space
    for (int i=0; i<N; i++){		
        final[i][0] = (FftWField[i][0]*FftLambda[i][0]-FftWField[i][1]*FftLambda[i][1]);
        final[i][1] = (FftWField[i][0]*FftLambda[i][1]+FftWField[i][1]*FftLambda[i][0]);
    }

    // inverse transform to get c, the convolution of a and b;
    //NOTE: this has the side effect of overwriting C 
    p = fftw_plan_dft_c2r_1d(N, final, energy,FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);


    // Write the result to the output:
    for (int i=0; i<LengthWake; i++){
        OutEnergy[i] = -(charge) * K * energy[i] ;
        OutEnergy[i] /= (2*LengthWake);
    }	

    // Free the Memory
    fftw_free(energy);
}

/**
 * @brief	Calculate the energy of the Wakefunction with the lambda
 *
 *
 * @param[in] 	K a constant
 * @param[in] 	charge a constant
 * @param[in] 	LengthWake length of the Wake function
 * @param[in] 	FftWField FFT with zero padding of the wake function
 * @param[in] 	lambda the distribution of the Particles
 * @param[out]	OutEnergy this is the Output
 *
 */
void GreenWakeFunction::compEnergy(const double K, const double charge, const int LengthWake, const fftw_complex* FftWField, vector<double> lambda, double* OutEnergy)
{
    int N = 2*LengthWake -1;
    // Allocate Space for the Result
    double *energy = (double*)fftw_malloc(N*sizeof(double));

    // Allocate Space for the zero padded lambda and its Fourier Transformed
    double *pLambda = (double*)fftw_malloc(N*sizeof(double));
    fftw_complex* FftLambda = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex* final = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)* N);

    // fill the arrays with data
    for (int i=0;i<LengthWake;i++){
        pLambda[i] = lambda[i];
    }
    // make the Zero padding
    for (int i=LengthWake;i<N;i++){
        pLambda[i] = 0;
    }

    //FFT of the lambda
    fftw_plan p;
    p =  fftw_plan_dft_r2c_1d(N, pLambda, FftLambda, FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);

    // Convolution -> just a multiplication in Fourier space
    for (int i=0; i<N; i++){		
        final[i][0] = (FftWField[i][0]*FftLambda[i][0]-FftWField[i][1]*FftLambda[i][1]);
        final[i][1] = (FftWField[i][0]*FftLambda[i][1]+FftWField[i][1]*FftLambda[i][0]);
    }

    // IFFT
    p = fftw_plan_dft_c2r_1d(N, final, energy,FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);

    // Write the result to the output:
    for (int i=0; i<LengthWake; i++){
        OutEnergy[i] = -(charge ) * K * energy[i] ;
        OutEnergy[i] /= (2*LengthWake);
    }	

    // Free the Memory
    fftw_free(energy);
}

/**
 * @brief 	Calculate the FFT of the Wakefunction 
 *
 * @param[in]	Z0 impedance of the tube
 * @param[in]	radius radius of the tube
 * @param[in]	sigma material constant 
 * @param[in]	c speed of light
 * @param[in]	mode 1 for AC and 2 for DC
 * @param[in]	tau material constant
 * @param[in]	direction 0 for transversal and 1 for Longitudinal
 * @param[in]	spacing distance between 2 slice in the line distribution
 * @param[in]	Lbunch length of the Wake function
 * @param[in] 	FftWField FFT from the zero padded wake function
 *
 */
void GreenWakeFunction::CalcWakeFFT(double Z0, double radius, double sigma, double c, double mode, double tau, double direction, double spacing, int Lbunch, fftw_complex *FftWField)
{
    // Set integration properties
    double  a, b;
    unsigned int N;
    a=1;
    b=1000000;
    N=1000000;

    double wake[Lbunch];	

    // Allocate Memory for the FFT
    double *wfield;
    wfield = (double*) fftw_malloc(sizeof(double) * (2*Lbunch-1));

    ofstream file;

#ifdef ENABLE_WAKE_TESTS
    if (Ippl::myNode() == 0) {
        file.open("wake.dat"); 
        file << "# Wake calculated in Opal" << endl;
        file << "# Z0 = " << Z0 << endl;
        file << "# radius = " << radius << endl;
        file << "# sigma = " << sigma << endl;
        file << "# c = " << c << endl;
        file << "# mode = " << mode << endl;
        file << "# tau = " << tau << endl;
        file << "# direction = " << direction << endl;
        file << "# spacing = " << spacing << endl;  
        file << "# Lbunch = " << Lbunch << endl;
    }
#endif

    const pair<int,int> myDist = distrIndices(Lbunch);
    const int lowIndex = myDist.first;
    const int hiIndex  = myDist.second;

    for (int i=0;i<Lbunch;i++)
        wfield[i] = 0.0;

    /** 
      Calculate the Wakefield on all processors	 
      */

    if (Ippl::myNode() != Ippl::getNodes()-1) {
        for (int i=lowIndex; i <= hiIndex; i++) {    
            Wake w(i*spacing, Z0, radius, sigma, c, mode, tau, direction);
            wfield[i] = simpson(w,a,b,N);
        }
    } 
    else {
        //IFF: changed  to <= with new distr
        for (int i=lowIndex; i <= hiIndex; i++) {    
            Wake w(i*spacing, Z0, radius, sigma, c, mode, tau, direction);
            wfield[i] = simpson(w,a,b,N);
        }
    }
    /**
      Reduce the results
      */
    reduce(&(wfield[0]),&(wfield[0]) + Lbunch, &(wfield[0]), OpAddAssign());


#ifdef ENABLE_WAKE_TESTS
    if (Ippl::myNode() == 0) {
        for (int i=0;i<Lbunch;i++) {
            file << i+1 << "   " << wfield[i] << endl;
        }
        (file).close(); 
    }
#endif

    // make the Zero padding
    for (int i=Lbunch;i<(2*Lbunch-1);i++){
        wfield[i] = 0;
    }

    // calculate the FFT of the Wakefield
    fftw_plan p;
    p = fftw_plan_dft_r2c_1d((2*Lbunch-1), wfield, FftWField, FFTW_ESTIMATE);
    fftw_execute(p); 
    fftw_destroy_plan(p);

#ifdef ENABLE_WAKE_TESTS_FFT_OUT
    ofstream  f2("FFTwake.dat"); 
    f2 << "# FFT of the Wake calculated in Opal" << endl;
    f2 << "# Z0 = " << Z0<<endl;
    f2 << "# radius = " << radius<<endl;
    f2 << "# sigma = " << sigma<<endl;
    f2 << "# c = " << c<<endl;
    f2 << "# mode = " << mode<<endl;
    f2 << "# tau = " << tau<<endl;
    f2 << "# direction = " << direction<<endl;
    f2 << "# spacing = " << spacing<<endl;
    f2 << "# Lbunch = " << Lbunch<<endl;
    for (int i=0;i<Lbunch;i++){
        f2 << i+1 <<": "<< FftWField[i][0] << " + " << FftWField[i][1] << " i" << endl; 
    }
    (f2).close(); 
#endif
    // Free the memory
    fftw_free(wfield);

}

/**
 * @brief	reads in the wakefield from file
 */
void GreenWakeFunction::setWakeFromFile(int NBin, double spacing)
{
    Inform msg ("Read sdds wake function ");
    std::ifstream fs;
    fs.open(filename_m.c_str());

    if(fs.fail()) {
        throw OpalException("GreenWakeFunction::setWake",
                "Open file operation failed, please check if \""
                + filename_m +  "\" really exists.");

        msg << "readSDDS1 "<< "Open file operation failed, please check if " << filename_m <<  " really exists." << endl;
        return;
    }

    string name;
    fs >> name;
    msg << "readSDDS" << " SSDS1 read = " <<name<< endl;
    if(name.compare("SDDS1") != 0) {
        throw OpalException("GreenWakeFunction::setWake",
                " No SDDS1 File. A SDDS1 file should start with a SDDS1 String. Check file \""
                + filename_m +  "\" ");

        msg << "readSDDS1 "<< "No SDDS1 File. A SDDS1 file should start with a SDDS1 String. Check file \"" << filename_m <<  "\" "<< endl;
        return;
    }

    char temp[256];
    for (int i=0; i<6; i++) {
        fs.getline(temp,256);
        msg << "readSDDS" << "line " << i << " :   " << temp << endl;
    }

    int Np;
    fs >> Np;
    msg << "readSDDS" << " header read" << endl;
    if ( Np <= 0 ) {
        throw OpalException("GreenWakeFunction::setWake",
                " The particle number should be bigger than zero! Please check the first line of file \""
                + filename_m +  "\".");
        msg << "readSDDS" << " The particle number should be bigger than zero! Please check the first line of file " << filename_m << endl;
    }

    msg  << "readSDDS" << " Np = " <<Np << endl;
    double wake[Np], dist[Np], t[Np];

    // read the wakefunction
    for(unsigned int i=0; i<Np; i++) {
        if( !fs.eof() ) {
            fs >> dist[i] >> wake[i] >> t[i];
        }
        if(fs.eof()) {
            throw OpalException("GreenWakeFunction::setWake",
                    " End of file reached before the whole wakefield is imported, please check file \""
                    + filename_m +  "\".");

            msg << "readSDDS" <<" End of file reached before the whole wakefield is imported, please check file " << filename_m << endl;
            return;
        }
    }
    // if needed interpolate the wake in a way that the wake form the file fits to the wake needs in the code (??)
    double wakefield[NBin];
    int j;
    for(int i=0; i<NBin; i++) {
        j=0;
        while(dist[j]<i*spacing)
            j++; 
        // linear interpolation
        wakefield[i] = wake[j] + ((wake[j+1] - wake[j]) / (dist[j+1] - dist[j]) * (i*spacing - dist[j]));
    }
    setWake(wakefield, NBin);
}



/**
 * @brief	Set the new Wakefield
 *
 * set constLength to true, so that this wakefield is not overwritten afterwards
 *
 * @param[in]	wakefield array with the data of the new Wakefield
 * @param[in]	N length of the wakefield
 *
 */
void GreenWakeFunction::setWake(double* wakefield, double N)
{
    //FIXME: this is an evil side-effect. WHY?
    constLength = true;
    NBin = N;
    // allocate new space for the wakefield
    if (FftWField != 0)
        fftw_free(FftWField);
    FftWField = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (2*NBin-1));


    // Allocate Memory for the FFT
    double *wfield;
    wfield = (double*) fftw_malloc(sizeof(double) * (2*NBin-1));

    // fill the wake with the new wake data
    for (int i=0;i<NBin;i++)
        wfield[i] = wakefield[i];
    
    // make the Zero padding
    for (int i=NBin;i<(2*NBin-1);i++)
        wfield[i] = 0;

    // calculate the FFT of the Wakefield
    fftw_plan p; 
    p = fftw_plan_dft_r2c_1d((2*NBin-1), wfield, FftWField, FFTW_ESTIMATE);
    fftw_execute(p); 
    fftw_destroy_plan(p);

    // Free the memory
    fftw_free(wfield);
}

const string GreenWakeFunction::getType() const
{
    return "GreenWakeFunction";
}

