#include "wake.h"

using namespace std;

#define pi = M_PI;



class Wake {

	public:
	/**
	* @brief	Constructor
	*/
	Wake() 
	: Z0_(120. * M_PI), a_(5.), sigma_(6.45337 * 10000000), c_(299792458), s_(0), acMode_(1), tau_(2.70187 * 0.00000000000001), direction_(1)
	{}
	
	/**
	* @brief	Constructor
	*
	* @param[in]	s position in the wake
	*
	*/
    	Wake(double s) 
	: Z0_(120. * M_PI), a_(5.), sigma_(6.45337 * 10000000), c_(299792458), s_(s), acMode_(1), tau_(2.70187 * 0.00000000000001), direction_(1)
	{}
	
	/**
	* @brief	Constructor
	*
	* @param[in]	s position in the wake
	* @param[in]	Z0 impedance of the beam pipe
	* @param[in]	sigma material constant 
	* @param[in]	c speed of light
	* @param[in]	acMode 1 for AC and 2 for DC
	* @param[in]	tau material constant
	* @param[in]	direction 0 for transversal and 1 for Longitudinal
	*/
	Wake(double s, double Z0, double a, double sigma, double c, double acMode, double tau, double direction) 
	: Z0_(Z0), a_(a), sigma_(sigma), c_(c), s_(s), acMode_(acMode), tau_(tau), direction_(direction)
	{}
	
	/**
	* @brief	Constructor
	*
	* @param[in]	material can be either "Copper" or "Alu"	
	* @param[in]	direction can be either "Transversal" or "Longitudinal"
	* @param[in]	Mode can be either "AC" or "DC"
	* @param[in]	a radius of the tube
	* @param[in]	s position in the wake
	*
	*/	
	Wake(string material, string direction, string Mode, double a, double s)
	{
		a_ = a;
		Z0_ = 120. * M_PI;
		c_ = 299792458;
		s_ = s;
		if (material.compare("Copper") ==0)
		{
			sigma_ = 6.45337 * 10000000;
			tau_ = 2.70187 * 0.00000000000001;
		}
		else if (material.compare("Alu") ==0)
		{
			sigma_ = 4.22807 * 10000000;
			tau_ = 8.00554 * 0.000000000000001;
		}
		if (direction.compare("Transversal") == 0)
		{
			direction_ = 0;
		}
		else if (direction.compare("Longitudinal") == 0)
		{
			direction_ = 1;
		}
		if (Mode.compare("AC")==0)
		{
			acMode_ = 1;
		}
		else if (Mode.compare("DC")==0)
		{
			acMode_ = 2;
		}
	}
	
	/**
	* @brief	Used to integrate the function 
	*
	* @param[in]	k parameter
	*
	* @return	the function value at position k	
	*/	
	double operator()(double k) {
		
		complex <double> i(0, 1);
		complex <double>  Z(0,0);
		double signK;
		signK = (k>0 ? 1 : -1);
		
	    if (acMode_ ==1){	// AC 
			Z = ( Z0_/(2 * M_PI*0.001*a_))* 1.0/(sqrt( Z0_*abs(k)/2)* sqrt(sigma_ /(1.0-i*c_*k*tau_)) *(i+signK)/k -(i*k*0.001*a_)/2.0);
	    }
	    if(acMode_ ==2){	//DC
	    	Z = ( Z0_/(2 * M_PI*0.001*a_))* 1.0/(sqrt(sigma_* Z0_*abs(k)/2) *(i+signK)/k -(i*k*0.001*a_)/2.0);
	    }
	    if (direction_==1){	// Longitudinal	
			return real(Z) *cos(k*s_) *2.0 *c_ / M_PI * 0.000000000001;  
	    }
	    else{	//Transversal
	    	return real(Z)*c_/k *cos(k*s_) *2.0 *c_ / M_PI * 0.000000000001; 
	    }
	}
	
	private: 
	/// impedanz
	double Z0_;
	/// radius
	double a_;
	/// material constant
	double sigma_;
	/// speed of light 
	double c_;
	/// distanc from the particle
	double s_;
	/// conductivity either "AC" or "DC"
	double acMode_;
	/// material constant
	double tau_;
	/// direction either "Transversal" or "Longitudinal"
	double direction_;
	
};



/**
* @brief	Calculate meany diffrent Wakefield and the FFT's those Wakfields
*
*/
void calcWakeSerie(){

	cout.precision(30);
	
	// Set integration propartys
  	double  a, b;
  	unsigned int N;
	a=1;
	b=1000000;
	N=1000000;
    	int Lbunch = 256;
    	
    	double radius = 5.0;
    	string material ("Copper");
    	string direction("Longitudinal");
    	string mode ("AC");
    	ofstream W1( "wake_Lo_Circ_Cu_AC_5_CPP.dat"); 
    	ofstream FFT1( "wakeFFT_Lo_Circ_Cu_AC_5_CPP.dat"); 
 	CalcWakeFFT(&W1, &FFT1, Lbunch, material, direction, mode, radius) ; 
 	
 	
 	material.erase();
 	material.insert(0,"Copper");
     	direction.erase();
   	direction.insert(0,"Longitudinal");
    	mode.erase();
    	mode.insert(0,"DC");
    	
    	ofstream W2( "wake_Lo_Circ_Cu_DC_5_CPP.dat"); 
    	ofstream FFT2( "wakeFFT_Lo_Circ_Cu_DC_5_CPP.dat"); 
 	CalcWakeFFT(&W2, &FFT2, Lbunch, material, direction, mode, radius) ; 
 	
 	
 	material.erase();
 	material.insert(0,"Copper");
     	direction.erase();
   	direction.insert(0,"Transversal");
    	mode.erase();
    	mode.insert(0,"AC");
    	
    	ofstream W3( "wake_Tr_Circ_Cu_AC_5_CPP.dat"); 
    	ofstream FFT3( "wakeFFT_Tr_Circ_Cu_AC_5_CPP.dat");  	
 	CalcWakeFFT(&W3, &FFT3, Lbunch, material, direction, mode, radius) ; 
 	
 	
 	material.erase();
 	material.insert(0,"Copper");
     	direction.erase();
   	direction.insert(0,"Transversal");
    	mode.erase();
    	mode.insert(0,"DC");
    	
    	ofstream W4( "wake_Tr_Circ_Al_DC_5_CPP.dat"); 
    	ofstream FFT4( "wakeFFT_Tr_Circ_Al_DC_5_CPP.dat"); 
 	CalcWakeFFT(&W4, &FFT4, Lbunch, material, direction, mode, radius) ; 
 	
 	
 	material.erase();
 	material.insert(0,"Alu");
     	direction.erase();
   	direction.insert(0,"Longitudinal");
    	mode.erase();
    	mode.insert(0,"AC");
    	ofstream W5( "wake_Lo_Circ_Al_AC_5_CPP.dat"); 
    	ofstream FFT5( "wakeFFT_Lo_Circ_Al_AC_5_CPP.dat"); 
 	CalcWakeFFT(&W5, &FFT5, Lbunch, material, direction, mode, radius) ; 
 	
 	
 	material.erase();
 	material.insert(0,"Alu");
     	direction.erase();
   	direction.insert(0,"Longitudinal");
    	mode.erase();
    	mode.insert(0,"DC");
    	
    	ofstream W6( "wake_Lo_Circ_Al_DC_5_CPP.dat"); 
    	ofstream FFT6( "wakeFFT_Lo_Circ_Al_DC_5_CPP.dat"); 
 	CalcWakeFFT(&W6, &FFT6, Lbunch, material, direction, mode, radius) ; 
 	
 	
 	material.erase();
 	material.insert(0,"Alu");
     	direction.erase();
   	direction.insert(0,"Transversal");
    	mode.erase();
    	mode.insert(0,"AC");
    	
    	ofstream W7( "wake_Tr_Circ_Al_AC_5_CPP.dat"); 
    	ofstream FFT7( "wakeFFT_Tr_Circ_Al_AC_5_CPP.dat");  	
 	CalcWakeFFT(&W7, &FFT7, Lbunch, material, direction, mode, radius) ; 
 	
 	
 	material.erase();
 	material.insert(0,"Alu");
     	direction.erase();
   	direction.insert(0,"Transversal");
    	mode.erase();
    	mode.insert(0,"DC");
    	
    	ofstream W8( "wake_Tr_Circ_Al_DC_5_CPP.dat"); 
    	ofstream FFT8( "wakeFFT_Tr_Circ_Al_DC_5_CPP.dat"); 
 	CalcWakeFFT(&W8, &FFT8, Lbunch, material, direction, mode, radius) ; 
}



/**
* @brief	Calculate out = FFT(in) using fftw3
*
*
* @param 	in: array of real values to be FFT transformed
* @param 	out: arry of complex relults of the FFT transformation
* @param 	N: length of Array in and the Array out
*
*
*/
void fft(double* in , fftw_complex  *out, int N){
	
	fftw_plan p;
	//in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

	p = fftw_plan_dft_r2c_1d(N, in, out, FFTW_ESTIMATE);


	fftw_execute(p); /* repeat as needed */
	

	fftw_destroy_plan(p);
	
  }
  
/**
* @brief	Calculate the Wakeunktion and the FFT of this Wake funktion and saves them in Files
*
* @param 	WakeStream: Save the Wakefield in this File
* @param 	FFTStream: Save the FFT of the Wakefiled in this File
* @param 	Lbunch: Number of points calculatet in the Wakefield
* @param 	material: The Material of the Wakefield. (either "Alu" or "Copper")
* @param 	direction: "Transversal" or Longitudinal"
* @param 	mode: "AC" or "DC"
* @param 	a_Wake: radius of the Tube [mm]
*/
void CalcWakeFFT(ofstream* WakeStream, ofstream* FFTStream, int Lbunch, string material, string direction, string mode, double a_Wake)
{
	// Set integration propartys
  	double  a, b;
  	unsigned int N;
	a=1;
	b=1000000;
	N=1000000;
    
    
    	// Allocate Memory for the FFT
	fftw_complex  *FftWField;
	double *wfield;
	wfield = (double*) fftw_malloc(sizeof(double) * N);
	FftWField = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
	
	// generate stream to write the Data to file
  	*WakeStream << "#Wakefield calculated in C++" << endl;
  	*WakeStream << "#Material: " << material << endl;
  	*WakeStream << "#Direction: " << direction << endl;
  	*WakeStream << "#Mode: " << mode << endl;
  	*WakeStream << "#Form: Circular" << endl;
  	*WakeStream << "#Tube radius: " << a_Wake << endl;
  	*WakeStream << "#Bunch Length: " << Lbunch << endl;
  	
  	// Calculate the Wakefield
  	#pragma omp parallel for 
  	for (int i=0;i<Lbunch;i++){
  		Wake w(material, direction, mode, a_Wake, i*0.000001);
  		//Wake w(i*0.000001);
   	 	wfield[i] = simpson(w,a,b,N);
  		*WakeStream <<i+1 << "    "<< wfield[i] << endl;
  	}
  	(*WakeStream).close(); 
  
  	// stream to write the Data to file
	*FFTStream << "#FFT of the Wakefield calculated in C++"<< endl;
	*FFTStream << "#Wakefield parameter:"<< endl;
  	*FFTStream << "#Material: " << material << endl;
  	*FFTStream << "#Direction: " << direction << endl;
  	*FFTStream << "#Mode: " << mode << endl;
  	*FFTStream << "#Form: Circular" << endl;
  	*FFTStream << "#Tube radius: " << a_Wake << endl;
  	*FFTStream << "#Bunch Length: " << Lbunch << endl;
	
	// calculate the FFT of the Wakefield
	fft(wfield, FftWField, Lbunch);

	// Write the results in a file
  	for (int i=0;i<Lbunch;i++){
  		*FFTStream <<i << "   "<< sqrt(FftWField[i][0]*FftWField[i][0] + FftWField[i][1]*FftWField[i][1])<< endl; 
  	}
  	(*FFTStream).close(); 
  
  	// Free the memory
	fftw_free(wfield); fftw_free(FftWField);
}


/**
* @brief	Calculate the Wakeunktion and write the result in a array
*
* @param	Lbunch: Number of points calculatet in the Wakefield
* @param 	material: The Material of the Wakefield. (either "Alu" or "Copper")
* @param 	direction: "Transversal" or Longitudinal"
* @param 	mode: "AC" or "DC"
* @param  	a_Wake: radius of the Tube [mm]
* @param 	wake: write the calculatet wake funktion in this array
*/
void CalcWake(int Lbunch, string material, string direction, string mode, double a_Wake, double *wake)
{
	// Set integration propartys
  	double  a, b;
  	unsigned int N;
	a=1;
	b=1000000;
	N=1000000;
    
  	// Calculate the Wakefield
  	#pragma omp parallel for 
  	for (int i=0;i<Lbunch;i++){
  		Wake w(material, direction, mode, a_Wake, i*0.000001);
  		//Wake w(i*0.000001);
   	 	wake[i] = simpson(w,a,b,N);
  	}
}

  

