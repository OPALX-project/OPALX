#include "wake.h"
#include <iostream>

//#define debug

/**
* @brief	Callculate the energy of the Wakefunction with the profile
*
*
* @param[in] 	K a constant
* @param[in] 	charge a constant
* @param[in] 	wake the wake function
* @param[in] 	LengthWake length of the Wake function
* @param[in] 	profile the distibution of the Particles
* @param[out]	OutEnergy this is the Output
*/

void compEnergy(const double K, const double charge, const double* wake, const int LengthWake, const double* profile, double* OutEnergy){


	fftw_plan p;
	int N = 2*LengthWake -1;
	// Allocate Space for the Result
	double *energy = (double*)fftw_malloc(N*sizeof(double));
	
	// Allocate Space for the zerro paded Wakefield and its Fourier Transformed
	double *wfield = (double*)fftw_malloc(N*sizeof(double));
	fftw_complex* FftWField = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
	
	// Allocate Space for the zerro pades profile and its Fourier Transformed
	double *pProfile = (double*)fftw_malloc(N*sizeof(double));
	fftw_complex* FftProfile = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
	
	fftw_complex* final = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)* N);
	
	#ifdef debug
  		cout << "profile" << endl;
  	#endif
	// fill the arrays with data
	for (int i=0;i<LengthWake;i++){
   	 	wfield[i] = wake[i];
  		pProfile[i] = profile[i];
  		#ifdef debug
  		cout << profile[i] << endl;
  		#endif
  	}
  	// make the Zerro pading
	for (int i=LengthWake;i<N;i++){
   	 	wfield[i] = 0;
  		pProfile[i] = 0;
  	}
	
	// calculate the FFTs
	//FFT of the Wakefield
	p = fftw_plan_dft_r2c_1d(N, wfield, FftWField, FFTW_ESTIMATE);
     	fftw_execute(p); 
	fftw_destroy_plan(p);
	
	//FFT of the profile
	p =  fftw_plan_dft_r2c_1d(N, pProfile, FftProfile, FFTW_ESTIMATE);
	fftw_execute(p);
	fftw_destroy_plan(p);
	
	
	#ifdef debug
  		cout << "FFT of the Profile" << endl;
  		
  		for(int i=0; i< N; i++)
  		{
  			cout << i <<": "<< FftProfile[i][0] << " + " << FftProfile[i][1] << " i" << endl; 
  		}
  		cout << "FFT of the Wakefield" << endl;
  		
  		for(int i=0; i< N; i++)
  		{
  			cout << i <<": "<< FftWField[i][0] << " + " << FftWField[i][1] << " i" << endl; 
  		}
  	#endif
	
	
	
	// calculate the .....
	// till here it is consitent with matlab
	for (int i=0; i<N; i++){		
		final[i][0] = (FftWField[i][0]*FftProfile[i][0]-FftWField[i][1]*FftProfile[i][1]);
		final[i][1] = (FftWField[i][0]*FftProfile[i][1]+FftWField[i][1]*FftProfile[i][0]);
	}
	
	#ifdef debug
	  	cout << "final" << endl;
  		
  		for(int i=0; i< N; i++)
  		{
  			cout << i <<": "<< final[i][0] << " + " << final[i][1] << " i" << endl; 
  		}
  	#endif
	
	// inverse transform to get c, the convolution of a and b;
	//this has the side effect of overwriting C 
	p = fftw_plan_dft_c2r_1d(N, final, energy,FFTW_ESTIMATE);
	fftw_execute(p);
	fftw_destroy_plan(p);
	

	// Write the result to the output:
	for (int i=0; i<LengthWake; i++){
		OutEnergy[i] = -0.000001 * (charge * 1000) * K * energy[i] ;
		OutEnergy[i] /= (2*1000*LengthWake);
	}	

	// Free the Memory
	fftw_free(wfield); fftw_free(FftWField); fftw_free(energy);
} 


/**
* @brief	Calculate the energy of a Wake funktion and the profile and saves them in Files
*
* @param[in]	Stream write the result in this file 	
* @param[in]	K a constant
* @param[in]	charge a constant
* @param[in]	wake the wake function
* @param[in]	Lbunch length of the Wake function
* @param[in]	profile the distibution of the Particles
* @param[in]	material can be either "Copper" or "Alu"	
* @param[in]	direction can be either "Transversal" or "Longitudinal"
* @param[in]	Mode can be either "AC" or "DC"
* @param[in]	radius radius of the tube
*
*/
void compEnergyFile(ofstream* Stream, const double K, const double charge, const double* wake, const int Lbunch, const double* profile, string material, string direction, string mode, double radius)
{
	double OutEnergy[Lbunch];
	compEnergy(K, charge, wake, Lbunch, profile, OutEnergy);
	
	// generate stream to write the Data to file
  	*Stream << "#Wakefield calculated in C++" << endl;
  	*Stream << "#Material: " << material << endl;
  	*Stream << "#Direction: " << direction << endl;
  	*Stream << "#Mode: " << mode << endl;
  	*Stream << "#Form: Circular" << endl;
  	*Stream << "#Tube radius: " << radius << endl;
  	*Stream << "#Bunch Length: " << Lbunch << endl;
  	
  	// Calculate the Wakefield
  	for (int i=0;i<Lbunch;i++){
  		*Stream <<i+1 << "    "<< OutEnergy[i] << endl;
  	}
  	(*Stream).close(); 
}

/**
* @brief	Calculate meany diffrent energys of Wakfields
*
*/
void calcEnergySerie(const double* profile){

	cout.precision(30);
	
	// Set integration propartys
  	double  a, b;
  	unsigned int N;
	a=1;
	b=1000000;
	N=1000000;
    	int Lbunch = 256;
    	double wake[Lbunch];
    	double K = 0.20536314319923724;
	double charge = 0.8;
	double radius = 5;
    	
    	
    	string material ("Copper");
    	string direction("Longitudinal");
    	string mode ("AC");
    	ofstream E1( "energy_Lo_Circ_Cu_AC_5_CPP.dat"); 
    	CalcWake(Lbunch, material, direction, mode,radius, wake);
    	compEnergyFile(&E1, K, charge, wake, Lbunch, profile, material, direction, mode, radius);


 	
 	material.erase();
 	material.insert(0,"Copper");
     	direction.erase();
   	direction.insert(0,"Longitudinal");
    	mode.erase();
    	mode.insert(0,"DC");
 
     	ofstream E2( "energy_Lo_Circ_Cu_DC_5_CPP.dat");    	
    	CalcWake(Lbunch, material, direction, mode,radius, wake);
    	compEnergyFile(&E2, K, charge, wake, Lbunch, profile, material, direction, mode, radius);
 	
 	
 	material.erase();
 	material.insert(0,"Copper");
     	direction.erase();
   	direction.insert(0,"Transversal");
    	mode.erase();
    	mode.insert(0,"AC");
    	
     	ofstream E3( "energy_Tr_Circ_Cu_AC_5_CPP.dat");    	
    	CalcWake(Lbunch, material, direction, mode,radius, wake);
    	compEnergyFile(&E3, K, charge, wake, Lbunch, profile, material, direction, mode, radius);
    
 	
 	material.erase();
 	material.insert(0,"Copper");
     	direction.erase();
   	direction.insert(0,"Transversal");
    	mode.erase();
    	mode.insert(0,"DC");
    	
     	ofstream E4( "energy_Tr_Circ_Cu_DC_5_CPP.dat");    	
    	CalcWake(Lbunch, material, direction, mode,radius, wake);
    	compEnergyFile(&E4, K, charge, wake, Lbunch, profile, material, direction, mode, radius);
    	
 	
 	material.erase();
 	material.insert(0,"Alu");
     	direction.erase();
   	direction.insert(0,"Longitudinal");
    	mode.erase();
    	mode.insert(0,"AC");
    	
     	ofstream E5( "energy_Lo_Circ_Al_AC_5_CPP.dat");    	
    	CalcWake(Lbunch, material, direction, mode,radius, wake);
    	compEnergyFile(&E5, K, charge, wake, Lbunch, profile, material, direction, mode, radius);
    	
 	
 	material.erase();
 	material.insert(0,"Alu");
     	direction.erase();
   	direction.insert(0,"Longitudinal");
    	mode.erase();
    	mode.insert(0,"DC");
    	
     	ofstream E6( "energy_Lo_Circ_Al_DC_5_CPP.dat");    	
    	CalcWake(Lbunch, material, direction, mode,radius, wake);
    	compEnergyFile(&E6, K, charge, wake, Lbunch, profile, material, direction, mode, radius);
    	
 	
 	material.erase();
 	material.insert(0,"Alu");
     	direction.erase();
   	direction.insert(0,"Transversal");
    	mode.erase();
    	mode.insert(0,"AC");
    	
     	ofstream E7( "energy_Tr_Circ_Al_AC_5_CPP.dat");    	
    	CalcWake(Lbunch, material, direction, mode,radius, wake);
    	compEnergyFile(&E7, K, charge, wake, Lbunch, profile, material, direction, mode, radius);
    	
 	
 	material.erase();
 	material.insert(0,"Alu");
     	direction.erase();
   	direction.insert(0,"Transversal");
    	mode.erase();
    	mode.insert(0,"DC");
    	
     	ofstream E8( "energy_Tr_Circ_Al_DC_5_CPP.dat");    	
    	CalcWake(Lbunch, material, direction, mode,radius, wake);
    	compEnergyFile(&E8, K, charge, wake, Lbunch, profile, material, direction, mode, radius);
    
}

