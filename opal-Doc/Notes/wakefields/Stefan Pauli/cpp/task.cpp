#include "wake.h"
#include "twoParticle.h"
#include "profile.h"
#include <fstream>
#include <string>
#include <vector>
#include <istream>
#include <iostream>  // Neeeded for stream I/O
#include <fstream>   // Needed for file I/O
#include <iomanip>   // Needed for I/O manipulators

//#include "Utilities/OpalException.h"

using namespace std;

void readSDDS1(string filename);

/**
* @brief	This will generate a serie of wakefield with diffrent parameters, and calculate the resulting kick for the particles with a given bunch linedistribution.
*
*/
int main() {
	// Calculate a serie of wakefields with diffrent parameters 
	// and save the results in files
	calcWakeSerie();


	// Calculate a serie of kiks of wakefields with diffrent parameters 
	// and save the results in files
 	calcEnergySerie(profile);


	// read in a SDDS1 file
	readSDDS1("test.sdds");
	return 0;
}



/**
* @brief	This will generate a serie of wakefield with diffrent parameters, and calculate the resulting kick for the particles with a given bunch linedistribution.
*
* @param[in]	filename name of the file in which the wake field is stored in a SDDS1 format 
*/
void readSDDS1(string filename)
{
std::ifstream fs;
fs.open(filename.c_str());

if(fs.fail()){
/* throw OpalException("Distribution::Create()",
                          "Open file operation failed, please check if \""
                          + filename +  "\" really exists.");*/
  cout << "readSDDS1 "<< "Open file operation failed, please check if \""
                          << filename <<  "\" really exists." << endl;
  return;
}

string name;
fs >> name;
cout << "readSDDS" << " SSDS1 read = " <<name<< endl;
if(name.compare("SDDS1") != 0){
  cout << "readSDDS1 "<< "No SDDS1 File. A SDDS1 file should start with a SDDS1 String. Check file \""
                          << filename <<  "\" "<< endl;
  return;
}

char temp[256];
for (int i=0; i<6; i++)
{
 fs.getline(temp,256);
 cout << "readSDDS" << "line " <<i<< " :   "<<temp << endl;
}

int Np;
fs >> Np;
cout << "readSDDS" << " header read" << endl;
if ( Np <= 0 )
   {
      /*throw OpalException("Distribution::Create()",
                          " The particle number should be bigger than zero! Please check the first line of file \""
                          + filename +  "\".");*/
      cout << "readSDDS"<<
                          " The particle number should be bigger than zero! Please check the first line of file \""
                          << filename <<  "\"."<< endl;
    }
    
    cout << "readSDDS" << " Np = " <<Np << endl;
   double wake[Np], dist[Np], t[Np];
    
    // read the wakefunction
    for(unsigned int i=0; i<Np; i++) 
    {      
      if ( !fs.eof() )
      {
        fs>> dist[i] >> wake[i] >> t[i];
      }
      if (fs.eof() )
      {
        /*throw OpalException("Distribution::Create()",
                            "End of file reached before the whole wakefield is imported, please check file \""
                            + filename +  "\".");*/
        cout << "readSDDS" <<
                            "End of file reached before the whole wakefield is imported, please check file \""
                            << filename <<  "\"."<< endl;
        return;
      }
    }
    
    cout << wake  <<endl;
    for (int i=0; i<Np; i++)
    {
      cout << wake[i] << endl;
    }
}
