/*
 *  Output.h
 *  Initializer
 *
 *  Created by Zarija on 1/27/10.
 *
 */

#ifndef Output_Header_Included
#define Output_Header_Included

#include <string>
#include "InputParser.h"
#include <memory>

#ifdef USENAMESPACE
namespace initializer {
#endif
using namespace std;
   
class Output {

public:   
   void grid2phys(std::unique_ptr<real[]> &pos_x, std::unique_ptr<real[]> &pos_y, std::unique_ptr<real[]> &pos_z, 	 
                  std::unique_ptr<real[]> &vel_x, std::unique_ptr<real[]> &vel_y, std::unique_ptr<real[]> &vel_z, 	 
                  int Npart, int np, float rL);

   void write_hcosmo_ascii(std::unique_ptr<real[]> &pos_x, std::unique_ptr<real[]> &pos_y, std::unique_ptr<real[]> &pos_z,
                           std::unique_ptr<real[]> &vel_x, std::unique_ptr<real[]> &vel_y, std::unique_ptr<real[]> &vel_z,
                           int Npart, string outBase);

   void write_hcosmo_serial(std::unique_ptr<real[]> &pos_x, std::unique_ptr<real[]> &pos_y, std::unique_ptr<real[]> &pos_z, 	 
                            std::unique_ptr<real[]> &vel_x, std::unique_ptr<real[]> &vel_y, std::unique_ptr<real[]> &vel_z, 	 
                            std::unique_ptr<integer[]> &id, int Npart, string outBase);

   void write_hcosmo_parallel(std::unique_ptr<real[]> &pos_x, std::unique_ptr<real[]> &pos_y, std::unique_ptr<real[]> &pos_z, 	 
                              std::unique_ptr<real[]> &vel_x, std::unique_ptr<real[]> &vel_y, std::unique_ptr<real[]> &vel_z, 	 
                              std::unique_ptr<integer[]> &id, int Npart, string outBase);
#ifdef H5PART

   void write_hcosmo_h5(std::unique_ptr<real[]> &pos_x, std::unique_ptr<real[]> &pos_y, std::unique_ptr<real[]> &pos_z,
                        std::unique_ptr<real[]> &vel_x, std::unique_ptr<real[]> &vel_y, std::unique_ptr<real[]> &vel_z,
                        std::unique_ptr<integer[]> &id, int Npart, string outBase, 
                        h5part_int64_t ng, h5part_int64_t ng2d, h5part_int64_t np, 
                        h5part_int64_t rL, string indatName);
#endif

   void write_gadget(InputParser& par, std::unique_ptr<real[]> &pos_x, std::unique_ptr<real[]> &pos_y, std::unique_ptr<real[]> &pos_z,
                     std::unique_ptr<real[]> &vel_x, std::unique_ptr<real[]> &vel_y, std::unique_ptr<real[]> &vel_z, std::unique_ptr<integer[]> &id, 
                     int Npart, string outBase);

private:
   
   struct gadget_header
   {
      int      npart[6];
      double   mass[6];
      double   time;
      double   redshift;
      int      flag_sfr;
      int      flag_feedback;
      unsigned int  npartTotal[6];
      int      flag_cooling;
      int      num_files;
      double   BoxSize;
      double   Omega0;
      double   OmegaLambda;
      double   HubbleParam;
      int      flag_stellarage;
      int      flag_metals;
      unsigned int  npartTotalHighWord[6];
      int      flag_entropy_instead_u;
      char     fill[60];  /* fills to 256 Bytes */
   };
};

#ifdef USENAMESPACE
}   
#endif
   
#endif
