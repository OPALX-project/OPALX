/* Definitions of real and integer variables: 
       if the code should be in double precision, here's 
       the only place where change should be made. 
   Also, some F90 intrinsics are defined here.
 
                      Zarija Lukic, February 2009
                           zarija@lanl.gov
*/

#ifndef TypesDefs_Header_Included
#define TypesDefs_Header_Included

#ifdef INITIALIZER_NAMESPACE
namespace initializer {
#endif

typedef float real;
typedef int integer;

typedef struct{
	double re;
	double im;
} my_fftw_complex;

inline int MOD(int x, int y) { return (x - y*(integer)(x/y));}

#ifdef INITIALIZER_NAMESPACE
}
#endif

#endif
