#include <gtest/gtest.h>

//// Standard includes. /////////////////////////////////////////////////////

#include <iostream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <vector>
#include <cassert>
#include <cmath>


//// Toolkit includes. //////////////////////////////////////////////////////

#include "wfgProblems/Toolkit/ExampleProblems.h"
#include "wfgProblems/Toolkit/TransFunctions.h"


//// Used namespaces. ///////////////////////////////////////////////////////

using namespace WFG::Toolkit;
using namespace WFG::Toolkit::Examples;
using std::vector;
using std::string;

//** Using a uniform random distribution, generate a number in [0,bound]. ***

double next_double( const double bound = 1.0 )
{
  assert( bound > 0.0 );

  return bound * rand() / static_cast< double >( RAND_MAX );
}

//** Create a random Pareto optimal solution for WFG2-WFG7. *****************

vector< double > WFG_2_thru_7_random_soln( const int k, const int l )
{
  vector< double > result;  // the result vector


  //---- Generate a random set of position parameters.

  for( int i = 0; i < k; i++ )
  {
    result.push_back( next_double() );
  }


  //---- Set the distance parameters.

  for( int i = k; i < k+l; i++ )
  {
    result.push_back( 0.35 );
  }


  //---- Scale to the correct domains.

  for( int i = 0; i < k+l; i++ )
  {
    result[i] *= 2.0*(i+1);
  }


  //---- Done.

  return result;
}

//** Convert a double vector into a string. *********************************

string make_string( const vector< double >& v )
{
  std::ostringstream result;

  if( !v.empty() )
  {
    result << v.front();
  }

  for( int i = 1; i < static_cast< int >( v.size() ); i++ )
  {
    result << " " << v[i];
  }

  return result.str();
}

// processes each front from the file 
TEST(WFGPROBLEMSTest, Alive) {
    const int M = 2;          // the number of objectives
    const int k_factor = 1;   // k (# position parameters) = k_factor*( M-1 ) 
    const int l_factor = 1;   // l (# distance parameters) = l_factor*2

    srand( 0 );  // seed the random number generator

    {
      const int k = k_factor*( M-1 );
      const int l = l_factor*2;

      const vector< double >& z = WFG_2_thru_7_random_soln( k, l);
      const vector< double >& f = Problems::WFG2(z, k, M);

      std::cout << make_string( f ) << std::endl;
    }
}