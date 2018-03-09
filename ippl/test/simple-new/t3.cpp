/***************************************************************************
 *
 * 
 * 
 *
 ***************************************************************************/

#include <iostream>
#include <vector>

#include "Ippl.h"

int main(int argc, char *argv[])
{
  Ippl ippl(argc,argv);
  Inform msg("t3 ", INFORM_ALL_NODES);

  Inform msgmaster("master ");

  int Parent = 0;
  int tag = 101;

  int vnodes, sizeX, sizeY, centerX, centerY, iterations;

  std::vector<int> data(10);

  int arrayData[10];

  size_t abig;


  if( Ippl::Comm->myNode() == Parent ) {
    vnodes = 100; 
    sizeX = 64;
    sizeY = 128;
    centerX = 0;
    centerY = 1;
    abig = 800000000000000;
    iterations = 1000;
    for (int i=0;i<10;i++) {
      data[i] = (3*i);
      arrayData[i]=3*i;
    }

    // now broadcast data to other nodes
    Message *mess = new Message();
    putMessage( *mess, abig );
    putMessage( *mess, vnodes );
    putMessage( *mess, sizeX );
    putMessage( *mess, sizeY );
    putMessage( *mess, centerX );
    putMessage( *mess, centerY );
    putMessage( *mess, iterations );
    putMessage( *mess, arrayData, arrayData + 10 );
    putMessage( *mess, &data[0], &data[10] );
    Ippl::Comm->broadcast_all(mess, tag);

  }

  // now each node receives the data
  Message *mess = Ippl::Comm->receive_block(Parent, tag);
  PAssert(mess);
  getMessage( *mess, abig );
  getMessage( *mess, vnodes );
  getMessage( *mess, sizeX );
  getMessage( *mess, sizeY );
  getMessage( *mess, centerX );
  getMessage( *mess, centerY );
  getMessage( *mess, iterations );
  getMessage( *mess, arrayData );
  getMessage( *mess, &data[0] );

  msg << "std::vector content: " ;
  for (int i=0;i<10;i++)
      msg << data[i] << " "; 
  msg << endl;
  
  msg << "array content:       " ;
  for (int i=0;i<10;i++)
      msg << arrayData[i] << " ";
  msg << endl;

  delete mess;
   
  msg << "abig= " << abig << endl;

  msg << "received a message on node " << Ippl::Comm->myNode();
  msg << " from node " << Parent << " with tag " << tag << endl;


  reduce(abig,abig,OpAddAssign());
  msg << "sum(abig_= " << abig << endl;


  return 0;
}
/***************************************************************************
 * $RCSfile: t3.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/03/28 12:47:26 $
 * IPPL_VERSION_ID: $Id: t3.cpp,v 1.1.1.1 2003/03/28 12:47:26 adelmann Exp $ 
 ***************************************************************************/
