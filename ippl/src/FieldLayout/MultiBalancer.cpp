// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 * 
 * This program was prepared by PSI. 
 * All rights in the program are reserved by PSI.
 * Neither PSI nor the author(s)
 * makes any warranty, express or implied, or assumes any liability or
 * responsibility for the use of this software
 *
 * Visit www.amas.web.psi for more details
 *
 ***************************************************************************/

// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 * 
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

// include files
#include "FieldLayout/MultiBalancer.h"
#include "Utility/PAssert.h"

// Standard STL names this header "algorithm"
#include <algorithm>
using namespace std;

//////////////////////////////////////////////////////////////////////

//
// Construct a MultiBalancer given:
// The number of processors
// The number of vnodes.
//
// We don't actually allocate space in the arrays here.
//

MultiBalancer::MultiBalancer(int procs, int vnodes)
: 
  // Record the number of procs and vnodes.
  m_procs(procs) , 
  m_vnodes(vnodes) ,
  // The number of materials starts out at zero.
  m_materials(0) ,
  // It starts out with no materials.
  m_phase( noMaterials ),
  // Set the size of the container of vnode destination procs.
  m_vnodeProcs(vnodes),
  // Set the size of the container of processor weights.
  m_procWeights(procs)
{
}

//////////////////////////////////////////////////////////////////////

//
// The dtor for MultiBalancer
//

MultiBalancer::~MultiBalancer()
{
  // Delete the nested vectors
  for (unsigned int i=0; i<m_materials; ++i)
    {
      delete m_inputWeights[i];
      delete m_materialProcWeights[i];
    }
}
  
//////////////////////////////////////////////////////////////////////

//
// Start recording information for a new material.
//

void MultiBalancer::newMaterial()
{
  // Make sure we're not in the wrong mode.
  PAssert( m_phase != distributed );
  
  // Do a sanity check on the previous material,
  // provided that this is not the first one.
  if ( m_phase == someMaterials )
    // Make sure we got the right number of vnodes.
    // Note that braces are added since PAssert can be empty
    {PAssert( m_inputWeights.back()->size() == m_vnodes );}

  // Add a new vector of weights to the materials container.
  m_inputWeights.push_back( new VnodeWeights_t );

  // Make it can hold the right number of vnodes.
  m_inputWeights.back()->reserve( m_vnodes );

  // Increment the number of materials.
  ++m_materials;

  // Set the new mode.
  m_phase = someMaterials;
}

//////////////////////////////////////////////////////////////////////

//
// Record weights from the vnodes on the next processor.
//

void MultiBalancer::appendWeights(double *begin, double *end)
{
  // Make sure we have a nonzero number of materials.
  PAssert(m_materials>0);

  // Get a reference to the container where we'll be putting these.
  VnodeWeights_t& weights = *m_inputWeights.back();

  // Loop over the inputs, appending to the container.
  for (; begin!=end; ++begin)
    weights.push_back(*begin);

  // Make sure we didn't make it longer than it should be.
  PAssert( weights.size() <= m_vnodes );
}

//////////////////////////////////////////////////////////////////////

//
// Calculate the redistribution.
//

void MultiBalancer::distribute()
{
  // Make sure we're in the right mode.
  PAssert( m_phase != noMaterials );

  // If we've already distributed, we don't need to do anything.
  if ( m_phase == distributed )
    return;

  unsigned int mat_index;

  // Build a vector of the total weight per vnode.
  VnodeWeights_t vnodeTotalWeights(m_vnodes,0.0);
  calcTotalWeights(vnodeTotalWeights);

  // Build an indirection vector so we can treat the vnodes
  // in descending order of total weight.
  vector<unsigned int> vnodeIndirect(m_vnodes);
  calcIndirectSort(vnodeIndirect,vnodeTotalWeights);

  // Keep track of the maximum weight for each material as we go along.
  m_materialMaxWeights.resize( m_materials );
  for (mat_index=0; mat_index<m_materials; ++mat_index)
    m_materialMaxWeights[mat_index] = 0.0;

  // Keep track of the total weight of each material on each processor.
  m_materialProcWeights.reserve(m_materials);
  for (mat_index=0; mat_index<m_materials; ++mat_index)
    m_materialProcWeights.push_back(new VnodeWeights_t(m_procs));

  // Loop over the vnodes in decreasing order of weight,
  // finding the processor for each one.

  for (unsigned int i=m_vnodes; i-->0;) 
    {
      // Which vnode is this?
      int vnode = vnodeIndirect[i];

      // Find the processor for which this weight
      // minimizes the maximum weight across processors.
      int proc = findProcForVnode(vnode);

      // Record that proc.
      m_vnodeProcs[vnode] = proc;

      // Record the total weight for that processor.
      m_procWeights[proc] += vnodeTotalWeights[vnode];

      // Record the weights for each material.
      for (mat_index=0; mat_index<m_materials; ++mat_index)
        {
          // Accumulate the weight into the weight for this proc.
          m_materialProcWeights[mat_index][0][proc] += m_inputWeights[mat_index][0][vnode];

          // See if this proc is now the maximum weight for this material.
          if ( m_materialMaxWeights[mat_index] < m_materialProcWeights[mat_index][0][proc] )
            m_materialMaxWeights[mat_index] = m_materialProcWeights[mat_index][0][proc];
        }
    }

  // Set the new mode.
  m_phase = distributed;
}

//////////////////////////////////////////////////////////////////////

//
// Find the processor we should add this vnode to.
// This is the processor that will minimize the increase 
// in the total weight. 
//

unsigned int MultiBalancer::findProcForVnode( unsigned int v )
{
  // The minimum of the weights on a processer when this vnode is included.
  double minWeight = 1e30;

  // Which proc has the minimum weight when this vnode is included.
  unsigned int minProc = m_procs;

  // Try putting this vnode on each processor in turn.
  for (unsigned int proc_index = 0; proc_index<m_procs; ++proc_index) 
    {
      // Find the estimated total weight if this vnode is put on this proc.
      double totalWeight = 0.0;
      for (unsigned int mat_index=0; mat_index<m_materials; ++mat_index)
        {
          // Add the weight from this material on this vnode
          // to the weight from this material already on this proc.
          double w = m_materialProcWeights[mat_index][0][proc_index]+m_inputWeights[mat_index][0][v];

          // If the weight is less than the maximum weight on other
          // procs, then the other proc determines the weight.
          if ( w < m_materialMaxWeights[mat_index] )
            w = m_materialMaxWeights[mat_index];

          // Accumulate the weight into the total.
          totalWeight += w;
        }

      // If the new weight for this proc is less than the previous
      // least weight, then this is the new best candidate.
      if ( totalWeight < minWeight )
        {
          minWeight = totalWeight;
          minProc = proc_index;
        }
    }

  // Return the best candidate processor found.
  // We had better have found something.
  PAssert(minProc<m_procs);
  return minProc;
}

//////////////////////////////////////////////////////////////////////

//
// Calculate the total weight per vnode.
//
void 
MultiBalancer::calcTotalWeights(VnodeWeights_t& vnodeTotalWeights)
{
  MaterialVnodeWeights_t::iterator matp=m_inputWeights.begin();
  for (; matp != m_inputWeights.end(); ++matp)
    {
      // Do a sanity check.
      // Make sure that all of the input materials have the same
      // number of vnodes.
      PAssert( (*matp)->size() == m_vnodes );

      // Accumulate the weight of each vnode into the 
      // array of total weights.
      for (unsigned int i = 0; i<m_vnodes; ++i)
        vnodeTotalWeights[i] += (**matp)[i];
    }
}

//////////////////////////////////////////////////////////////////////

//
// Calculate a list of integers we can use to 
// access the vnodes in sorted order.
//
void
MultiBalancer::calcIndirectSort(vector<unsigned int>& vnodeIndirect,
                                const VnodeWeights_t& vnodeTotalWeights)
{
  // Initialize the indirection list with increasing integers.
  for (unsigned int i=0; i<m_vnodes; ++i)
    vnodeIndirect[i] = i;

  // Sort the indirection list using a comparator
  // which looks up values in the weight array.
  sort(
    vnodeIndirect.begin(), 
    vnodeIndirect.end(), 
    IndirectComparator(vnodeTotalWeights.begin()));
}


/***************************************************************************
 * $RCSfile: MultiBalancer.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:27 $
 * IPPL_VERSION_ID: $Id: MultiBalancer.cpp,v 1.1.1.1 2003/01/23 07:40:27 adelmann Exp $ 
 ***************************************************************************/
