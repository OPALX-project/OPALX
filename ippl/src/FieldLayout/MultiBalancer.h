// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 * 
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

//////////////////////////////////////////////////////////////////////
// Class MultiBalancer
//////////////////////////////////////////////////////////////////////

#ifndef MULTI_BALANCER_H
#define MULTI_BALANCER_H

//////////////////////////////////////////////////////////////////////
//
// Description
// The guts of the load balancing algorithm for ConejoBalancer.
// Inputs weights for vnodes and figures out where the vnodes should go.
//
//////////////////////////////////////////////////////////////////////

// include files
#ifdef IPPL_STDSTL
#include <vector>
#else
#include <vector.h>
#endif // IPPL_STDSTL

/*

  This is a workhorse class for ConejoBalancer.
  It inputs the weights for all the vnodes in the system and
  calculates what processor they will go to.

  This works only at the vnode level.  The idea is that this class is
  used by a higher level class which knows about Fields and actually
  does the redistribution.

  **************************************************
  General Problem Outline
  **************************************************

  There are M materials
  There are V vnodes
  There are P processors

  Each time step proceeds in a series of phases, one for each
  material.  Each of those phases syncronizes often, so each phase
  must be load balanced.

  Each vnode has a weight associated with each material.  It can
  have nonzero weight associated with any number of materials. 
  The weights represent the amount of computation that must be
  performed for that vnode in each phase.

  The total time for a time step is the sum of the times for each
  material.  The time for each material is the maximum over the 
  processors of the total weight on that processor.  The total weight 
  on each processor is the sum of the weights on that processor. 
  Algebraically that looks like:
    T = sum_m max_p sum_v T_{mpv}
  
  Our goal is to chose a distribution of vnodes on processors that
  minimize T.  This problem is bound to be NP complete, so there are
  two ways to find a solution:
  1. Try all combinations. This is expensive.
  2. Find a good heuristic.  This is the approach taken here.

  **************************************************
  General Algorithm Outline
  **************************************************

  There are three phases.  

  1. For each material input a series of blocks of vnodes.  The idea
  is that the user will input the vnodes from one processor after
  another for each material in turn since the high level interface is
  expected to input a Field of weights for each material.

  2. Tell the balancer to figure out what the new distribution should
  be.  It does this by looping over the vnodes, assigning each one to
  a processor.  It figures out the processor by finding the processor
  that minimizes T above when adding each vnode.

  3. Read out the new destination processor for each vnode.

  The skeleton of how MultiBalancer would be used is then something
  like:

  // Construct with the number of procs and vnodes.
  MultiBalancer balancer(procs,vnodes);

  // Tell the balancer about the vnodes.
  // Loop over materials
  for (int m=0; m<materials; ++m_)
    {
      // Loop over input from each processor
      for (int p=0; p<procs; ++p)
        {
	  // Tell the balancer the weights for the vnodes
	  // on the processor where they start out.
	  balancer.appendWeights(...);
        }
    }
  
  // Tell the balancer to figure out the new distribution.
  balancer.distribute();

  // Read out the new distributions.  
  for (MultiBalancer::iterator newProc=balancer.begin();
       newProc != balancer.end();
       ++newProc)
    {
      // Get the new processor for this vnode.
      ... *newProc ...
    }

*/

// class definition for MultiBalancer
class MultiBalancer
{
public:

  // Constructor.
  // Tell it the number of processors.
  MultiBalancer(int procs, int vnodes);

  // Destructor.
  ~MultiBalancer();

  // Tell it to start recording data for a new material.
  void newMaterial();
  
  // Append weights for the current material on the next processor.
  // Inputs are:
  // begin, end: iterators to loop over the weights.
  void appendWeights(double *begin, double* end);

  // Distribute the vnodes to the processors.
  void distribute();

  // Interface for reading the data out again.
  // What we output is a processor for each
  // vnode, so just hold a vector of integers.
  typedef std::vector<int>::iterator iterator;

  // Begin and end iterators for iterating over the vnodes
  // for outputting the destination processor.
  // The vnodes are read out in the order in which they were inserted.
  iterator begin() { return m_vnodeProcs.begin(); }
  iterator end()   { return m_vnodeProcs.end(); }

private:

  // This thing knows the number of processors, vnodes, materials it has.
  unsigned int m_procs;
  unsigned int m_vnodes;
  unsigned int m_materials;

  // Record what phase it is in.
  int m_phase;
  enum { noMaterials, someMaterials, distributed };
  
  // A vector with one weight per vnode.
  typedef std::vector<double> VnodeWeights_t;

  // A vector with one weight per material.
  typedef std::vector<double> MaterialWeights_t;

  // A vector with one weight per processor.
  typedef std::vector<double> ProcWeights_t;

  // A vector of vectors to hold weights for all materials on each vnode.
  typedef std::vector<VnodeWeights_t*> MaterialVnodeWeights_t;

  // A vector of vectors to hold weights for all materials on each proc.
  typedef std::vector<VnodeWeights_t*> MaterialProcWeights_t;

  // A vector of integers to hold the destination processors.
  typedef std::vector<int> VnodeDestProcs_t;

  // Record of the weights when the come in from other processors.
  MaterialVnodeWeights_t m_inputWeights;

  // The list of processors where each vnode will end up.
  VnodeDestProcs_t m_vnodeProcs;

  // The weight assigned to each processor.
  ProcWeights_t m_procWeights;

  // The maximum weights for each material.
  MaterialWeights_t m_materialMaxWeights;

  // The weight of each material on each processor.
  MaterialProcWeights_t m_materialProcWeights;

  //
  // We need to do an indirect sort.
  // That is, we need to sort a list of integers
  // where the integers are offsets in an array of doubles,
  // and the sort criterion is the size of the doubles.
  //
  // So, we need to build a comparator object to be passed to
  // the sort algorithm.
  //
  struct IndirectComparator
  {
    // A vector with one weight per vnode.
    typedef std::vector<double> VnodeWeights_t;

    // Store a pointer to the weights.
    VnodeWeights_t::const_iterator m_weights;

    // Construct with a pointer to the weights.
    IndirectComparator(VnodeWeights_t::const_iterator p)
      : m_weights(p) {}

    // Given two integers, compare the weights at those offsets.
    bool operator()(unsigned int i, unsigned int j) 
      { return m_weights[i] < m_weights[j]; }
  };

  //
  // Accumulate the total weight per vnode into a container.
  //
  void calcTotalWeights(VnodeWeights_t& vnodeTotalWeights);

  //
  // Build an indirection vector that can be used to access the
  // vnodes in order of total weight.
  //
  void calcIndirectSort(std::vector<unsigned int>& vnodeIndirect,
                        const VnodeWeights_t& vnodeTotalWeights);

  //
  // Find what proc a vnode should be put on.
  //
  unsigned int findProcForVnode(unsigned int vnode);

  // Allow the tester to see in here.
  friend void multiBalancerTester();

};

//////////////////////////////////////////////////////////////////////
#endif //MULTI_BALANCER_H

/***************************************************************************
 * $RCSfile: MultiBalancer.h,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:27 $
 * IPPL_VERSION_ID: $Id: MultiBalancer.h,v 1.1.1.1 2003/01/23 07:40:27 adelmann Exp $ 
 ***************************************************************************/
