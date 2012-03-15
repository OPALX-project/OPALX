// ------------------------------------------------------------------------
// $RCSfile: PartBunch.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class PartBunch
//   Interface to a particle bunch.
//   Can be used to avoid use of a template in user code.
//
// ------------------------------------------------------------------------
// Class category: Algorithms
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 18:57:53 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include "Algorithms/PartBunch.h"
#include "FixedAlgebra/FMatrix.h"
#include "FixedAlgebra/FVector.h"
#include "FixedAlgebra/FVps.h"
#include "FixedAlgebra/LinearMap.h"
#include <iostream>
#include <cfloat>
#include <fstream>
#include "Physics/Physics.h"
#include "AbstractObjects/OpalData.h"

#include "ListElem.h"

using Physics::c;
using Physics::pi;

// Class PartBunch
// ------------------------------------------------------------------------

PartBunch::PartBunch(const PartData *ref):
  reference(ref)
{
  addAttribute(X);
  addAttribute(P);
  addAttribute(Q);
  addAttribute(Ef);
  addAttribute(Eftmp);

  addAttribute(Bf);
  addAttribute(Bin);
  addAttribute(dt);
  addAttribute(LastSection);
  selfFieldTimer_m = IpplTimings::getTimer("computeSelfField");
  boundpTimer_m = IpplTimings::getTimer("Compute bounding box");
  statParamTimer_m = IpplTimings::getTimer("Compute Statistics");
  compPotenTimer_m  = IpplTimings::getTimer("Compute Potential");

  distrCreate_m = IpplTimings::getTimer("creatDistr");
  distrReload_m = IpplTimings::getTimer("loadDistr");

  dh_m = 0;
  fs_m = 0;
  dt = 0.0;
  pbin_m = 0;
  trackStep_m = 0;
  stepsPerTurn_m = 0;
  numBunch_m = 1;

  lineDensity_m = 0;
  nBinsLineDensity_m = 0;

  partPerNode_m = new double[Ippl::getNodes()];
  globalPartPerNode_m = new double[Ippl::getNodes()];
}


PartBunch::PartBunch(const PartBunch &rhs):
  reference(rhs.reference)
{}


PartBunch::PartBunch(const std::vector<Particle> &rhs, const PartData *ref):
  reference(ref)
{}

bool PartBunch::hasFieldSolver() {
  if  (fs_m)
    return fs_m->hasValidSolver();
  else
    return false;
}

PartBunch::~PartBunch()
{
  if(bingamma_m)
    delete bingamma_m;

  if(binemitted_m)
    delete binemitted_m;

  if(lineDensity_m)
    delete lineDensity_m;


  if(partPerNode_m)
    delete partPerNode_m;

  if(partPerNode_m)
    delete globalPartPerNode_m;
}

bool PartBunch::hasZeroNLP() {
  /**
     Check if a node has no particles
  */
  Inform m("hasZeroNLP() ", INFORM_ALL_NODES);
  int minnlp = 0;
  int nlp=getLocalNum();
  minnlp=100000;
  reduce(nlp,minnlp,OpMinAssign());
  return (minnlp==0);
}

void PartBunch::calcLineDensity() {
  /**
   * \method calcLineDensity()
   * \brief calculates the 1d line density (not normalized) and append it to a file.
   * \see ParallelTTracker
   * \warning none yet
   *
   * DETAILED TODO
   *
   */
  e_dim_tag decomp[3];
  list<ListElem> listz;

  NDIndex<Dim> domain = getFieldLayout().getDomain();

  for (int d=0; d < 3; ++d) {
    decomp[d] = getFieldLayout().getRequestedDistribution(d);
  }
  FieldLayout_t *FL  = new FieldLayout_t(getMesh(), decomp);

  if (lineDensity_m == 0) {
    if (nBinsLineDensity_m == 0)
      nBinsLineDensity_m = nr_m[2];
    lineDensity_m = new double[nBinsLineDensity_m];
  }

  for (int i=0; i < nBinsLineDensity_m; ++i)
    lineDensity_m[i] = 0.0;

  tmpFieldZ_m.initialize(getMesh(), *FL);
  tmpFieldZ_m[domain] = this->rho_m[domain];

  this->Q.scatter(this->tmpFieldZ_m, this->R, IntrplCIC_t());

  NDIndex<Dim> idx = FL->getLocalNDIndex();
  NDIndex<Dim> idxdom = FL->getDomain();
  NDIndex<Dim> elem;

  int tag = Ippl::Comm->next_tag(IPPL_APP_TAG0, IPPL_APP_CYCLE);
  if (Ippl::myNode() == 0) {
    for (int i=idx[2].min(); i<=idx[2].max(); ++i) {
      double localquantsum=0.0;
      elem[2]=Index(i,i);
      for (int j=idx[1].min(); j<=idx[1].max(); ++j) {
        elem[1]=Index(j,j);
        for (int k=idx[0].min(); k<=idx[0].max(); ++k) {
          elem[0]=Index(k,k);
          localquantsum += tmpFieldZ_m.localElement(elem);
        }
      }
      listz.push_back(ListElem(get_sPos(),getT(),i+1,i+1,localquantsum));

      // wait for msg from all processors (EXEPT NODE 0)
      int notReceived =  Ippl::getNodes() - 1;
      unsigned int dataBlocks=0;
      unsigned  coor;
      double projVal;
      while (notReceived > 0) {
        int node = COMM_ANY_NODE;
        Message* rmsg =  Ippl::Comm->receive_block(node, tag);
        if (rmsg == 0) {
          ERRORMSG("Could not receive from client nodes in main." << endl);
        }
        notReceived--;
        rmsg->get(&dataBlocks);
        for (unsigned i=0; i < dataBlocks; i++) {
          rmsg->get(&coor);
          rmsg->get(&projVal);
          listz.push_back(ListElem(get_sPos(),getT(),coor,coor,projVal));
        }
        delete rmsg;
      }
    }
    listz.sort();
    /// copy line density to array of double
    fillArray(lineDensity_m,listz);
  }
  else {
    Message* smsg = new Message();
    smsg->put(idx[2].max()-idx[2].min()+1);
    for (int i=idx[2].min(); i<=idx[2].max(); ++i) {
      double localquantsum = 0.0;
      elem[2] = Index(i,i);
      for (int j=idx[1].min(); j<=idx[1].max(); ++j) {
        elem[1] = Index(j,j);
        for (int k=idx[0].min(); k<=idx[0].max(); ++k) {
          elem[0] = Index(k,k);
          localquantsum +=  tmpFieldZ_m.localElement(elem);
        }
      }
      smsg->put(i+1);
      smsg->put(localquantsum);
    }
    bool res = Ippl::Comm->send(smsg, 0, tag);
    if (! res)
      ERRORMSG("Ippl::Comm->send(smsg, 0, tag) failed " << endl;);
  }
  /// finaly make sure all nodes have the line density.
  reduce(lineDensity_m, lineDensity_m + (int)nBinsLineDensity_m, lineDensity_m, OpAddAssign());
}

void PartBunch::fillArray(double *lineDensity, const list<ListElem> &l) {
  unsigned int mmax=0;
  unsigned int nmax=0;
  unsigned int count=0;

  for(list<ListElem>::const_iterator it=l.begin(); it!=l.end() ;++it)  {
    if (it->m > mmax) mmax=it->m;
    if (it->n > nmax) nmax=it->n;
  }
  for(list<ListElem>::const_iterator it=l.begin();it!=l.end();++it)
    if ((it->m < mmax)&&(it->n < nmax)) {
      lineDensity[count] = it->den;
      count++;
    }
}

void PartBunch::getLineDensity(vector<double> &lineDensity)
{
  if (lineDensity_m)
    {
      if (lineDensity.size() != nBinsLineDensity_m)
        lineDensity.resize(nBinsLineDensity_m,0.0);
      for (int i  = 0; i < nBinsLineDensity_m; ++i)
        lineDensity[i] = lineDensity_m[i];
    }
}

void PartBunch::updateBinStructure() {

  pbin_m->updateBinStructure();

  if (Ippl::myNode() !=0) {
    for (int i=0; i<=pbin_m->getLastemittedBin(); i++)
      binemitted_m[i]= 0;
  }
  // only node 0 is "emitting"
  for (int i=0; i<=pbin_m->getLastemittedBin(); i++)
    reduce(binemitted_m[i],binemitted_m[i],OpAddAssign());
}

void PartBunch::calcGammas() {

  const int emittedBins = pbin_m->getLastemittedBin();

  for (int i=0; i<emittedBins+1; i++)
    bingamma_m[i] = 0.0;

  for (int n=0; n< getLocalNum(); n++)
    bingamma_m[this->Bin[n]] += sqrt(1.0 + dot(this->P[n],this->P[n]));

  for (int i=0; i<emittedBins+1; i++) {
    reduce(bingamma_m[i],bingamma_m[i],OpAddAssign());
    bingamma_m[i] /= binemitted_m[i];
    INFOMSG("Bin " << i << " gamma = " << bingamma_m[i] << endl);
  }

  if (emittedBins >= 2){
    for (int i=1; i<emittedBins+1; i++) {
      INFOMSG("dE= " << Physics::m_e*1.0E6*(bingamma_m[i-1]-bingamma_m[i]) << " [keV] of Bin " << i-1 << " and " << i << endl);
    }
  }
}


void PartBunch::calcGammas_cycl() {

  const int emittedBins = pbin_m->getLastemittedBin();

  for (int i=0; i<emittedBins; i++)
    bingamma_m[i] = 0.0;

  for (int n=0; n< getLocalNum(); n++)
    bingamma_m[this->Bin[n]] += sqrt(1.0 + dot(this->P[n],this->P[n]));

  for (int i=0; i<emittedBins; i++) {
    reduce(bingamma_m[i],bingamma_m[i],OpAddAssign());
    bingamma_m[i] /= pbin_m->getBinCont(i);

    INFOMSG("Bin " << i <<" : particle number="<<pbin_m->getBinCont(i) <<" gamma = " << bingamma_m[i] << endl);
  }

}


double PartBunch::getMaxdEBins() {

  const int emittedBins = pbin_m->getLastemittedBin();

  double maxdE = DBL_MIN;
  double mindE = DBL_MAX;
  double maxdEGlobal = DBL_MIN;
  double mindEGlobal = DBL_MAX;
  if (emittedBins >= 2){
    for (int i=1; i<emittedBins+1; i++) {
      double de = fabs(Physics::m_e * 1.0E6 * (bingamma_m[i-1] - bingamma_m[i]));
      if (de > maxdE)
        maxdE=de;
    }

    reduce(maxdE,maxdEGlobal,OpMaxAssign());

    return maxdEGlobal;
  }
  else
    return DBL_MAX;
}


void PartBunch::computeSelfFields(int binNumber)
{
  IpplTimings::startTimer(selfFieldTimer_m);

  /// Set initial charge density to zero. Create image charge
  /// potential field.
  rho_m = 0.0;
  Field_t imagePotential = rho_m;

  /// Set initial E field to zero.
  eg_m = Vector_t(0.0);

  if (fs_m->hasValidSolver()) {

    /// Scatter charge onto space charge grid.
    this->Q *= this->dt;
    this->Q.scatter(this->rho_m, this->R, IntrplCIC_t());
    this->Q /= this->dt;
    this->rho_m /= getdT();

    /// Calculate mesh-scale factor and get gamma for this specific slice (bin).
    double scaleFactor = Physics::c * getdT();
    double gammaz = getBinGamma(binNumber);

    /// Scale charge density to get charge density in real units. Account for
    /// Lorentz transformation in longitudinal direction.
    double tmp2 = 1/hr_m[0] * 1/hr_m[1] * 1/hr_m[2] / (scaleFactor*scaleFactor*scaleFactor) / gammaz;
    rho_m *= tmp2;

    /// Scale mesh spacing to real units (meters). Lorentz transform the
    /// longitudinal direction.
    Vector_t hr_scaled = hr_m * Vector_t(scaleFactor);
    hr_scaled[2] *= gammaz;

    /// Find potential from charge in this bin (no image yet) using Poisson's
    /// equation (without coefficient: -1/(eps)).
    IpplTimings::startTimer(compPotenTimer_m);
    imagePotential = rho_m;
    fs_m->solver_m->computePotential(rho_m, hr_scaled);
    IpplTimings::stopTimer(compPotenTimer_m);

    /// Scale mesh back (to same units as particle locations.)
    rho_m *= hr_scaled[0]*hr_scaled[1]*hr_scaled[2];

    /// The scalar potential is given back in rho_m
    /// and must be converted to the right units.
    rho_m *= getCouplingConstant();

    /// IPPL Grad numerical computes gradient to find the
    /// electric field (in bin rest frame).
    eg_m = -Grad(rho_m,eg_m);

    /// Scale field. Combine Lorentz transform with conversion to proper field
    /// units.
    eg_m *= Vector_t(gammaz/(scaleFactor), gammaz/(scaleFactor), 1.0/(scaleFactor*gammaz));

    /// Interpolate electric field at particle positions.  We reuse the
    /// cached information about where the particles are relative to the
    /// field, since the particles have not moved since this the most recent
    /// scatter operation.
    Eftmp.gather(eg_m, this->R, IntrplCIC_t());

    /** Magnetic field in x and y direction induced by the electric field.
     *
     *  \f[ B_x = \gamma(B_x^{'} - \frac{beta}{c}E_y^{'}) = -\gamma \frac{beta}{c}E_y^{'} = -\frac{beta}{c}E_y \f]
     *  \f[ B_y = \gamma(B_y^{'} - \frac{beta}{c}E_x^{'}) = +\gamma \frac{beta}{c}E_x^{'} = +\frac{beta}{c}E_x \f]
     *  \f[ B_z = B_z^{'} = 0 \f]
     *
     */
    double betaC = sqrt(gammaz*gammaz - 1.0)/gammaz/Physics::c;

    Bf(0) = Bf(0) - betaC * Eftmp(1);
    Bf(1) = Bf(1) + betaC * Eftmp(0);

    Ef += Eftmp;

    /// Now compute field due to image charge. This is done separately as the image charge
    /// is moving to -infinity (instead of +infinity) so the Lorentz transform is different.

    /// Find z shift for shifted Green's function.
    Vector_t rmax, rmin;
    get_bounds(rmin,rmax);
    double zshift = - (rmax(2) + rmin(2))*gammaz*scaleFactor;

    /// Find potential from image charge in this bin using Poisson's
    /// equation (without coefficient: -1/(eps)).
    IpplTimings::startTimer(compPotenTimer_m);
    fs_m->solver_m->computePotential(imagePotential, hr_scaled, zshift);
    IpplTimings::stopTimer(compPotenTimer_m);

    /// Scale mesh back (to same units as particle locations.)
    imagePotential *= hr_scaled[0]*hr_scaled[1]*hr_scaled[2];

    /// The scalar potential is given back in rho_m
    /// and must be converted to the right units.
    imagePotential *= getCouplingConstant();

    /// IPPL Grad numerical computes gradient to find the
    /// electric field (in rest frame of this bin's image
    /// charge).
    eg_m = -Grad(imagePotential,eg_m);

    /// Scale field. Combine Lorentz transform with conversion to proper field
    /// units.
    eg_m *= Vector_t(gammaz/(scaleFactor), gammaz/(scaleFactor), 1.0/(scaleFactor*gammaz));

    /// Interpolate electric field at particle positions.  We reuse the
    /// cached information about where the particles are relative to the
    /// field, since the particles have not moved since this the most recent
    /// scatter operation.
    Eftmp.gather(eg_m, this->R, IntrplCIC_t());

    /** Magnetic field in x and y direction induced by the image charge electric field. Note that beta will have
     *  the opposite sign from the bunch charge field, as the image charge is moving in the opposite direction.
     *
     *  \f[ B_x = \gamma(B_x^{'} - \frac{beta}{c}E_y^{'}) = -\gamma \frac{beta}{c}E_y^{'} = -\frac{beta}{c}E_y \f]
     *  \f[ B_y = \gamma(B_y^{'} - \frac{beta}{c}E_x^{'}) = +\gamma \frac{beta}{c}E_x^{'} = +\frac{beta}{c}E_x \f]
     *  \f[ B_z = B_z^{'} = 0 \f]
     *
     */
    Bf(0) = Bf(0) + betaC * Eftmp(1);
    Bf(1) = Bf(1) - betaC * Eftmp(0);

    Ef += Eftmp;
  }
  IpplTimings::stopTimer(selfFieldTimer_m);
}


void PartBunch::computeSelfFields()
{
  IpplTimings::startTimer(selfFieldTimer_m);
  rho_m = 0.0;
  eg_m = Vector_t(0.0);

  if (fs_m->hasValidSolver()) {

    this->Q *= this->dt;
    this->Q.scatter(this->rho_m, this->R, IntrplCIC_t());
    this->Q /= this->dt;
    this->rho_m /= getdT();

    //calculating mesh-scale factor
    double gammaz = sum(this->P)[2]/getTotalNum();
    gammaz *= gammaz; gammaz = sqrt(gammaz + 1.0);
    double scaleFactor = Physics::c * getdT();
    double tmp2 = 1/hr_m[0] * 1/hr_m[1] * 1/hr_m[2] / (scaleFactor*scaleFactor*scaleFactor) / gammaz;

    rho_m *= tmp2;

    Vector_t hr_scaled = hr_m * Vector_t(scaleFactor);
    hr_scaled[2] *= gammaz;

    // charge density is in rho_m

    fs_m->solver_m->computePotential(rho_m, hr_scaled);

    //scale mesh back
    rho_m *= hr_scaled[0]*hr_scaled[1]*hr_scaled[2];

    // the scalar potential is given back in rho_m
    // and must be converted to the right units
    rho_m *= getCouplingConstant();

    //write out rho
#ifdef DBG_SCALARFIELD
    *gmsg << "*** START DUMPING SCALAR FIELD ***" << endl;
    ostringstream oss;
    MPI_File file;
    MPI_Status status;
    MPI_Info fileinfo;
    MPI_File_open(MPI_COMM_WORLD, "rho_scalar", MPI_MODE_WRONLY | MPI_MODE_CREATE, fileinfo, &file);
    NDIndex<3> myidx = getFieldLayout().getLocalNDIndex();
    for(int x=myidx[0].first(); x<=myidx[0].last(); x++) {
      for(int y=myidx[1].first(); y<=myidx[1].last(); y++) {
        for(int z=myidx[2].first(); z<=myidx[2].last(); z++) {
          oss << x+1 << " " << y+1 << " " << z+1 << " " <<  rho_m[x][y][z].get() << endl;
        }
      }
    }
    MPI_File_write_shared(file, (char*)oss.str().c_str(), oss.str().length(), MPI_CHAR, &status);
    MPI_File_close(&file);
    *gmsg << "*** FINISHED DUMPING SCALAR FIELD ***" << endl;
#endif

    // IPPL Grad divides by hr_m
    eg_m = -Grad(rho_m,eg_m);

    eg_m *= Vector_t(gammaz/(scaleFactor), gammaz/(scaleFactor), 1.0/(scaleFactor*gammaz));

    // interpolate electric field at particle positions.  We reuse the
    // cached information about where the particles are relative to the
    // field, since the particles have not moved since this the most recent
    // scatter operation.
    Ef.gather(eg_m, this->R,  IntrplCIC_t());

    /** Magnetic field in x and y direction induced by the eletric field
     *
     *  \f[ B_x = \gamma(B_x^{'} - \frac{beta}{c}E_y^{'}) = -\gamma \frac{beta}{c}E_y^{'} = -\frac{beta}{c}E_y \f]
     *  \f[ B_y = \gamma(B_y^{'} - \frac{beta}{c}E_x^{'}) = +\gamma \frac{beta}{c}E_x^{'} = +\frac{beta}{c}E_x \f]
     *  \f[ B_z = B_z^{'} = 0 \f]
     *
     */
    double betaC = sqrt(gammaz*gammaz - 1.0)/gammaz/Physics::c;

    Bf(0) = Bf(0) - betaC * Ef(1);
    Bf(1) = Bf(1) + betaC * Ef(0);

  }
  IpplTimings::stopTimer(selfFieldTimer_m);
}

void PartBunch::computeSelfFields_cycl(double gamma)
{
  IpplTimings::startTimer(selfFieldTimer_m);

  /// set initial charge density to zero.
  rho_m = 0.0;

  /// set initial E field to zero
  eg_m = Vector_t(0.0);

  if (fs_m->hasValidSolver()) {

    /// scatter particles charge onto grid.
    this->Q.scatter(this->rho_m, this->R, IntrplCIC_t());

    /// from charge to charge density.
    double tmp2 = 1.0/ gamma/(hr_m[0] * hr_m[1] * hr_m[2]);
    rho_m *= tmp2;

    /// Lorentz transformation
    /// In particle rest frame, the longitudinal length enlarged
    Vector_t hr_scaled = hr_m ;
    hr_scaled[1] *= gamma;

    /// now charge density is in rho_m
    /// calculate Possion equation (without coefficient: -1/(eps))
    fs_m->solver_m->computePotential(rho_m, hr_scaled);

    /// additional work of FFT solver
    /// now the scalar potential is given back in rho_m
    rho_m *= hr_scaled[0]*hr_scaled[1]*hr_scaled[2];

    /// retrive coefficient: -1/(eps)
    rho_m *= getCouplingConstant();

    /// calculate electric field vectors from field potential
    eg_m = -Grad(rho_m,eg_m);

    /// back Lorentz transformation
    eg_m *= Vector_t(gamma, 1.0, gamma);

    /*
    //debug
    // output field on the grid points

    int m1 = (int)nr_m[0]-1;
    int m2 = (int)nr_m[0]/2;

    for (int i=0; i<m1; i++ )
    *gmsg << "Field along x axis E = " << eg_m[i][m2][m2] << " Pot = " << rho_m[i][m2][m2]  << endl;

    for (int i=0; i<m1; i++ )
    *gmsg << "Field along y axis E = " << eg_m[m2][i][m2] << " Pot = " << rho_m[m2][i][m2]  << endl;

    for (int i=0; i<m1; i++ )
    *gmsg << "Field along z axis E = " << eg_m[m2][m2][i] << " Pot = " << rho_m[m2][m2][i]  << endl;
    // end debug
    */

    /// interpolate electric field at particle positions.
    Ef.gather(eg_m, this->R,  IntrplCIC_t());

    /// calculate coefficient
    double betaC = sqrt(gamma*gamma - 1.0)/gamma/Physics::c;

    /// calculate B field from E field
    Bf(0) =  betaC * Ef(2);
    Bf(2) = -betaC * Ef(0);

  }
  // *gmsg<<"gamma ="<<gamma<<endl;
  // *gmsg<<"dx,dy,dz =("<<hr_m[0]<<", "<<hr_m[1]<<", "<<hr_m[2]<<") [m] "<<endl;
  // *gmsg<<"max of bunch is ("<<rmax_m(0)<<", "<<rmax_m(1)<<", "<<rmax_m(2)<<") [m] "<<endl;
  // *gmsg<<"min of bunch is ("<<rmin_m(0)<<", "<<rmin_m(1)<<", "<<rmin_m(2)<<") [m] "<<endl;
  IpplTimings::stopTimer(selfFieldTimer_m);
}



void PartBunch::computeSelfFields_cycl(int bin)
{
  IpplTimings::startTimer(selfFieldTimer_m);

  /// set initial charge dentsity to zero.
  rho_m = 0.0;

  /// set initial E field to zero
  eg_m = Vector_t(0.0);

  /// get gamma of this bin
  double gamma = getBinGamma(bin);

  if (fs_m->hasValidSolver()) {

    /// scatter particles charge onto grid.
    this->Q.scatter(this->rho_m, this->R, IntrplCIC_t());

    /// from charge to charge density.
    double tmp2 = 1.0/ gamma/(hr_m[0] * hr_m[1] * hr_m[2]);
    rho_m *= tmp2;

    /// Lorentz transformation
    /// In particle rest frame, the longitudinal length enlarged
    Vector_t hr_scaled = hr_m ;
    hr_scaled[1] *= gamma;

    /// now charge density is in rho_m
    /// calculate Possion equation (without coefficient: -1/(eps))
    fs_m->solver_m->computePotential(rho_m, hr_scaled);

    /// additional work of FFT solver
    /// now the scalar potential is given back in rho_m
    rho_m *= hr_scaled[0]*hr_scaled[1]*hr_scaled[2];

    /// retrive coefficient: -1/(eps)
    rho_m *= getCouplingConstant();

    /// calculate electric field vectors from field potential
    eg_m = -Grad(rho_m,eg_m);

    /// back Lorentz transformation
    eg_m *= Vector_t(gamma, 1.0, gamma);

    /*
    //debug
    // output field on the grid points

    int m1 = (int)nr_m[0]-1;
    int m2 = (int)nr_m[0]/2;

    for (int i=0; i<m1; i++ )
    *gmsg << "Field along x axis E = " << eg_m[i][m2][m2] << " Pot = " << rho_m[i][m2][m2]  << endl;

    for (int i=0; i<m1; i++ )
    *gmsg << "Field along y axis E = " << eg_m[m2][i][m2] << " Pot = " << rho_m[m2][i][m2]  << endl;

    for (int i=0; i<m1; i++ )
    *gmsg << "Field along z axis E = " << eg_m[m2][m2][i] << " Pot = " << rho_m[m2][m2][i]  << endl;
    // end debug
    */

    /// interpolate electric field at particle positions.
    Eftmp.gather(eg_m, this->R,  IntrplCIC_t());

    /// calculate coefficient
    double betaC = sqrt(gamma*gamma - 1.0)/gamma/Physics::c;

    /// calculate B_bin field from E_bin field accumulate B and E field
    Bf(0) = Bf(0) + betaC * Eftmp(2);
    Bf(2) = Bf(2) - betaC * Eftmp(0);

    Ef += Eftmp;
  }
  // *gmsg<<"gamma ="<<gamma<<endl;
  // *gmsg<<"dx,dy,dz =("<<hr_m[0]<<", "<<hr_m[1]<<", "<<hr_m[2]<<") [m] "<<endl;
  // *gmsg<<"max of bunch is ("<<rmax_m(0)<<", "<<rmax_m(1)<<", "<<rmax_m(2)<<") [m] "<<endl;
  // *gmsg<<"min of bunch is ("<<rmin_m(0)<<", "<<rmin_m(1)<<", "<<rmin_m(2)<<") [m] "<<endl;
  IpplTimings::stopTimer(selfFieldTimer_m);
}

void PartBunch::boundp() {
  /*
    Assume rmin_m < 0.0
  */
  Vector_t len;
  const int dimIdx = 3;
  IpplTimings::startTimer(boundpTimer_m);

  //  if (hasZeroNLP())
  //    BinaryRepartition(*this);

  NDIndex<3> domain = getFieldLayout().getDomain();
  for(int i=0; i<Dim; i++)
    nr_m[i] = domain[i].length();

  get_bounds(rmin_m, rmax_m);
  len = rmax_m - rmin_m;

  for(int i=0; i<dimIdx; i++) {
    rmax_m[i] += dh_m*abs(rmax_m[i]-rmin_m[i]);
    rmin_m[i] -= dh_m*abs(rmax_m[i]-rmin_m[i]);
    hr_m[i]    = (rmax_m[i]-rmin_m[i]) / (nr_m[i]-1);
  }

  // rescale mesh
  getMesh().set_meshSpacing(&(hr_m[0]));
  getMesh().set_origin( rmin_m - Vector_t(hr_m[0] / 2.0, hr_m[1] / 2.0, hr_m[2] / 2.0) );

  rho_m.initialize(getMesh(),
                   getFieldLayout(),
                   GuardCellSizes<Dim>(1),
                   bc_m);
  eg_m.initialize(getMesh(),
                  getFieldLayout(),
                  GuardCellSizes<Dim>(1),
                  vbc_m);
  update();
  IpplTimings::stopTimer(boundpTimer_m);
}

void PartBunch::boundpNoRep() {
  /*
    Assume rmin_m < 0.0
  */
  Vector_t len;
  const int dimIdx = 3;
  IpplTimings::startTimer(boundpTimer_m);

  NDIndex<3> domain = getFieldLayout().getDomain();
  for(int i=0; i<Dim; i++)
    nr_m[i] = domain[i].length();

  get_bounds(rmin_m, rmax_m);
  len = rmax_m - rmin_m;

  for(int i=0; i<dimIdx; i++) {
    rmax_m[i] += dh_m*abs(rmax_m[i]-rmin_m[i]);
    rmin_m[i] -= dh_m*abs(rmax_m[i]-rmin_m[i]);
    hr_m[i]    = (rmax_m[i]-rmin_m[i]) / (nr_m[i]-1);
  }

  // rescale mesh
  getMesh().set_meshSpacing(&(hr_m[0]));
  getMesh().set_origin( rmin_m );

  rho_m.initialize(getMesh(),
                   getFieldLayout(),
                   GuardCellSizes<Dim>(1),
                   bc_m);
  eg_m.initialize(getMesh(),
                  getFieldLayout(),
                  GuardCellSizes<Dim>(1),
                  vbc_m);
  update();
  IpplTimings::stopTimer(boundpTimer_m);
}

void PartBunch::calcWeightedAverages(Vector_t &CentroidPosition, Vector_t &CentroidMomentum) const
{
  double gamma;
  double cent[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  const double N =  static_cast<double>(this->getTotalNum());

  for (int i = 0; i < this->getLocalNum(); i++)
    {
      gamma = sqrt(1.0 + dot(this->P[i],this->P[i]));
      cent[0] += this->R[i](0);
      cent[1] += this->R[i](1);
      cent[2] += this->R[i](2);
      cent[3] += this->P[i](0) / gamma;
      cent[4] += this->P[i](1) / gamma;
      cent[5] += this->P[i](2) / gamma;

    }
  reduce(&(cent[0]),&(cent[0]) + 6, &(cent[0]), OpAddAssign());

  CentroidPosition(0) = cent[0] / N;
  CentroidPosition(1) = cent[1] / N;
  CentroidPosition(2) = cent[2] / N;
  CentroidMomentum(0) = cent[3] / N;
  CentroidMomentum(1) = cent[4] / N;
  CentroidMomentum(2) = cent[5] / N;
}

void PartBunch::rotateAbout(const Vector_t &Center, const Vector_t &EulerAngles)
{
  Vector_t tmpR, tmpP;
  double cos0 = cos(EulerAngles(0)),
    sin0 = sin(EulerAngles(0)),
    cos1 = cos(EulerAngles(1)),
    sin1 = sin(EulerAngles(1)),
    cos2 = cos(EulerAngles(2)),
    sin2 = sin(EulerAngles(2));
  double T00 = cos2 * cos0 - sin2 * cos1 * sin0,
    T01 = - sin2 * sin1,
    T02 = cos2 * sin0 + sin2 * cos1 * cos0,
    T10 = -sin1 * sin0,
    T11 = cos1,
    T12 = sin1 * cos0,
    T20 = -sin2 * cos0 - cos2 * cos1 * sin0,
    T21 = -cos2 * sin1,
    T22 = -sin2 * sin0 + cos2 * cos1 * cos0;

  for (int i = 0; i < this->getLocalNum(); i++)
    {
      tmpR = Vector_t(this->R[i](0) - Center(0), this->R[i](1) - Center(1), this->R[i](2) - Center(2));
      tmpP = this->P[i];
      R[i](0) = T00 * tmpR(0) + T01 * tmpR(1) + T02 * tmpR(2) + Center(0);
      R[i](1) = T10 * tmpR(0) + T11 * tmpR(1) + T12 * tmpR(2) + Center(1);
      R[i](2) = T20 * tmpR(0) + T21 * tmpR(1) + T22 * tmpR(2) + Center(2);
      P[i](0) = T00 * tmpP(0) + T01 * tmpP(1) + T02 * tmpP(2);
      P[i](1) = T10 * tmpP(0) + T11 * tmpP(1) + T12 * tmpP(2);
      P[i](2) = T20 * tmpP(0) + T21 * tmpP(1) + T22 * tmpP(2);
    }
}

void PartBunch::moveBy(const Vector_t &Center)
{
  for (int i = 0; i < this->getLocalNum(); i++)
    {
      R[i] += Center;
    }
}

void PartBunch::ResetLocalCoordinateSystem(const int &i, const Vector_t &Orientation, const double &origin)
{
  Vector_t temp(R[i](0), R[i](1), R[i](2) - origin);
  if (fabs(Orientation(0)) > 1e-6 || fabs(Orientation(1)) > 1e-6)
    {
      double cosa = cos(Orientation(0));
      double sina = sin(Orientation(0));
      double cosb = cos(Orientation(1));
      double sinb = sin(Orientation(1));

      X[i](0) = cosa *         temp(0)                  - sina *        temp(2);
      X[i](1) = -sina * sinb * temp(0) + cosb * temp(1) - cosa * sinb * temp(2);
      X[i](2) = sina * cosb *  temp(0) + sinb * temp(1) + cosa * cosb * temp(2);
    }
  else
    {
      X[i] = temp;
    }
}


void PartBunch::beamEllipsoid(FVector<double,6>   &centroid,
                              FMatrix<double,6,6> &moment)
{
  for (int i = 0; i < 6; ++i) {
    centroid(i) = 0.0;
    for (int j = 0; j <= i; ++j) {
      moment(i,j) = 0.0;
    }
  }

  //  PartBunch::const_iterator last = end();
  // for (PartBunch::const_iterator part = begin(); part != last; ++part) {

  Particle part;

  for (int ii = 0; ii < this->getLocalNum(); ii++) {
    part = get_part(ii);
    for (int i = 0; i < 6; ++i) {
      centroid(i) += part[i];
      for (int j = 0; j <= i; ++j) {
        moment(i,j) += part[i] * part[j];
      }
    }
  }

  double factor = 1.0 / double(this->getTotalNum());
  for (int i = 0; i < 6; ++i) {
    centroid(i) *= factor;
    for (int j = 0; j <= i; ++j) {
      moment(j,i) = moment(i,j) *= factor;
    }
  }
}


void PartBunch::gatherLoadBalanceStatistics()
{
  for (int i=0; i<Ippl::getNodes(); i++) 
    partPerNode_m[i] = globalPartPerNode_m[i] = 0.0;

  partPerNode_m[Ippl::myNode()] = this->getLocalNum();

  reduce(partPerNode_m,partPerNode_m+Ippl::getNodes(),globalPartPerNode_m,OpAddAssign());

}

void PartBunch::calcMoments() {

  double part[2*Dim];
  double loc_centroid[2*Dim];
  double loc_moment[2*Dim][2*Dim];
  double moments[2*Dim][2*Dim];

  for (int i = 0; i < 2*Dim; ++i) {
    loc_centroid[i] = 0.0;
    for (int j = 0; j <= i; ++j) {
      loc_moment[i][j] = 0.0;
    }
  }

  for (unsigned long k=0; k< this->getLocalNum(); ++k){
    part[1] = this->P[k](0);
    part[3] = this->P[k](1);
    part[5] = this->P[k](2);
    part[0] = this->R[k](0);
    part[2] = this->R[k](1);
    part[4] = this->R[k](2);

    for (int i = 0; i < 2*Dim; ++i) {
      loc_centroid[i] += part[i];
      for (int j = 0; j <= i; ++j) {
        loc_moment[i][j] += part[i] * part[j];
      }
    }
  }

  for (int i = 0; i < 2*Dim; ++i) {
    for (int j = 0; j < i; ++j) {
      loc_moment[j][i] = loc_moment[i][j];
    }
  }

  reduce(&(loc_moment[0][0]), &(loc_moment[0][0]) + 2*Dim*2*Dim,
         &(moments[0][0]), OpAddAssign());

  reduce(&(loc_centroid[0]), &(loc_centroid[0]) + 2*Dim,
         &(centroid_m[0]), OpAddAssign());

  for (int i = 0; i < 2*Dim; ++i) {
    for (int j = 0; j <= i; ++j) {
      moments_m(i,j) = moments[i][j];
      moments_m(j,i) = moments_m(i,j);
    }
  }
}

void PartBunch::calcMomentsInitial() {

  double part[2*Dim];

  for (int i = 0; i < 2*Dim; ++i) {
    centroid_m[i] = 0.0;
    for (int j = 0; j <= i; ++j) {
      moments_m(i,j) = 0.0;
    }
  }

  for (size_t k = 0; k < pbin_m->getNp(); k++) {
    for (int binNumber = 0; binNumber < pbin_m->getNBins(); binNumber++) {
      vector<double> p;

      if (pbin_m->getPart(k, binNumber, p)) {
        part[0] = p.at(0);
        part[1] = p.at(3);
        part[2] = p.at(1);
        part[3] = p.at(4);
        part[4] = p.at(2);
        part[5] = p.at(5);

        for (int i = 0; i < 2*Dim; ++i) {
          centroid_m[i] += part[i];
          for (int j = 0; j <= i; ++j) {
            moments_m(i,j) += part[i] * part[j];
          }
        }
      }
    }
  }

  for (int i = 0; i < 2*Dim; ++i) {
    for (int j = 0; j < i; ++j) {
      moments_m(j,i) = moments_m(i,j);
    }
  }
}

void PartBunch::calcBeamParameters()
{
  using Physics::c;
  using Physics::m_e;

  Vector_t eps2,fac,rsqsum,psqsum,rpsum;

  IpplTimings::startTimer(statParamTimer_m);

  const size_t locNp = this->getLocalNum();

  const double zero = 0.0;
  const double N =  static_cast<double>(this->getTotalNum());
  if (N == 0)
    {
      for (unsigned int i=0 ; i<Dim; i++) {
        rmean_m(i) = 0.0;
        pmean_m(i) = 0.0;
        rrms_m(i) = 0.0;
        prms_m(i) = 0.0;
        eps_norm_m(i)  = 0.0;
      }
      rprms_m = 0.0;
      mean_energy_m = 0.0;
      eps_m = 0.0;
      IpplTimings::stopTimer(statParamTimer_m);
      return;
    }
  calcMoments();

  for (unsigned int i=0 ; i<Dim; i++) {
    rmean_m(i) = centroid_m[2*i]/N;
    pmean_m(i) = centroid_m[(2*i)+1]/N;
    rsqsum(i) = moments_m(2*i,2*i)-N*rmean_m(i)*rmean_m(i);
    psqsum(i) = moments_m((2*i)+1,(2*i)+1)-N*pmean_m(i)*pmean_m(i);
    if(psqsum(i) < 0)
      psqsum(i) = 0;
    rpsum(i) = moments_m((2*i),(2*i)+1) - N*rmean_m(i)*pmean_m(i);
  }
  eps2 = (rsqsum * psqsum - rpsum * rpsum)/(N*N);
  rpsum /= N;

  for (unsigned int i=0 ; i<Dim; i++) {
    rrms_m(i) = sqrt( rsqsum(i)/N );
    prms_m(i) = sqrt( psqsum(i)/N );
    eps_norm_m(i)  = sqrt( max( eps2(i), zero ) );
    double tmp = rrms_m(i) * prms_m(i);
    fac(i) = (tmp == 0) ? zero : 1.0 / tmp;
  }
  rprms_m = rpsum * fac;

  /*
    double rmax = sqrt(dot(rmax_m,rmax_m));
    fplasma_m = sqrt(2.0*get_perverance())*get_beta()*c/rmax;
    budkerp_m = (get_perverance()/2.0)*pow(get_gamma(),3.0)*pow(get_beta(),2.0);
  */

  // calculate mean energy

  mean_energy_m = 0.0;
  for(int k=0; k < locNp; k++)
    mean_energy_m += (sqrt(dot(P[k], P[k]) + 1.0) - 1.0) * m_e * 1e6 / 1.0e3;
  reduce(mean_energy_m,mean_energy_m,OpAddAssign());
  mean_energy_m /= N;

  // Find unnormalized emittance.
  double betagamma = 0.0;
  for(size_t i=0; i<locNp;i++)
    betagamma += sqrt(1.0 + dot(P[i],P[i]));
  reduce(betagamma,betagamma,OpAddAssign());
  betagamma /= N;
  betagamma *= sqrt(1.0 - (1/betagamma)*(1/betagamma));
  eps_m = eps_norm_m / Vector_t(betagamma);

  IpplTimings::stopTimer(statParamTimer_m);
}

void PartBunch::calcBeamParametersInitial()
{
  using Physics::c;
  using Physics::m_e;

  if (Ippl::myNode() == 0) {
    // Only use one node as this function will get called only once before
    // particles have been emitted (at least in principle).
    Vector_t eps2,fac,rsqsum,psqsum,rpsum;

    const double zero = 0.0;
    const double  N =  static_cast<double>(pbin_m->getNp());

    calcMomentsInitial();

    for (unsigned int i=0 ; i<Dim; i++) {
      rmean_m(i) = centroid_m[2*i]/N;
      pmean_m(i) = centroid_m[(2*i)+1]/N;
      rsqsum(i) = moments_m(2*i,2*i)-N*rmean_m(i)*rmean_m(i);
      psqsum(i) = moments_m((2*i)+1,(2*i)+1)-N*pmean_m(i)*pmean_m(i);
      if(psqsum(i) < 0)
	psqsum(i) = 0;
      rpsum(i) =  moments_m((2*i),(2*i)+1)-N*rmean_m(i)*pmean_m(i);
    }
    eps2 = (rsqsum * psqsum - rpsum * rpsum)/(N*N);
    rpsum /= N;

    for (unsigned int i=0 ; i<Dim; i++) {

      rrms_m(i) = sqrt( rsqsum(i)/N );
      prms_m(i) = sqrt( psqsum(i)/N );
      eps_m(i)  = sqrt( max( eps2(i), zero ) );
      double tmp = rrms_m(i) * prms_m(i);
      fac(i) = (tmp == 0) ? zero : 1.0 / tmp;
    }
    rprms_m = rpsum * fac;
  }
}

void PartBunch::setSolver(FieldSolver *fs) {
  fs_m = fs;
  fs_m->initSolver(*this);
  /**
     CAN not re-inizialize ParticleLayout
     this is an IPPL issue
  */
  if (!OPAL.hasBunchAllocated())
    initialize(&(fs_m->getParticleLayout()));
}

void PartBunch::maximumAmplitudes(const FMatrix<double,6,6> &D,
                                  double &axmax, double &aymax)
{
  axmax = aymax = 0.0;
  Particle part;

  for (int ii = 0; ii < getLocalNum(); ii++) {

    part = get_part(ii);

    double xnor =
      D(0,0)*part.x()  + D(0,1)*part.px() + D(0,2)*part.y() +
      D(0,3)*part.py() + D(0,4)*part.t()  + D(0,5)*part.pt();
    double pxnor =
      D(1,0)*part.x()  + D(1,1)*part.px() + D(1,2)*part.y() +
      D(1,3)*part.py() + D(1,4)*part.t()  + D(1,5)*part.pt();
    double ynor =
      D(2,0)*part.x()  + D(2,1)*part.px() + D(2,2)*part.y() +
      D(2,3)*part.py() + D(2,4)*part.t()  + D(2,5)*part.pt();
    double pynor =
      D(3,0)*part.x()  + D(3,1)*part.px() + D(3,2)*part.y() +
      D(3,3)*part.py() + D(3,4)*part.t()  + D(3,5)*part.pt();

    axmax = std::max(axmax, (xnor*xnor + pxnor*pxnor));
    aymax = std::max(aymax, (ynor*ynor + pynor*pynor));
  }
}

/**
 * Here we emit particles from the cathode. All particles in a new simulation (not a restart) initially reside in the bin
 container "pbin_m" and are not part of the beam bunch (so they cannot "see" fields, space charge etc.). In pbin_m, particles
 drift in a field free region with their initial properties unchanged until they reach the cathode position (z = 0.0). Once
 they cross the cathode position, particles are moved from pbin_m into the beam bunch "PartBunch," where there are affected
 by fields, space charge and so on.

 A typical integration time step, \f$\Delta t\f$, is broken down into 3 sub-steps:

    1) Drift particles for \f$\frac{\Delta t}{2}\f$.

    2) Calculate fields and advance momentum.

    3) Drift particles for \f$\frac{\Delta t}{2}\f$ at the new momentum to complete the
       full time step.

 The difficulty for emission is that at the cathode position there is a step function discontinuity in the  fields. If we
 apply the typical integration time step across this boundary, we get an artificial numerical bunching of the beam, especilly
 at very high accelerating fields. This function takes the cathode position boundary into account in order to achieve
 smoother particle emission.

 Particles can cross the cathode boundary in either step 1) or step 3) from above. If they cross in step 1), we break the total
 integration step into two parts:
 \f[
 \Delta t_{full-timesetp} = \Delta t_{emit} + \Delta t_{temp}
\f]
 \f$\Delta t_{emit}\f$ is the time it takes the particle to reach the cathode position at its intial velocity. So, we drift the
 particle to the cathode position and emit it into PartBunch and set its time step to the temporary value \f$\Delta t_{temp}\f$.
 The particle is then integrated with this reduced time step during its first step inside the real simulation space. For,
 subsequent integrations its time step is reset back to \f$\Delta t_{full-timestep}\f$.

 If particles cross the cathode boundary in step 3), we again break the total integration step into two parts:
\f[
\Delta t_{full-timestep} = -\Delta t_{emit} + \Delta t_{temp}
\f]
 Now \f$\Delta t_{emit}\f$ is equal to how long ago the particle passed the cathode position moving at its initial velocity. In
 this case, we drift the particle back to the cathode position and set the time step to the temporary value \f$\Delta t_{temp}\f$.
 The particle is then integrated with this increased time step during its first step inside the real simulation space. For,
 subsequent integrations its time step is reset back to \f$\Delta t_{full-timestep}\f$.
*/
size_t PartBunch::emitParticles(int binNumber) {

  /// Get number of particles that have been emitted and are on this processor.
  size_t numberOfEmittedParticles = this->getLocalNum();
  size_t oldNumberOfEmittedParticles = numberOfEmittedParticles;

  /// Check if bin has already been emitted. Make sure it is a legal bin number.
  if (!pbin_m->getBinEmitted(binNumber) && binNumber < pbin_m->getNBins()) {

    /// Loop over particles in bin.
    for (size_t particleNumber = 0; particleNumber < pbin_m->getNp(); particleNumber++) {

      /// Get particles.
      vector<double> particleCoordinates;
      if (pbin_m->getPart(particleNumber, binNumber, particleCoordinates)) {
	if (!pbin_m->isEmitted(particleNumber, binNumber)) {

	  /// If beam has not been emitted proceed.
	  double particleGamma = sqrt(1.0 + pow(particleCoordinates.at(3), 2.0) + pow(particleCoordinates.at(4), 2.0) + pow(particleCoordinates.at(5), 2.0));
	  double particleSVelocity = particleCoordinates.at(5) * Physics::c / particleGamma;
	  double sDrift = particleSVelocity * getdT() / 2.0;

	  if (particleCoordinates.at(2) + sDrift >= 0.0) {
	    /// Particle is emitted. Update position and change time step appropriately.

	    create(1); /// Emit particle into PartBunch.

	    /** Find the time when the particle was emitted. This will either be in the future (particle emitted in step 1 of the full
		integration step) or in the past (particle emitted in step 3 of the full integration step). So, if we define
		\f[
		t_{emit} = \frac{z_{particle}}{v_{particle}}
		\f]
		then we change the time step for the particle to
		\f[
		\Delta t_{particle} = \Delta t_{simulation} + t_{emit}
		\f]
		This time step will be changed back to its original value after the first step after emission.

		This call is made after the first half step in the general integration step for PartBunch in the beam tracker.
		So, we end by advancing the particles one half of their new, temporary time step.

		Also note that at this point in the general integration algorithm the positions should be divided by /f$c /Delta t/f$.
	    */

	    this->dt[numberOfEmittedParticles] = particleCoordinates.at(2) / particleSVelocity + getdT();
	    //	    this->dt[numberOfEmittedParticles] = getdT();

	    double oneOverCDt = 1.0 / (Physics::c * this->dt[numberOfEmittedParticles]);

	    this->R[numberOfEmittedParticles] = Vector_t(oneOverCDt * (particleCoordinates.at(0)
								       + particleCoordinates.at(3) * this->dt[numberOfEmittedParticles] * Physics::c / (2.0 * particleGamma)),
							 oneOverCDt * (particleCoordinates.at(1)
								       + particleCoordinates.at(4) * this->dt[numberOfEmittedParticles] * Physics::c / (2.0 * particleGamma)),
							 //							 oneOverCDt * (particleCoordinates.at(2)
							 //								       + particleCoordinates.at(5) * this->dt[numberOfEmittedParticles] * Physics::c / (2.0 * particleGamma)));
							 oneOverCDt * (particleCoordinates.at(5) * this->dt[numberOfEmittedParticles] * Physics::c / (2.0 * particleGamma)));

	    /// Fill in the rest of the particle propeties.
	    this->P[numberOfEmittedParticles] = Vector_t(particleCoordinates.at(3),
							 particleCoordinates.at(4),
							 particleCoordinates.at(5));
	    this->Bin[numberOfEmittedParticles] = binNumber;
	    this->Q[numberOfEmittedParticles] = getChargePerParticle();
	    this->LastSection[numberOfEmittedParticles] = 0;
	    this->Ef[numberOfEmittedParticles] = Vector_t(0.0); /// Initialize fields to zero.
	    this->Bf[numberOfEmittedParticles] = Vector_t(0.0);

	    pbin_m->setEmitted(particleNumber, binNumber);     /// Set particle as emitted.
	    numberOfEmittedParticles++;                        /// Iterate number of particles that have been emitted.
	    binemitted_m[binNumber]++;                         /// Iterate number of particles in this bin that have been emitted.

	  } else {
	    /// Particle not emitted. Drift with initial parameters.
	    pbin_m->updatePartPos(particleNumber, binNumber, particleCoordinates.at(2) + 2.0 * sDrift);
	  }
	}
      }
    }

    /// Check if all particles in bin have been emitted. Set bin to emitted if so.
    if (binemitted_m[binNumber] == pbin_m->getBinCont(binNumber)) {
      pbin_m->setBinEmitted(binNumber);
      *gmsg << "* Bin number: " << binNumber << " has emitted all particles." << endl;
    }

    /// Set this as the last bin to have emitted.
    if ((numberOfEmittedParticles - oldNumberOfEmittedParticles) != 0) {
      pbin_m->setActualemittedBin(binNumber);
    }
  }

  /// Return number of particles emitted.
  return numberOfEmittedParticles - oldNumberOfEmittedParticles;
}

/** When a particle crosses
where they become part
 of the simulation in that theyat which point they areonly after they crosssimulation space (i.e. inside the beam line  that would step over the cathode (z=0) in the first half step of the integration scheme.
 This particles are added to the bunch and reset their timestep to the remaining distance they would cover until after the second half step diveded by two.
 The timestep for these particles will be shorter in this step than the full timestep of the simulation.
 The timestep of emitted particles is reset to the simluation timestep after the full timestep has been applied.
 \f[
 \Delta t = \frac{P_z}{\beta c} + \frac{\Delta t_{full-timestep}}{2}
 \f]
 The particles that do not step over the cathode update their z position.
*/
// size_t PartBunch::emitParticles(int bin) {
//   /**
//      We copy all particles in the given
//      bin from pbin_m to the particle
//      container.
//   */
//   size_t nloc = this->getLocalNum();
//   size_t oldtot = nloc;

//   if(!pbin_m->getBinEmitted(bin) && bin < pbin_m->getNBins()) {

//     //*gmsg << "* ************** E M I T *********************************************************** " << endl;
//     //*gmsg << "* BIN= " << bin << " out of " << pbin_m->getNBins() << endl;
//     //  *gmsg << "* rmin_i= " << rmin_m << " rmax_i= " << rmax_m << " h_i= " << hr_m << endl;

//     double bingamma = 0.0;
//     int nonemittedp = 0;
//     //FIXME: DONT DO EVERY EMIT CALL!
//     for (size_t i=0; i<pbin_m->getNp(); i++) {
//       vector<double> p;
//       if (pbin_m->getPart(i,bin,p)) {
//         //         if( !pbin_m->isEmitted(i,bin) ) {
//         bingamma += sqrt(1.0 + p[3]*p[3] + p[4]*p[4] + p[5]*p[5]);
//         //           ++nonemittedp;
//         //         }
//       }
//     }
//     //     bingamma /= nonemittedp;
//     bingamma /= pbin_m->getBinCont(bin);
//     const double binbetac = sqrt(1.0 - (1.0/(bingamma*bingamma))) * Physics::c;
//     const double tempz = getdT()/2.0 * binbetac;
//     const double recpcdt = 1.0 / (Physics::c * getdT());

//     for (size_t i=0; i<pbin_m->getNp(); i++) {
//       vector<double> p;
//       if (pbin_m->getPart(i,bin,p)) {
//         if( !pbin_m->isEmitted(i,bin) ) {
// 	  const double recpgamma = 1.0 / sqrt(1.0 + p[3]*p[3] + p[4]*p[4] + p[5]*p[5]);
// 	  double betaTemp = p[5] * recpgamma;
// 	  double vTemp = betaTemp * Physics::c;
// 	  double zTemp = vTemp * getdT() / 2.0;
//           if(p[2] + tempz > 0.0) {
//             create(1);
//             /** calculate the time when the particle is emitted, \f$t_{emit}\f$, and subtract this from \f$dT\f$ to get the proper \f$dt\f$ for the particle.
//                 Set then the particle back to the surface of the cathode and push it from there to the right position with its own \f$dt\f$.
//                 \f[
//                 t_{emit} = -\frac{z}{\beta c} \; (z<0)\\
//                 dt = dT - t_{emit} = dT + \frac{z}{\beta c} \\
//                 z' = \frac{dt}{2} \beta c
//                 \f]
//                 The position should be dimensionless at this stage of the algorithm therefore we devide by c*dt.
//             */

// 	    this->dt[nloc] = p[2]/(binbetac) + getdT();
// 	    //            this->R[nloc] = Vector_t(p[0]*recpcdt, p[1]*recpcdt, 0.25 * p[5] * recpgamma * (p[2]/tempz + 2.0));
//             this->R[nloc] = Vector_t(p[0]*recpcdt, p[1]*recpcdt, 0.25 * p[5] * recpgamma * (p[2]/zTemp + 2.0));

//             this->P[nloc] = Vector_t(0.0,0.0,p[5]);
//             this->Bin[nloc]=bin;
//             this->Q[nloc] = getChargePerParticle();
//             this->LastSection[nloc]=0;
//             this->Ef[nloc] = Vector_t(0.0);
//             this->Bf[nloc] = Vector_t(0.0);


//             //             cerr << this->R[nloc](0)/recpcdt << "\t"
//             //                  << this->R[nloc](1)/recpcdt << "\t"
//             //                  << this->R[nloc](2)/recpcdt << endl;

//             pbin_m->setEmitted(i,bin);
//             nloc++;
//             binemitted_m[bin]++;
//           } else {
//             //update z position of particle 'i' in bin 'bin'
//             pbin_m->updatePartPos(i, bin, p[2] + 2. * tempz);
//           }
//         }
//       }
//     }

//     if(binemitted_m[bin] == pbin_m->getBinCont(bin)) {
//       pbin_m->setBinEmitted(bin);
//       *gmsg << "* Bin " << bin << " has emitted all particles" << endl;
//     }
//     /*
//        Because only node 0 is doing that we can not
//        do an boundp here!
//        boundp();
//     */

//     if ((nloc-oldtot) != 0) {
//       pbin_m->setActualemittedBin(bin);
//     }
//   }
//   return (nloc-oldtot);
// }

double PartBunch::calcTimeDelay(const double &jifactor)
{
  double gamma = pbin_m->getGamma();
  double beta = sqrt(1. - (1. / (gamma*gamma)));
  double xmin, xmax;
  pbin_m->getExtrema(xmin,xmax);
  double length = xmax - xmin;

  return length * jifactor / (Physics::c * beta);
}

void PartBunch::moveBunchToCathode(double &t)
{
  double avrg_betagamma = 0.0;
  double maxspos = -9999999.99;

  for (int bin = 0; bin < getNumBins(); ++bin)
    {
      for (size_t i = 0; i < pbin_m->getNp(); ++i)
        {
          vector<double> p;
          if (pbin_m->getPart(i,bin,p))
            {
              avrg_betagamma += sqrt(1.0 + p[3]*p[3] + p[4]*p[4] + p[5]*p[5]);
              if (p[2] > maxspos) maxspos = p[2];
            }
        }
    }
  avrg_betagamma /= pbin_m->getNp();
  double dist_per_step = sqrt(1.0 - (1.0/(avrg_betagamma * avrg_betagamma))) * Physics::c * getdT();
  if (maxspos < 0.0)
    {
      double num_steps = floor(-maxspos/dist_per_step);
      t += num_steps * getdT();

      Inform gmsg("PartBunch");
      gmsg << "move bunch by " << num_steps*dist_per_step << "; DT = " << num_steps * getdT() << endl;
      for (int bin = 0; bin < getNumBins(); ++bin)
        {
          for (size_t i = 0; i < pbin_m->getNp(); ++i)
            {
              vector<double> p;
              if (pbin_m->getPart(i,bin,p))
                {
                  pbin_m->updatePartPos(i,bin, p[2] + num_steps * dist_per_step);
                }
            }
        }
    }
}

void PartBunch::printBinHist()
{
  if (weHaveBins())
    {
      int *binhisto = new int[getNumBins()];
      double maxz = -999999999.99, minz = 999999999.99;


      for (int bin = 0; bin < getNumBins(); ++bin)
        {
          binhisto[bin] = 0;
        }
      for (size_t i = 0; i < getLocalNum(); ++i)
        ++binhisto[this->Bin[i]];
      reduce(&(binhisto[0]), &(binhisto[0]) + getNumBins(), &(binhisto[0]), OpAddAssign());


      if (Ippl::myNode() == 0)
        {
          for (int bin = 0; bin < getNumBins(); ++bin)
            {
              for (size_t i = 0; i < pbin_m->getNp(); ++i)
                {
                  vector<double> p;
                  if (pbin_m->getPart(i,bin,p))
                    {
                      if (p[2] > maxz) maxz = p[2];
                      if (p[2] < minz) minz = p[2];
                    }
                }
            }
          double dz = (maxz - minz) / 20;
          int minihist[20] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
          for (int bin = 0; bin < getNumBins(); ++bin)
            {
              for (size_t i = 0; i < pbin_m->getNp(); ++i)
                {
                  vector<double> p;
                  if (pbin_m->getPart(i,bin,p))
                    {
                      ++minihist[(int)floor((p[2] - minz) / dz)];
                    }
                }
            }
          cout << "particles are between " << minz << " and " << maxz << endl;
          for (int i = 0; i < 20; ++i)
            {
              cout << i << ": " << minihist[i] << endl;
            }

          cout << "\n\n"
               << "effective histogram:" << endl;;
          for (int bin = 0; bin < getNumBins(); ++bin)
            cout << bin << ": " << binhisto[bin] << endl;

          delete[] binhisto;
        }
    }
}


Inform &PartBunch::print(Inform &os)
{
  if (this->getTotalNum() != 0) { // to suppress Nan's
    os << "* ************** B U N C H ********************************************************* " << endl;
    os << "* NP= " << this->getTotalNum() << " Qtot= " << abs(sum(Q)*1.0E9) << " [nC]  Qi= " << abs(qi_m) << " [C] " << endl;
    //  INFOMSG(getMesh() << endl);
    os << "* rmax= " << rmax_m << endl;
    os << "* rmin= " << rmin_m << endl;
    os << "* rms beam size= " << rrms_m << endl;
    os << "* rms momenta= " << prms_m << endl;
    os << "* mean beam size= " << rmean_m << endl;
    os << "* mean momenta= " << pmean_m << endl;
    os << "* rms emmitance= " << eps_m << " (not normalized)" << endl;
    os << "* rms correlation= " << rprms_m << endl;
    os << "* dh = " << dh_m << "\t dT= " << getdT() << endl;
    os << "* hr = " << hr_m << endl;
    os << "* tEmission= " << getTEmission() << endl;
    os << "* spos= " << get_sPos() << endl;
    os << "* ********************************************************************************** " << endl;
  }
  return os;
}

void PartBunch::calcBeamParameters_cycl()
{
  using Physics::c;
  using Physics::m_p;

  Vector_t eps2,fac,rsqsum,psqsum,rpsum;

  const size_t locNp = this->getLocalNum();
  double localMeanEnergy = 0.0;

  const double zero = 0.0;
  const double TotalNp =  static_cast<double>(this->getTotalNum());

  // calculate centroid_m and moments_m
  calcMoments();

  for (unsigned int i=0 ; i<Dim; i++) {
    rmean_m(i) = centroid_m[2*i]/TotalNp;
    pmean_m(i) = centroid_m[(2*i)+1]/TotalNp;
    rsqsum(i) = moments_m(2*i,2*i)-TotalNp*rmean_m(i)*rmean_m(i);
    psqsum(i) = moments_m((2*i)+1,(2*i)+1)-TotalNp*pmean_m(i)*pmean_m(i);
    rpsum(i) =  moments_m((2*i),(2*i)+1)-TotalNp*rmean_m(i)*pmean_m(i);
  }
  eps2      = (rsqsum * psqsum - rpsum * rpsum)/(TotalNp*TotalNp);
  rpsum /= TotalNp;

  for (unsigned int i=0 ; i<Dim; i++) {
    rrms_m(i) = sqrt( rsqsum(i)/TotalNp );
    prms_m(i) = sqrt( psqsum(i)/TotalNp );
    //eps_m(i)  = sqrt( max( eps2(i), zero ) );
    eps_norm_m(i)  = sqrt( max( eps2(i), zero ) );
    double tmp    = rrms_m(i) * prms_m(i);
    fac(i)  = (tmp == 0) ? zero : 1.0 / tmp;
  }

  rprms_m = rpsum * fac;

  // calculate mean energy
  mean_energy_m = 0.0;
  for(int k=0; k < locNp; k++)
    mean_energy_m += (sqrt(dot(P[k], P[k]) + 1.0) - 1.0) * m_p * 1000.0;  // Unit: MeV
  localMeanEnergy = mean_energy_m / locNp;
  // sum energy of all nodes
  reduce(mean_energy_m,mean_energy_m,OpAddAssign());
  mean_energy_m /= TotalNp;

  double meanLocalBetaGamma = sqrt( pow( 1 + localMeanEnergy/(1000.0*m_p), 2.0 ) -1 );

  double betagamma = meanLocalBetaGamma*locNp;
  // sum the betagamma of all nodes
  reduce(betagamma,betagamma,OpAddAssign());
  betagamma /= TotalNp;

  // obtain the global RMS emmitance, it make no sense for multi-bunch simulation
  eps_m = eps_norm_m / Vector_t(betagamma);
}

void PartBunch::boundp_destroy() {

  Inform gmsgAll("boundp ",INFORM_ALL_NODES);

  Vector_t len;
  const int dimIdx = 3;
  IpplTimings::startTimer(boundpTimer_m);

  size_t *countLost;
  if(weHaveBins()){
    const int tempN = pbin_m->getLastemittedBin();
    countLost = new size_t[tempN];
    for (int ii=0; ii< tempN; ii++) countLost[ii] = 0;
  }

  NDIndex<3> domain = getFieldLayout().getDomain();
  for(int i=0; i<Dim; i++)
    nr_m[i] = domain[i].length();

  get_bounds(rmin_m, rmax_m);
  len = rmax_m - rmin_m;

  calcBeamParameters_cycl();

  double beamRadius = sqrt( dot( rrms_m, rrms_m) );

  if ( sqrt( dot(len,len)) > 10*sqrt( beamRadius) ){
    double distance = 0;
    for(int ii = 0; ii < this->getLocalNum(); ii++){
      distance = sqrt ( dot( ( R[ii] -  rmean_m ),( R[ii] -  rmean_m ) ) );
      if ( distance > 3*beamRadius ){
        // put particle onto deletion list .
        destroy(1,ii);

        //update bin parameter
        if(weHaveBins()) countLost[Bin[ii]] += 1 ;

        gmsgAll<<"PARTICLE DELETION(too far away): ID = "<<ID[ii]<<", distance="<<distance<<", beam rms radius = "<<beamRadius<<endl;
      }
    }
  }


  for(int i=0; i<dimIdx; i++) {
    rmax_m[i] += dh_m*abs(rmax_m[i]-rmin_m[i]);
    rmin_m[i] -= dh_m*abs(rmax_m[i]-rmin_m[i]);
    hr_m[i]    = (rmax_m[i]-rmin_m[i]) / (nr_m[i]-1);
  }

  // rescale mesh
  getMesh().set_meshSpacing(&(hr_m[0]));
  getMesh().set_origin( rmin_m );

  rho_m.initialize(getMesh(),
                   getFieldLayout(),
                   GuardCellSizes<Dim>(1),
                   bc_m);
  eg_m.initialize(getMesh(),
                  getFieldLayout(),
                  GuardCellSizes<Dim>(1),
                  vbc_m);


  if(weHaveBins()){
    pbin_m->updatePartInBin(countLost);
    if (countLost!=NULL) delete[] countLost;
  }

  update();
  IpplTimings::stopTimer(boundpTimer_m);

}

// angle range [0~2PI) degree
double PartBunch::calculateAngle(double x, double y)
{
  double thetaXY;

  if (x<0)                   thetaXY=pi+atan(y/x);
  else if ((x>0) && (y>=0))  thetaXY=atan(y/x);
  else if ((x>0) && (y<0))   thetaXY=2.0*pi+atan(y/x);
  else if ((x==0) && (y> 0)) thetaXY=pi/2.0;
  else if ((x==0) && (y< 0)) thetaXY=3.0/2.0*pi;

  return thetaXY;

}

double PartBunch::calcMeanPhi() {

  const int emittedBins = pbin_m->getLastemittedBin();

  double phi[emittedBins];
  double px[emittedBins];
  double py[emittedBins];

  for (int ii=0; ii < emittedBins; ii++){
    phi[ii] = 0.0;
    px[ii] = 0.0;
    py[ii] = 0.0;
  }

  for (int ii=0; ii < getLocalNum(); ii++){

    px[Bin[ii]] += P[ii](0);
    py[Bin[ii]] += P[ii](1);
  }

  for (int ii=0; ii < emittedBins; ii++) {

    reduce(px[ii], px[ii], OpAddAssign());
    reduce(py[ii], py[ii], OpAddAssign());

    phi[ii] = calculateAngle(px[ii], py[ii]);
    INFOMSG("Bin " << ii  <<"  mean phi = " << phi[ii]*180.0/pi-90.0 <<"[degree]"<<endl );
  }

  double meanPhi = 0.0;
  for (int ii=0; ii < emittedBins; ii++) meanPhi += phi[ii];

  meanPhi /= emittedBins;

  return meanPhi;
}

size_t PartBunch::getNumPartInBin(int BinID)
{
  if (weHaveBins())
    return pbin_m->getBinCont(BinID);
  else
    return this->getTotalNum();
}

// this function reset the BinID for each particles according to its current gamma
// for the time it is designed for cyclotron where energy gain per turn may be changing.
bool PartBunch::resetPartBinID()
{

  size_t partInBin[numBunch_m];

  if( numBunch_m != pbin_m->getLastemittedBin() )
    {
      ERRORMSG("Bunch number does NOT equal to bin number! Not implemented yet!!!" <<endl);
      return false;
    }

  for(int ii=0; ii<numBunch_m; ii++) partInBin[ii] = 0 ;

  // update mass gamma for each bin first
  INFOMSG("Before reset Bin: "<<endl);
  calcGammas_cycl();

  // reset bin index for each particle and
  // calculate total particles number for each bin
  for (unsigned long int n=0; n < this->getLocalNum(); n++)
    {
      double deltgamma[numBunch_m];
      double gamma = sqrt(1.0 + dot(P[n], P[n]));

      int index = 0;
      for(int ii=0; ii<numBunch_m; ii++)
        deltgamma[ii] = abs( bingamma_m[ii] - gamma );

      for(int ii=0; ii<numBunch_m; ii++)
        if( *(deltgamma+index) > *(deltgamma+ii) )
          index = ii;

      Bin[n]=index;
      partInBin[index]++;
    }

  pbin_m->resetPartInBin( partInBin );

  // after reset Particle Bin ID, update mass gamma of each bin again
  INFOMSG("After reset Bin: "<<endl);
  calcGammas_cycl();

  return true;

}

