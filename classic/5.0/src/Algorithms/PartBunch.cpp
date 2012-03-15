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
#include "Physics/Physics.h"

using Physics::c;

// Class PartBunch
// ------------------------------------------------------------------------

PartBunch::PartBunch()
{
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

  dh_m = 0;
  fs_m = 0;
  dt = 0.0;
  pbin_m = 0;
  trackStep_m = 0;
  numBunch_m = 0;
}


PartBunch::PartBunch(const PartBunch &rhs)
{}


PartBunch::PartBunch(const std::vector<Particle> &rhs)
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

}



void PartBunch::calcLineDensity() {
  /*


  e_dim_tag decomp2[3]; 
  for (int d=0; d < 3; ++d) {
  decomp[d] = layout_m->getRequestedDistribution(d);
  }




   *gmsg << "Q= " << sum(this->Q) << end;

    this->Q.scatter(this->rho_m, this->R, IntrplCIC_t());
    calcZDensity(rho_m,decomp,string("DensityZ.dat");

    rho_m = 0;
    this->R(0).scatter(this->rho_m, this->R, IntrplCIC_t());
    calcZDensity(rho_m,decomp,string("DensityXs.dat");
    
  */



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

double PartBunch::getMaxdEBins() {

  const int emittedBins = pbin_m->getLastemittedBin();

  double maxdE = DBL_MIN;
  if (emittedBins >= 2){
    for (int i=1; i<emittedBins+1; i++) {	
      double de = abs(Physics::m_e*1.0E6*(bingamma_m[i-1]-bingamma_m[i]));
      if (de > maxdE)
	maxdE=de;
    }
    return maxdE;
  }
  else 
    return DBL_MAX;
}




void PartBunch::computeSelfFields(int bin) 
{ 
  IpplTimings::startTimer(selfFieldTimer_m);

  rho_m = 0.0;
  eg_m = Vector_t(0.0);

  if (fs_m->hasValidSolver()) {
    this->Q *= this->dt;
    this->Q.scatter(this->rho_m, this->R, IntrplCIC_t());
    this->Q /= this->dt;
    this->rho_m /= getdT();

    // calculating mesh-scale factor, gamma now from a specific slice
    double scaleFactor = Physics::c * getdT();
    double gammaz = getBinGamma(bin);
    
    double tmp2 = 1/hr_m[0] * 1/hr_m[1] * 1/hr_m[2] / (scaleFactor*scaleFactor*scaleFactor) / gammaz;

    rho_m *= tmp2;
    
    Vector_t hr_scaled = hr_m * Vector_t(scaleFactor);
    hr_scaled[2] *= gammaz;

    Vector_t rmax, rmin;
    get_bounds(rmin,rmax);
    double zshift = -(rmax(2) + rmin(2))*gammaz*scaleFactor;

    IpplTimings::startTimer(compPotenTimer_m);
    // charge density is in rho_m
    fs_m->solver_m->computePotential(rho_m, hr_scaled, zshift);
    IpplTimings::stopTimer(compPotenTimer_m);

    //scale mesh back
    rho_m *= hr_scaled[0]*hr_scaled[1]*hr_scaled[2];

    // the scalar potential is given back in rho_m
    // and must be converted to the right units
    rho_m *= getCouplingConstant();
 
    // IPPL Grad divides by hr_m
    eg_m = -Grad(rho_m,eg_m);
  
    eg_m *= Vector_t(gammaz/(scaleFactor), gammaz/(scaleFactor), 1.0/(scaleFactor*gammaz));

    // interpolate electric field at particle positions.  We reuse the
    // cached information about where the particles are relative to the
    // field, since the particles have not moved since this the most recent
    // scatter operation.

    Eftmp.gather(eg_m, this->R,  IntrplCIC_t());

    /** Magnetic field in x and y direction induced by the eletric field
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
    /*double sum_boundary = 0.0;
      ofstream fstr2;
      fstr2.precision(9);
      fstr2.open("rho",ios::out);
      for(int k=0; k<32; k++) {
      for(int i=0; i<32; i++) {
      for(int j=0; j<32; j++) {
      fstr2 << k+1 << " " << i+1 << " " << j+1 << " " << rho_m[k][i][j] << endl;
      }
      }
      }
      fstr2.close();
    */
    
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

  /// set initial charge dentsity to zero.
  rho_m = 0.0;

  /// set initial E field to zero
  eg_m = Vector_t(0.0);
  
  if (fs_m->hasValidSolver()) {
    
    /// scatter particles charge onto grid.
    this->Q.scatter(this->rho_m, this->R, IntrplCIC_t());

    double tmp2 = 1/hr_m[0] * 1/hr_m[1] * 1/hr_m[2] / gamma;

    /// from charge to charge density.
    rho_m *= tmp2;

    /// Lorentz transformation
    /// In particle rest frame, the longitudinal length enlarged
    Vector_t hr_scaled = hr_m ;
    hr_scaled[2] *= gamma;

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
    eg_m *= Vector_t(gamma, gamma, 1.0);

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

    /** Magnetic field in x and y direction induced by the eletric field
     *
     *  \f[ B_x = \gamma(B_x^{'} - \frac{beta}{c}E_y^{'}) = -\gamma \frac{beta}{c}E_y^{'} = -\frac{beta}{c}E_y \f]
     *  \f[ B_y = \gamma(B_y^{'} - \frac{beta}{c}E_x^{'}) = +\gamma \frac{beta}{c}E_x^{'} = +\frac{beta}{c}E_x \f]
     *  \f[ B_z = B_z^{'} = 0 \f]
     *
     */
    double betaC = sqrt(gamma*gamma - 1.0)/gamma/Physics::c;

    // realy compact!
    Bf(0) = Bf(0) - betaC * Ef(1);
    Bf(1) = Bf(1) + betaC * Ef(0);
    
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

void PartBunch::calcBeamParameters()
{
  using Physics::c;
  using Physics::m_e;
  
  Vector_t eps2,fac,rsqsum,psqsum,rpsum;

  IpplTimings::startTimer(statParamTimer_m);

  const size_t locNp = this->getLocalNum();
  
  const double zero = 0.0;
  const double N =  static_cast<double>(this->getTotalNum());

  calcMoments();  

  for (unsigned int i=0 ; i<Dim; i++) {
    rmean_m(i) = centroid_m[2*i]/N;
    pmean_m(i) = centroid_m[(2*i)+1]/N;
    rsqsum(i) = moments_m(2*i,2*i)-N*rmean_m(i)*rmean_m(i);
    psqsum(i) = moments_m((2*i)+1,(2*i)+1)-N*pmean_m(i)*pmean_m(i);
    if(psqsum(i) < 0)
      psqsum(i) = 0;
    rpsum(i) =  moments_m((2*i),(2*i)+1)-N*rmean_m(i)*pmean_m(i);
  }
  eps2      = (rsqsum * psqsum - rpsum * rpsum)/(N*N);
  rpsum /= N;

  for (unsigned int i=0 ; i<Dim; i++) {
    rrms_m(i) = sqrt( rsqsum(i)/N );
    prms_m(i) = sqrt( psqsum(i)/N );
    eps_m(i)  = sqrt( max( eps2(i), zero ) );
    double tmp    = rrms_m(i) * prms_m(i);
    fac(i)  = (tmp == 0) ? zero : 1.0 / tmp;
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
  
  //normalize emittance
  double betagamma = 0.0;
  for(size_t i=0; i<locNp;i++)
    betagamma += sqrt(1.0 + dot(P[i],P[i]));
  reduce(betagamma,betagamma,OpAddAssign());
  betagamma /= N;
  betagamma *= sqrt(1.0 - (1/betagamma)*(1/betagamma));
  eps_norm_m = eps_m / Vector_t(betagamma);

  IpplTimings::stopTimer(statParamTimer_m);
}


void PartBunch::setSolver(FieldSolver &fs) {
  fs_m = &fs;
  fs_m->initSolver(*this);
  initialize(&fs_m->getParticleLayout());
  /*
  if(!fs_m->hasValidSolver()) 
    fs_m = 0;
  else {
    //: update/initalize the particle container with Layout information

  }
  */
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
 * Here we emit all particles from bin 'bin' that would step over the cathode (z=0) in the first half step of the integration scheme. 
 This particles are added to the bunch and reset their timestep to the remaining distance they would cover until after the second half step diveded by two. 
 The timestep for these particles will be shorter in this step than the full timestep of the simulation. 
 The timestep of emitted particles is reset to the simluation timestep after the full timestep has been applied.
 \f[
    \Delta t = \frac{P_z}{\beta c} + \frac{\Delta t_{full-timestep}}{2}
 \f]
 The particles that do not step over the cathode update their z position.
 */
size_t PartBunch::emitFirstHalfStep(int bin) {
  /**
     We copy all particles in the given
     bin from pbin_m to the particle 
     container.
  */  
  size_t nloc = this->getLocalNum();
  size_t oldtot = nloc;
    
  if(!pbin_m->getBinEmitted(bin) && bin < pbin_m->getNBins()) {

    //*gmsg << "* ************** E M I T *********************************************************** " << endl;
    //*gmsg << "* BIN= " << bin << " out of " << pbin_m->getNBins() << endl;
    //  *gmsg << "* rmin_i= " << rmin_m << " rmax_i= " << rmax_m << " h_i= " << hr_m << endl; 

    double bingamma = 0.0;
    
    //FIXME: DONT DO EVERY EMIT CALL!
    for (size_t i=0; i<pbin_m->getNp(); i++) {
      vector<double> p;
      if (pbin_m->getPart(i,bin,p)) {
        bingamma += sqrt(1.0 + p[3]*p[3] + p[4]*p[4] + p[5]*p[5]);
      }
    }
    bingamma /= pbin_m->getBinCont(bin); 
    double betac = sqrt(1.0 - (1.0/(bingamma*bingamma)))*Physics::c;

    for (size_t i=0; i<pbin_m->getNp(); i++) {
      vector<double> p;
      double tempz;
      if (pbin_m->getPart(i,bin,p)) {
        if( !pbin_m->isEmitted(i,bin) ) {
//           *gmsg << "p[2]= " << endl;// << p[2] << endl;
          tempz = getdT()/2.0 * betac;
          p[2] += tempz;
          if(p[2] > 0.0) {
            create(1);
            /** calculate the time when the particle is emitted, \f$t_{emit}\f$, and subtract this from \f$dT\f$ to get the proper \f$dt\f$ for the particle. 
		Set then the particle back to the surface of the cathode and push it from there to the right position with its own \f$dt\f$.
                \f[
                   t_{emit} = \frac{dT}{2} - \frac{z}{\beta c}\\
                   dt = dT - t_{emit} = \frac{dT}{2} + \frac{z}{\beta c} \\
                   z' = \frac{dt}{2} \beta c = \frac{dT}{4}\beta c + \frac{z}{2}
                \f]

            */
            this->R[nloc] = Vector_t(p[0],p[1],getdT()/4.0 * betac + p[2]/2.0);
            this->P[nloc] = Vector_t(p[3],p[4],p[5]);
            this->Bin[nloc]=bin;
            this->Q[nloc] = getChargePerParticle();
            this->dt[nloc] = p[2]/betac + getdT()/2.0; //getdT() + p[2] - getdT()/2.0*betac;
	    this->LastSection[nloc]=0;

            pbin_m->setEmitted(i,bin);
            nloc++;
            binemitted_m[bin]++;
            } else {
              //update z position of particle 'i' in bin 'bin'
              pbin_m->updatePartPos(i, bin, p[2]);
//             }
          }
        }
      }
    }
    
    if(binemitted_m[bin] == pbin_m->getBinCont(bin)) {
      pbin_m->setBinEmitted(bin);
      *gmsg << "* Bin " << bin << " has emitted all particles" << endl;
    }
    /* 
      Because only node 0 is doing that we can not
      do an boundp here!
      boundp();
    */

    if ((nloc-oldtot) != 0) {
      pbin_m->setActualemittedBin(bin); 
    }
  }
  return (nloc-oldtot);
}

/**
 * Here we emit all particles from bin 'bin' that would step over the cathode (z=0) in the second half step of the integration scheme. 
   This particles are added to the bunch with a z position of 0 since they dont get the kick from the second timestep. 
   Their timestep is reset to a full timestep plus the additional lost time due to not getting the kick in the second half step. 
   The timestep in the next step will be larger the the full timestep of the simulation. 
   The timestep of emitted particles is reset to the simluation timestep after the full timestep has been applied.
 \f[
    \Delta t = \Delta t_{full-timestep} + \frac{P_z}{\beta c}
 \f]
 The particles that do not step over the cathode update their z position.
 */
size_t PartBunch::emitSecondHalfStep(int bin) {
  /**
     We copy all particles in the given
     bin from pbin_m to the particle 
     container.
  */  
  
  size_t nloc = this->getLocalNum();
  size_t oldtot = nloc;

  if(!pbin_m->getBinEmitted(bin) && bin < pbin_m->getNBins()) {

    //*gmsg << "* ************** E M I T *********************************************************** " << endl;
    //*gmsg << "* BIN= " << bin << " out of " << pbin_m->getNBins() << endl;
    //  *gmsg << "* rmin_i= " << rmin_m << " rmax_i= " << rmax_m << " h_i= " << hr_m << endl; 

    double bingamma = 0.0;

    //FIXME: DONT DO EVERY EMIT CALL!
    for (size_t i=0; i<pbin_m->getNp(); i++) {
      vector<double> p;
      if (pbin_m->getPart(i,bin,p)) {
        bingamma += sqrt(1.0 + p[3]*p[3] + p[4]*p[4] + p[5]*p[5]);
      }
    }
    bingamma /= pbin_m->getBinCont(bin); 
    double betac = sqrt(1.0 - (1.0/(bingamma*bingamma)))*Physics::c;

    for (size_t i=0; i<pbin_m->getNp(); i++) {
      vector<double> p;
      if (pbin_m->getPart(i,bin,p)) {
        if( !pbin_m->isEmitted(i,bin) ) {
          p[2] += getdT()/2.0 * betac;
          if(p[2] > 0.0) {
            create(1);
            /** calculate the time when the particle is emitted, \f$t_{emit}\f$, and subtract this from \f$dT\f$ to get the proper \f$dt\f$ for the particle. 
		Set then the particle back to the surface of the cathode and push it from there to the right position with its own \f$dt\f$.
                \f[
                   t_{emit} = dT - \frac{z}{\beta c}\\
                   dt = dT - t_{emit} = \frac{z}{\beta c} \\
                   z' = \frac{dt}{2} \beta c = \frac{z}{2}
                \f]

            */
            this->R[nloc] = Vector_t(p[0],p[1],p[2]/2.0);
            this->P[nloc] = Vector_t(p[3],p[4],p[5]);
            this->Bin[nloc]=bin;
            this->Q[nloc] = getChargePerParticle();
            this->dt[nloc] = p[2]/betac;
	    this->LastSection[nloc]=0;
            pbin_m->setEmitted(i,bin);
            nloc++;
            binemitted_m[bin]++;
          } else {
            //update z position of particle 'i' in bin 'bin'
            pbin_m->updatePartPos(i, bin, p[2]);
          }
        }
      }
    }
    if(binemitted_m[bin] == pbin_m->getBinCont(bin)){
      pbin_m->setBinEmitted(bin);
      *gmsg << "* Bin " << bin << " has emitted all particles" << endl;
    }
    /* 
       Because only node 0 is doing that we can not
       do an boundp here!
       boundp();
    */
  }
  return (nloc-oldtot);
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
  Inform msg ("bouldp  ");

  /*
    Assume rmin_m < 0.0 
  */
  //  if (getTotalNum() > 0) {
    Vector_t len;
    const int dimIdx = 3;
    IpplTimings::startTimer(boundpTimer_m);
    
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
        if ( distance > 5*beamRadius ){
          destroy(1,ii);
          msg <<"Particle "<<ID[ii]<<"terminated! distance="<<distance<<" beamRadius = "<<beamRadius<<endl;
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
    update();
    IpplTimings::stopTimer(boundpTimer_m);
    //  }
}
