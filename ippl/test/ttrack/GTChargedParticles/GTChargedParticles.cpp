
#ifdef GTHDF5
#include <hdf5.h>
#include "H5Part.hh"
#endif

template<class pl>
void ChargedParticles<pl>::openFiles(string baseFn, string title, bool restartMode)
{
  unsigned int N    = getTotalNum();
  double       qTot = getQTot();
  
  if (isRoot()) {
      if (restartMode) {
	  statof_m.open((baseFn+string(".stat.sdds")).c_str(),ios::app);
          statof_m.precision(15);
	  statof_m.setf(ios::scientific, ios::floatfield);
      }
      else {
	  statof_m.open((baseFn+string(".stat.sdds")).c_str(),ios::out);
	  statof_m.precision(15);
	  statof_m.setf(ios::scientific, ios::floatfield);
	  writeStatSDDSHeader(title,N,qTot);
      }
  }
}

template<class pl>
void ChargedParticles<pl>::closeFiles()
{
  if (isRoot()) 
    statof_m.close();
}

template<class pl>
void ChargedParticles<pl>::calcBeamParameters() {
  
  using Physics::c;
  
  Vector_t rsqsum,psqsum,rpsum,eps2,fac;

  const double N =  static_cast<double>(this->getTotalNum());
  const double zero = 0.0;
  
  calcMoments();  
  
  for (unsigned int i=0 ; i<3; i++) {
    rmean_m(i) = centroid_m[2*i]/N;
    pmean_m(i) = centroid_m[(2*i)+1]/N;
    rsqsum(i) = moments_m[2*i][2*i]-N*rmean_m(i)*rmean_m(i);
    psqsum(i) = moments_m[(2*i)+1][(2*i+1)]-N*pmean_m(i)*pmean_m(i);
    rpsum(i) =  moments_m[2*i][(2*i)+1]-N*rmean_m(i)*pmean_m(i);
  }
  eps2      = (rsqsum * psqsum - rpsum * rpsum)/(N*N);
  rpsum /= N;
  
  for (unsigned int i=0 ; i<3; i++) {
    rrms_m(i) = sqrt( rsqsum(i)/N );
    prms_m(i) = sqrt( psqsum(i)/N );
    eps_m(i)  = sqrt( max( eps2(i), zero ) );
    double tmp    = rrms_m(i) * prms_m(i);
    fac(i)  = (tmp == 0) ? zero : 1.0 / tmp;
  }  
  rprms_m = rpsum * fac;
  
  // calculate Courant Snider parameters
  for (unsigned int i=0 ; i<3; i++) {
    csbeta_m(i)  = pow(rrms_m(i),2.0)/eps_m(i);
    csgamma_m(i) = pow(prms_m(i),2.0)/eps_m(i);
    csalpha_m(i) = sqrt((pow(prms_m(i),2.0)*csbeta_m(i)/eps_m(i)) - 1.0);
    csalpha_m(i) = (rprms_m(i) < 0.0) ? csalpha_m(i) : -1.0*csalpha_m(i);
  }
}


template<class pl>
void ChargedParticles<pl>::calcMoments() {
  
  double part[2*3];
  double loc_centroid[2*3];
  double loc_moment[2*3][2*3];
  
  for (int i = 0; i < 2*3; ++i) {
    loc_centroid[i] = centroid_m[i] = 0.0;
    for (int j = 0; j <= i; ++j) {
      loc_moment[j][i] = loc_moment[i][j] = 0.0;
      moments_m[i][j]  = moments_m[j][i]  = 0.0; 
    }
  }
  
  for (unsigned long k=0; k< this->getLocalNum(); ++k){
    part[1] = P[k](0);
    part[3] = P[k](1);
    part[5] = P[k](2) * getGamma(); // LORENTZ TRANSFORMED!!!!!
   
    part[0] = R[k](0);
    part[2] = R[k](1);
    part[4] = R[k](2) * getGamma();// LORENTZ TRANSFORMED!!!!!
    
    for (int i = 0; i < 2*3; ++i) {
      loc_centroid[i] += part[i];
      for (int j = 0; j <= i; ++j) {
	loc_moment[i][j] += part[i] * part[j];
	loc_moment[j][i] = loc_moment[i][j];
      }
    }
  }

  reduce(&(loc_moment[0][0]), &(loc_moment[0][0]) + 2*3*2*3,
         &(moments_m[0][0]), OpAddAssign());
  
  reduce(&(loc_centroid[0]), &(loc_centroid[0]) + 2*3,
         &(centroid_m[0]), OpAddAssign());
}

template<class pl>
void ChargedParticles<pl>::writeStatSDDSHeader(string title, unsigned int N, double qTot)
{
  statof_m << "SDDS1" << endl;
  statof_m << "&description text=\"Statistics data " << title << "\" " << endl;
  statof_m << ", contents=\"stat parameters\" &end" << endl;
  
  statof_m << "&parameter name=Q,  type=double, units=C, ";
  statof_m << "   description=\"Total Bunch Charge\" &end" << endl;

  statof_m << "&parameter name=I0, type=double, units=A, ";
  statof_m << "   description=\"Initial Peak Current\" &end" << endl;

  statof_m << "&parameter name=N, type=long, ";
  statof_m << "   description=\"Number of Particles\" &end" << endl;

  statof_m << "&parameter name=N, type=double, ";
  statof_m << "   description=\"Interaction radius\" &end" << endl;
  

  statof_m << "&column name=T, type=double, units=s, ";
  statof_m << "   description=\"Actual time of the particle\" &end" << endl;

  statof_m << "&column name=gamma, type=double, units=1 , ";
  statof_m << "   description=\"Particle Energy\" &end" << endl;

  statof_m << "&column name=empty, type=double, units=1 , ";
  statof_m << "   description=\" \" &end" << endl;
  
  statof_m << "&column name=rmsx, type=double, units=m , ";
  statof_m << "   description=\"Rms Beamsize in x  \" &end" << endl;
  statof_m << "&column name=rmsy, type=double, units=m , ";
  statof_m << "   description=\"Rms Beamsize in y  \" &end" << endl;
  statof_m << "&column name=rmst, type=double, units=m , ";
  statof_m << "   description=\"Rms Beamsize in t  \" &end" << endl;
  
  statof_m << "&column name=rmspx, type=double, units=m/s , ";
  statof_m << "   description=\"Rms Momenta in x  \" &end" << endl;
  statof_m << "&column name=rmspy, type=double, units=m/s , ";
  statof_m << "   description=\"Rms Momenta in y  \" &end" << endl;
  statof_m << "&column name=rmspt, type=double, units=m/s , ";
  statof_m << "   description=\"Rms Momenta in t  \" &end" << endl;
  
  statof_m << "&column name=ex, type=double, units=m , ";
  statof_m << "   description=\"Emittance x  \" &end" << endl;
  statof_m << "&column name=ey, type=double, units=m , ";
  statof_m << "   description=\"Emittance y  \" &end" << endl;
  statof_m << "&column name=et, type=double, units=m , ";
  statof_m << "   description=\"Emittance t  \" &end" << endl;
  
  statof_m << "&column name=meanx, type=double, units=m , ";
  statof_m << "   description=\"Mean Beamsize in x  \" &end" << endl;
  statof_m << "&column name=meany, type=double, units=m , ";
  statof_m << "   description=\"Mean Beamsize in y  \" &end" << endl;
  statof_m << "&column name=meant, type=double, units=m , ";
  statof_m << "   description=\"Mean Beamsize in t  \" &end" << endl;
  
  statof_m << "&column name=maxx, type=double, units=m , ";
  statof_m << "   description=\"Max Beamsize in x  \" &end" << endl;
  statof_m << "&column name=maxy, type=double, units=m , ";
  statof_m << "   description=\"Max Beamsize in y  \" &end" << endl;
  statof_m << "&column name=maxt, type=double, units=m , ";
  statof_m << "   description=\"Max Beamsize in t  \" &end" << endl;
  
  /*  statof_m << "&column name=betax, type=double, units=m , ";
  statof_m << "   description=\"Beta function in x  \" &end" << endl;
  statof_m << "&column name=betay, type=double, units=m , ";
  statof_m << "   description=\"Beta function in y  \" &end" << endl;
  
  statof_m << "&column name=alphax, type=double, units=1 , ";
  statof_m << "   description=\"Alpha function in x  \" &end" << endl;
  statof_m << "&column name=alphay, type=double, units=1 , ";
  statof_m << "   description=\"Alpha function in y  \" &end" << endl;

  statof_m << "&column name=nIncidentElec, type=long, units=1 , ";
  statof_m << "   description=\"Number of incident electrons per timestep  \" &end" << endl;

  statof_m << "&column name=maxESecElec, type=double, units=eV , ";
  statof_m << "   description=\"Maximal incident electron energy per timestep \" &end" << endl;

  statof_m << "&column name=sumESecElec, type=double, units=eV , ";
  statof_m << "   description=\"Summed incident electron energy per timestep \" &end" << endl;

  statof_m << "&column name=nSecElec, type=long, units=1 , ";
  statof_m << "   description=\"Number of secondary electrons per timestep \" &end" << endl;
  */
  statof_m << "&column name=nColl, type=long, units=1 , ";
  statof_m << "   description=\"Particles collided via Touschek scattering \" &end" << endl;

  statof_m << "&column name=nLoss, type=long, units=1 , ";
  statof_m << "   description=\"Particles lost via Touschek scattering \" &end" << endl;

  statof_m << "&column name=nLoss2, type=long, units=1 , ";
  statof_m << "   description=\"Particles lost via Touschek scattering on x and y axis \" &end" << endl;

  statof_m << "&data mode=ascii &end" << endl;
  
  statof_m << qTot  << endl;
  statof_m << 0.0  << endl;
  statof_m << N  << endl;
  statof_m << interrad_m << endl;
}


template<class pl>
void ChargedParticles<pl>::writeStatistics()
{
  const int pwi = 8;
#ifdef TTProf
      IpplTimings::startTimer(statWritingTimer); 
#endif       
  bounds(R,rmin_m,rmax_m);
  calcBeamParameters();

  if(isRoot()) {
    statof_m << getTime()  << setw(pwi) << "\t";      // 1
    statof_m << getGamma() << setw(pwi)<< "\t";       // 2
    statof_m << 0.0 << setw(pwi)<< "\t";       // 3
    
    for (int i = 0 ; i < 3 ; ++i)
      statof_m << get_rrms()(i) << setw(pwi) << "\t";  // 4 5 6
    
    for (int i = 0 ; i < 3 ; ++i)
      statof_m << get_prms()(i) << setw(pwi) << "\t";  // 7 8 9
    
    for (int i = 0 ; i < 3 ; ++i)
      statof_m << get_emit()(i)   << setw(pwi)<< "\t"; // 10,11,12
    
    for (int i = 0; i < 3; ++i)
      statof_m << get_rmean()(i)  << setw(pwi)<< "\t"; // 13,14,15
    
    for (unsigned int i = 0 ; i < 3 ; ++i)                   // 16
      statof_m << get_maxExtend()(i) << setw(pwi) << "\t";
    
    // write out Courant Snyder parameter 
    /*  statof_m << get_csbeta()(0) << setw(pwi) << "\t";           // 19
    statof_m << get_csbeta()(1) << setw(pwi) << "\t";           // 20
    
    statof_m << get_csalpha()(0) << setw(pwi) << "\t";           // 21
    statof_m << get_csalpha()(1) << setw(pwi) << "\t";           // 22
    
    statof_m << neAtWall_m << setw(pwi) << "\t";                 // 24
    statof_m << maxEsecElec_m << setw(pwi) << "\t";             // 25
    statof_m << eSecElec_m << setw(pwi) << "\t";          // 26
    statof_m << nsTot_m << setw(pwi) << "\t";              // 27
    */
    statof_m << get_collision_num() << setw(pwi) << "\t";  // 28

    statof_m << get_lost_num() << setw(pwi) << "\t"; // 29

    statof_m << get_lost2_num() << setw(pwi) << "\t";
    
    for(unsigned j = 0; j < 3;j++)
      statof_m << get_csfault(j) << setw(pwi);


    statof_m << endl;
  }
#ifdef TTProf
      IpplTimings::stopTimer(statWritingTimer); 
#endif       
}

template<class pl>
void ChargedParticles<pl>::writePartSDDSHeader(string title, unsigned int N, double qTot)
{
  partof_m << "SDDS1" << endl;
  partof_m << "&description text=\"Phase space data " << title << "\" " << endl;
  partof_m << ", contents=\"part parameters\" &end" << endl;
  
  partof_m << "&parameter name=Q,  type=double, units=C, ";
  partof_m << "   description=\"Total Bunch Charge\" &end" << endl;

  partof_m << "&parameter name=I0, type=double, units=A, ";
  partof_m << "   description=\"Initial Peak Current\" &end" << endl;

  partof_m << "&parameter name=N, type=long, ";
  partof_m << "   description=\"Number of Particles\" &end" << endl;

  partof_m << "&column name=x, type=double, units=m, ";
  partof_m << "   description=\" \" &end" << endl;

  partof_m << "&column name=px, type=double, units=m/s , ";
  partof_m << "   description=\" \" &end" << endl;

  partof_m << "&column name=y, type=double, units=m, ";
  partof_m << "   description=\" \" &end" << endl;

  partof_m << "&column name=py, type=double, units=m/s , ";
  partof_m << "   description=\" \" &end" << endl;

  partof_m << "&column name=z, type=double, units=m, ";
  partof_m << "   description=\" \" &end" << endl;

  partof_m << "&column name=pz, type=double, units=m/s , ";
  partof_m << "   description=\" \" &end" << endl;

  partof_m << "&column name=id, type=double, units=1 , ";
  partof_m << "   description=\" unique particle id \" &end" << endl;

  partof_m << "&column name=rho, type=double, units=1 , ";
  partof_m << "   description=\" local density \" &end" << endl;

  partof_m << "&data mode=ascii &end" << endl;
  
  partof_m << qTot  << endl;
  partof_m << 0.0  << endl;
  partof_m << N  << endl;
}


template<class pl>
void ChargedParticles<pl>::writePhaseSpaceSDDS(string baseName) {
  
  int pwi           = 10;

  unsigned int N    = getTotalNum();  
  unsigned int nloc    = getLocalNum();

  double       qTot = getQTot();

  string title("Phase SPACE");
    
  int tag1 = Ippl::Comm->next_tag(IPPL_APP_TAG4, IPPL_APP_CYCLE);
  int tag2 = Ippl::Comm->next_tag(IPPL_APP_TAG5, IPPL_APP_CYCLE);
   
  string Fn;
  char numbuf[6];
  sprintf(numbuf, "%05d", partofCall_m);  
  Fn = baseName + "-" + string(numbuf) + ".part.sdds";
  partofCall_m++;

  if (isRoot()) {
    partof_m.open(Fn.c_str(),ios::out);
    partof_m.precision(15);
    partof_m.setf(ios::scientific, ios::floatfield);
    
    writePartSDDSHeader(title,N,qTot);
    
    for (unsigned i=0; i < nloc; i++)
      partof_m << R[i](0) << setw(pwi) << " \t" 
	       << P[i](0) << setw(pwi) << " \t" 
	       << R[i](1) << setw(pwi) << " \t" 
	       << P[i](1) << setw(pwi) << " \t" 
	       << R[i](2) << setw(pwi) << " \t" 
	       << P[i](2) << setw(pwi) << " \t" 
	       << ID[i]   << setw(pwi) << " \t" 
	       << getLocalDensityAnalytic(R[i],P[i]) << endl;

    partof_m.close();

    for (int node = 1; node < Ippl::getNodes(); node++) {
      Message* smsg = new Message();
      smsg->put(node);
      bool res = Ippl::Comm->send(smsg, node, tag1);	
      if (! res) 
	ERRORMSG("Ippl::Comm->send(smsg, 0, tag1) failed " << endl;);  
      // now wait
      Message* rmsg =  Ippl::Comm->receive_block(node, tag2);
      delete(rmsg);
    } 
  }
  else {    
    int node;
    int recNode = 0;
    Message* rmsg =  Ippl::Comm->receive_block(recNode, tag1);
    rmsg->get(&node);
    if (node == Ippl::myNode()) {
      partof_m.open(Fn.c_str(),ios::app);
      partof_m.precision(15);
      partof_m.setf(ios::scientific, ios::floatfield);
      for (unsigned i=0; i < nloc; i++)
	partof_m << R[i](0) << setw(pwi) << " \t" 
		 << P[i](0) << setw(pwi) << " \t" 
		 << R[i](1) << setw(pwi) << " \t" 
		 << P[i](1) << setw(pwi) << " \t" 
		 << R[i](2) << setw(pwi) << " \t" 
		 << P[i](2) << setw(pwi) << " \t" 
		 << ID[i] << " \t" 
		 << getLocalDensityAnalytic(R[i],P[i]) << endl;
      partof_m.close();
    }
    else
      ERRORMSG("Want save restart node<>root, got wrong node"<<endl;);
    Message* smsg = new Message();
    bool res = Ippl::Comm->send(smsg, 0, tag2);	
    if (! res) 
      ERRORMSG("Ippl::Comm->send(smsg, 0, tag2) failed " << endl;);
  }
}

/**
 * readRestartInfo writes out the restartinfo to the file given by basename.
 * also the parameter in which turn the machine is needs to be passed.
 * @param baseName file the restart information should be saved to.
 * @turn turn which is being computed
 */
template<class pl>
int ChargedParticles<pl>::readRestartInfo(string Fn) {

  Inform msg("Read restart ");

  ifstream rf;
  unsigned num;
  unsigned node;
  double x,y,t,px,py,pt;
  double spos;
  int lost = 0;
  int lost2 = 0;
  int turn;
  unsigned count = 0;


#ifdef GTHDF5

  Fn = Fn + string(".h5");
  msg << "Read restart file  " << Fn << endl;

  H5PartFile *file_m;
  
#ifdef PARALLEL_IO
  file_m=H5PartOpenFileParallel(Fn.c_str(),H5PART_READ,MPI_COMM_WORLD);
#else
  file_m=H5PartOpenFile(Fn.c_str(),H5PART_READ);
#endif
  if(!file_m) {
      INFOMSG("File open failed:  exiting!" << endl);
      exit(0);
  }   

  H5PartSetStep(file_m,1);      

  int nt,nds;
  char name[64];
  nt =H5PartGetNumSteps(file_m);     /* get number of steps in file */
  nds=H5PartGetNumDatasets(file_m);  /* get number of datasets in timestep 0 */ 
  unsigned  N=(unsigned)H5PartGetNumParticles(file_m);

  H5PartReadAttrib(file_m,"Step" ,&turn);
  H5PartReadAttrib(file_m,"Spos" ,&spos);
  H5PartReadAttrib(file_m,"lost2",&lost2);  
  H5PartReadAttrib(file_m,"tl" ,&lost);

  long long starti = 0;
  long long endi   = 0;
  long long  *globN = (long long*) malloc(Ippl::getNodes()*sizeof(long long));
  H5PartReadAttrib(file_m,"nloc",globN);
  
  if (Ippl::myNode() != 0) {
      for(int i=0; i< Ippl::myNode(); i++) 
	  starti += globN[i];
  }
  unsigned Nloc = globN[Ippl::myNode()];

  endi = starti + Nloc;

  H5PartSetView(file_m,starti,endi);
  
  void *varray = malloc(Nloc*8);
  double *farray = (double*)varray;
  long long int *larray = (long long int *)varray;

  create(Nloc);

  H5PartReadDataFloat64(file_m,"x",farray);
  for(unsigned i=0; i < Nloc; i++) 
      R[i] += Vector_t(farray[i],0.0,0.0); 
  
  H5PartReadDataFloat64(file_m,"y",farray);
  for(unsigned i=0; i < Nloc; i++) 
      R[i] += Vector_t(0.0,farray[i],0.0); 

  H5PartReadDataFloat64(file_m,"z",farray);
  for(unsigned i=0; i < Nloc; i++) 
      R[i] += Vector_t(0.0,0.0,farray[i]); 
  
  H5PartReadDataFloat64(file_m,"px",farray);
  for(unsigned i=0; i < Nloc; i++) 
      P[i](0) = farray[i];
  
  H5PartReadDataFloat64(file_m,"py",farray);
  for(unsigned i=0; i < Nloc; i++) 
      P[i](1) = farray[i];  

  H5PartReadDataFloat64(file_m,"pz",farray);
  for(unsigned i=0; i < Nloc; i++) 
      P[i](2) = farray[i];

  H5PartCloseFile(file_m);
  
  if(farray)
      free(farray);
#else
  char tmp[120];
  rf.open(Fn.c_str(),ios::in);

  rf.getline(tmp,120,'\n');
  rf.getline(tmp,120,'\n');
  rf.getline(tmp,120,'\n');
  
  rf >> turn;
  rf >> num;
  rf >> spos;
  rf >> lost;
  rf >> lost2;
  
  msg << "Restart mode using restartfile " << Fn
      << " total number of Particles to read: " << num << endl;
    
  for(unsigned i=0; i < num; i++) {
    rf >> x >> px >> y >> py >> t >> pt >> node;
    if (node == Ippl::myNode()) {
      create(1);
      R[count] = Vector_t(x,y,t); 
      P[count] = Vector_t(px,py,pt); 
      count++;
      }
    }
    rf.close();
#endif
    set_spos(spos);
    add_lost(lost);
    add_lost2(lost2);
    boundp();
    return turn;	
}

/**
 * writeRestartInfo writes out the restartinfo to the file given by basename.
 * also the parameter in which turn the machine is needs to be passed.
 * @param baseName file the restart information should be saved to.
 * @turn turn which is being computed
 */
template<class pl>
void ChargedParticles<pl>::writeRestartInfo(string Fn, unsigned turn) {
 
  int pwi       = 10;
  unsigned nloc = getLocalNum();
  unsigned myN  = Ippl::myNode();
  
#ifdef GTHDF5

  Inform m("H5PartWrite ");

  H5PartFile *file_m;
  
#ifdef PARALLEL_IO
  file_m=H5PartOpenFileParallel(Fn.c_str(),H5PART_WRITE,MPI_COMM_WORLD);
#else
  file_m=H5PartOpenFile(Fn.c_str(),H5PART_WRITE);
#endif
  if(!file_m) {
      INFOMSG("File open failed:  exiting!" << endl);
      exit(0);
  }   

  void *varray = malloc(nloc*sizeof(double));
  double *farray = (double*)varray;
  
  /* ------------------------------------------------------------------------ 
     Get the particle decomposition from all the nodes
  */
  long long *locN = (long long *) malloc(Ippl::getNodes()*sizeof(long long));
  long long  *globN = (long long*) malloc(Ippl::getNodes()*sizeof(long long));
    
  for(int i=0; i<Ippl::getNodes(); i++) {
      globN[i] = locN[i]=0;
  }
  locN[Ippl::myNode()] = nloc;
  reduce(locN, locN + Ippl::getNodes(), globN, OpAddAssign());
    
  /* ------------------------------------------------------------------------ */
    
  unsigned nTot           = getTotalNum();

  H5PartSetStep(file_m,turn);  
  H5PartSetNumParticles(file_m,nloc); 
  
  /* write scalar data i.e the header */
  long long step = turn;
  H5PartWriteStepAttrib(file_m,"tl", H5T_NATIVE_INT64,&lost_num,1);
  H5PartWriteStepAttrib(file_m,"Step", H5T_NATIVE_INT64,&step,1);
  H5PartWriteStepAttrib(file_m,"Ntot", H5T_NATIVE_INT64,&nTot,1);

  H5PartWriteStepAttrib(file_m,"lost2", H5T_NATIVE_INT64,&lost2_num,1);
  H5PartWriteAttrib(file_m,"Spos",     H5T_NATIVE_DOUBLE,&spos_m,1);
  
  H5PartWriteAttrib(file_m,"nloc",H5T_NATIVE_INT64, globN, Ippl::getNodes());

  m << "Write restart file turn " << turn << " s= " << spos_m << 
      " nTot= " << nTot << " lost= " << lost_num << " lost2= " << lost2_num << endl;

  for (long long i=0; i<nloc;i++)
      farray[i] =  R[i](0);
  H5PartWriteDataFloat64(file_m,"x",farray); 
  
  for (long long i=0; i<nloc;i++)
      farray[i] =  R[i](1);
  H5PartWriteDataFloat64(file_m,"y",farray); 
  
  for (long long i=0; i<nloc;i++)
      farray[i] =  R[i](2);
  H5PartWriteDataFloat64(file_m,"z",farray); 
    
  for (long long i=0; i<nloc;i++)
      farray[i] =  P[i](0);
  H5PartWriteDataFloat64(file_m,"px",farray); 

  for (long long i=0; i<nloc;i++)
      farray[i] =  P[i](1);
  H5PartWriteDataFloat64(file_m,"py",farray); 
  
  for (long long i=0; i<nloc;i++)
      farray[i] =  P[i](2);
  H5PartWriteDataFloat64(file_m,"pz",farray); 
  
  if(varray)  
      free(varray);

  H5PartCloseFile(file_m);
#else

  int tag1 = Ippl::Comm->next_tag(IPPL_APP_TAG4, IPPL_APP_CYCLE);
  int tag2 = Ippl::Comm->next_tag(IPPL_APP_TAG5, IPPL_APP_CYCLE);

  if (isRoot()) {
    partof_m.open(Fn.c_str(),ios::out);
    partof_m.precision(15);
    partof_m.setf(ios::scientific, ios::floatfield);
    
    partof_m << "& Restart Information" << endl;
    partof_m << "& turn, number, spos, lost, lost2" << endl;
    partof_m << "& particles: x, px, y, py, z, pz node" << endl;
    partof_m << turn << endl;
    partof_m << getTotalNum() << endl;
    partof_m << spos_m << endl;
    partof_m << lost_num << endl;
    partof_m << lost2_num << endl;
    //    gräfli isch doch schön im engadin!
    for (unsigned i=0; i < nloc; i++)
      partof_m << R[i](0) << setw(pwi) << " \t" 
	       << P[i](0) << setw(pwi) << " \t" 
	       << R[i](1) << setw(pwi) << " \t" 
	       << P[i](1) << setw(pwi) << " \t" 
	       << R[i](2) << setw(pwi) << " \t" 
	       << P[i](2) << setw(pwi) << " \t" << myN << endl;
    partof_m.close();

    for (int node = 1; node < Ippl::getNodes(); node++) {
      Message* smsg = new Message();
      smsg->put(node);
      bool res = Ippl::Comm->send(smsg, node, tag1);	
      if (! res) 
	ERRORMSG("Ippl::Comm->send(smsg, 0, tag1) failed " << endl;);  
      // now wait
      Message* rmsg =  Ippl::Comm->receive_block(node, tag2);
      delete(rmsg);
    } 
  }
  else {    
    int node;
    int recNode = 0;
    Message* rmsg =  Ippl::Comm->receive_block(recNode, tag1);
    rmsg->get(&node);
    if (node == Ippl::myNode()) {
      partof_m.open(Fn.c_str(),ios::app);
      partof_m.precision(15);
      partof_m.setf(ios::scientific, ios::floatfield);
      for (unsigned i=0; i < nloc; i++)
	partof_m << R[i](0) << setw(pwi) << " \t" 
		 << P[i](0) << setw(pwi) << " \t" 
		 << R[i](1) << setw(pwi) << " \t" 
		 << P[i](1) << setw(pwi) << " \t" 
		 << R[i](2) << setw(pwi) << " \t" 
		 << P[i](2) << setw(pwi) << " \t" << myN << endl;
      partof_m.close();
    }
    else
      ERRORMSG("Want save restart node<>root, got wrong node"<<endl;);
    Message* smsg = new Message();
    bool res = Ippl::Comm->send(smsg, 0, tag2);	
    if (! res) 
      ERRORMSG("Ippl::Comm->send(smsg, 0, tag2) failed " << endl;);
  }
#endif
}

template<class pl>
void ChargedParticles<pl>::send_lose_particle() {
  
  if(!isRoot()) {

    int tag = Ippl::Comm->next_tag(IPPL_APP_TAG4, IPPL_APP_CYCLE);
    unsigned lost_num = 0;
    Message *loss_msg = new Message();
    
    for(unsigned i = 0; i < getLocalNum(); i++)
      if(lost[i]) lost_num++;
    
    loss_msg->put(lost_num);
  
    for(unsigned i = 0; i < getLocalNum(); i++) {
      if(lost[i]) {
	for(unsigned j = 0; j < 3; j++) {
	  loss_msg->put(R[i](j));      
	  loss_msg->put(P[i](j));
	}
	loss_msg->put(get_spos());
	lost[i] = 0;
	destroy(1,i);
      }
    }
    
    bool res = Ippl::Comm->send(loss_msg, 0, tag);
    
    if(!res) ERRORMSG("Ippl::Comm->send(loss_msg, 0, tag) failed " << endl;); 
  }
}


template<class pl>
void ChargedParticles<pl>::write_lose_particle(string file) {

  int tag = Ippl::Comm->next_tag(IPPL_APP_TAG4, IPPL_APP_CYCLE);
  int pwi = 10;

  if(isRoot()) {
    partof_m.open(file.c_str(),ios::app);
    partof_m.precision(15);
    partof_m.setf(ios::scientific, ios::floatfield);    

    for(unsigned i = 0; i < getLocalNum(); i++) {
      if(lost[i]) {
	for(unsigned j = 0; j < 3; j++) {
	  partof_m  << R[i](j) << " \t"
		    << P[i](j) << " \t";      
	}
	partof_m << get_spos()  << endl;
	lost[i] = 0;
	destroy(1,i);
      }
    }
    
    int notReceived =  Ippl::getNodes()-1;
    while (notReceived > 0) {   
      int node = COMM_ANY_NODE;
      double x,y,z,px,py,pz,time;
      unsigned dataBlocks;
      Message* rmsg =  Ippl::Comm->receive_block(node, tag);
      if (rmsg == 0)
	ERRORMSG("Could not receive from client nodes in main." << endl);
      notReceived--;
      rmsg->get(&dataBlocks);
      for (unsigned i=0; i < dataBlocks; i+=7) {
	rmsg->get(&x);
	rmsg->get(&px);
	rmsg->get(&y);
	rmsg->get(&py);
	rmsg->get(&z);
	rmsg->get(&pz);
	rmsg->get(&time);
	partof_m  << x  << setw(pwi) << " \t"  
		  << px << setw(pwi) << " \t"  
		  << y  << setw(pwi) << " \t"             
		  << py << setw(pwi) << " \t"  
		  << z <<  setw(pwi) << " \t"  
		  << pz << setw(pwi) << " \t"             
		  << time  << endl;
      }
      delete(rmsg);
    } 
    partof_m.close();
  }
}

template<class pl>
void ChargedParticles<pl>::init_lose_particle(string file) {
  partof_m.open(file.c_str(),ios::out);
  partof_m.close();
}
