#ifndef SIMDATA
#define SIMDATA
class SimData
{
public:
  
  SimData(int argc, char *argv[])
  {
    Inform msg("SimData ");

    FileStream::setEcho(false);
    FileStream *is;
    
    try {
      is = new FileStream(string(argv[1]));
    } catch (...) {
      is = 0;
      ERRORMSG("No input file \"" << string(argv[1]) << "\" found." << endl);
    }
    skip = is->readToken().getInteger(); 
    jinit = is->readToken().getInteger(); 
    icbc = is->readToken().getInteger(); 
    norm = is->readToken().getInteger(); 
    trans = is->readToken().getInteger(); 
    initp = is->readToken().getInteger(); 
    norder = is->readToken().getInteger(); 
    frm = is->readToken().getInteger(); 
    pbchoice = is->readToken().getInteger(); 
    entest = is->readToken().getInteger(); 
    ng = is->readToken().getInteger(); 
    ng2d = is->readToken().getInteger(); 
    ns = is->readToken().getReal(); 
    iseed = is->readToken().getInteger(); 
    alpha = is->readToken().getReal(); 
    zin = is->readToken().getReal(); 
    zfin = is->readToken().getReal(); 
    zpr = is->readToken().getReal(); 
    np = is->readToken().getInteger(); 
    nsteps = is->readToken().getInteger(); 
    iprint = is->readToken().getInteger(); 
    irun = is->readToken().getInteger(); 
    hubble = is->readToken().getReal(); 
    omegadm = is->readToken().getReal(); 
    deut = is->readToken().getReal(); 
    rL = is->readToken().getReal(); 
    ss8 = is->readToken().getReal(); 
    qq = is->readToken().getReal(); 
    nxint = is->readToken().getInteger(); 
    nxdint = is->readToken().getInteger(); 
  }
  
  void print(ostream &os) const
  {
    os << "Simulation  input data: " << endl;
    os << "H        \t" << hubble << endl;
    os << "SKIP     \t" << skip << endl;     
    os << "JINIT    \t" << jinit << endl;     
    os << "ICBC     \t" << icbc << endl;     
    os << "NORM     \t" << norm << endl;     
    os << "TRANS    \t" << trans << endl;     
    os << "INITP    \t" << initp << endl;     
    os << "NORDER   \t" << norder << endl;     
    os << "FRM      \t" << frm << endl;     
    os << "PBCHOICE \t" << pbchoice << endl;     
    os << "ENTEST   \t" << entest << endl;     
    os << "NG       \t" << ng << endl;     
    os << "NG2D     \t" << ng2d << endl;     
    os << "NS       \t" << ns << endl;     
    os << "ISEED    \t" << iseed << endl;     
    os << "ALPHA    \t" << alpha << endl;     
    os << "ZIN      \t" << zin << endl;     
    os << "ZFIN     \t" << zfin << endl;     
    os << "ZPR      \t" << zpr << endl;     
    os << "NP       \t" << np << endl;     
    os << "NSTEPS   \t" << nsteps << endl;     
    os << "IPRINT   \t" << iprint << endl;     
    os << "IRUN     \t" << irun << endl;     
    os << "H        \t" << hubble << endl;     
    os << "OMEGADM  \t" << omegadm << endl;     
    os << "DEUT     \t" << deut << endl;     
    os << "HPHYS    \t" << rL << endl;     
    os << "SIGMA8   \t" << ss8 << endl;     
    os << "QQ       \t" << qq << endl;     
    os << "NXINT    \t" << nxint << endl;     
    os << "NXDINT   \t" << nxdint << endl;     
  }
  // read in from the input file
  int ng;
  int np;
  int nsteps;
  int norder;
  int iprint;
  int iseed;
  int jinit;
  int norm;
  int trans;
  int nxint;
  int nxdint;
  double ns;
  /// initial redshift
  double zin;
  /// final redshift
  double zfin;

  double zpr;
  /// Hubble constant/100 km/s/Mpc
  double hubble;
  double omegadm;
  double rL;
  double alpha;
  double deut;
  double qq;
  double ss8;
  int initp;
  int skip;
  int frm;
  int pbchoice;
  int irun;
  int icbc;
  double pmass;
  int entest;
  int ng2d;
  int nginit;  

  // derived quantities
  

};


// Output operator.
inline ostream &operator<<(ostream &os, const SimData &data)
{
  data.print(os);
  return os;
}
#endif














