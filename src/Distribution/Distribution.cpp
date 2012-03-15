// ------------------------------------------------------------------------
// $RCSfile: Distribution.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.3.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Distribution
//   The class for the OPAL Distribution command.
//
// ------------------------------------------------------------------------

#include "Distribution/Distribution.h"
#include "BasicActions/Option.h"
#include "AbstractObjects/OpalData.h"
#include <fstream>
#include <string>
#include <cmath>
#include <cfloat>
#include <vector>
#include <iostream>  // Neeeded for stream I/O
#include <fstream>   // Needed for file I/O
#include <iomanip>   // Needed for I/O manipulators


#include "Distribution/Bins.h"

using namespace Expressions;
using namespace Physics;
using namespace Attributes;

// Class Distribution
// ------------------------------------------------------------------------

// The attributes of class Distribution.

namespace {
    enum {
        // DESCRIPTION OF THE DISTRIBUTION:
        DISTRIBUTION,
        FNAME,
        XMULT,
        YMULT,
        TMULT,
        PXMULT,
        PYMULT,
        PTMULT,
        BETAX,
        BETAY,
        ALPHAX,
        ALPHAY,
        MX,
        MY,
        MT,
        DX,
        DDX,
        DY,
        DDY,
        R51,
        R52,
        R61,
        R62,
        PT,
        T,
        SIGMAX,
        SIGMAY,
        SIGMAT,
        TRANSVCUTOFF,
        SIGMAPX,
        SIGMAPY,
        SIGMAPT,
        RISETIME,
        FALLTIME,
        FLATTOPTIME,
        CORRX,
        CORRY,
        CORRT,
        TEMISSION,
        NBIN,
        DEBIN,
        ELASER,
        SIGLASER,
        W,
        FE,
        AG,
        SIZE
    };
}

Distribution::Distribution():
    Definition(SIZE,"DISTRIBUTION","The DISTRIBUTION statement defines data for the 6D particle distr.")
{
    itsAttr[DISTRIBUTION] = makeString("DISTRIBUTION", "Distribution type: GAUSS, BINOMIAL, ROTSYMBINOMIAL, FROMFILE,"
                                       "GUNGAUSS, GUNGAUSS3D, GUNUNIFORM, GUNGAUSSFLATTOP, GUNGAUSSFLATTOPTH, UNIFORMXYZ ", "GAUSS");

    itsAttr[FNAME] = makeString("FNAME", "File for read in 6D particle coordinates");

    itsAttr[XMULT] = makeReal("XMULT", "Multiplier for X",1.0);
    itsAttr[YMULT] = makeReal("YMULT", "Multiplier for Y",1.0);
    itsAttr[TMULT] = makeReal("TMULT", "Multiplier for T",1.0);
    itsAttr[TRANSVCUTOFF] = makeReal("TRANSVCUTOFF", "Transverse cut-off in units of sigma",3.0);

    itsAttr[PXMULT]= makeReal("PXMULT", "Multiplier for PX",1.0);
    itsAttr[PYMULT]= makeReal("PYMULT", "Multiplier for PY",1.0);
    itsAttr[PTMULT]= makeReal("PTMULT", "Multiplier for PT",1.0);

    itsAttr[ALPHAX]= makeReal("ALPHAX", "Courant Synder parameter",1.0);
    itsAttr[ALPHAY]= makeReal("ALPHAY", "Courant Synder parameter",1.0);

    itsAttr[BETAX] = makeReal("BETAX", "Courant Synder parameter",-1.0);
    itsAttr[BETAY] = makeReal("BETAY", "Courant Synder parameter",1.0);

    itsAttr[MX]    = makeReal("MX", "Defines the distribution in x, 0+eps .. inf",1.0);
    itsAttr[MY]    = makeReal("MY", "Defines the distribution in y, 0+eps .. inf",1.0);
    itsAttr[MT]    = makeReal("MT", "Defines the distribution in t, 0+eps .. inf",1.0);

    itsAttr[DX]    = makeReal("DX", "Dispersion in x (R16 in Transport notation",0.0);
    itsAttr[DDX]   = makeReal("DDX", "First derivative of Dx",0.0);


    itsAttr[DY]    = makeReal("DY", "DY",0.0);
    itsAttr[DDY]   = makeReal("DDY", "DDY",0.0);

    itsAttr[R51]    = makeReal("R51", "R51",0.0);
    itsAttr[R52]   = makeReal("R52", "R52",0.0);

    itsAttr[R61]    = makeReal("R61", "R61",0.0);
    itsAttr[R62]   = makeReal("R62", "R62",0.0);

    itsAttr[PT] = makeReal("PT", "average longitudinal momentum",0.0);
    itsAttr[T] = makeReal("T","average longitudinal position", 0.0);

    itsAttr[SIGMAX] = makeReal("SIGMAX", "SIGMAx (m)",1.0e-2);
    itsAttr[SIGMAY] = makeReal("SIGMAY", "SIGMAy (m)",1.0e-2);
    itsAttr[SIGMAT] = makeReal("SIGMAT", "SIGMAt (m)",1.0e-2);

    itsAttr[SIGMAPX] = makeReal("SIGMAPX", "SIGMApx",0.0);
    itsAttr[SIGMAPY] = makeReal("SIGMAPY", "SIGMApy",0.0);
    itsAttr[SIGMAPT] = makeReal("SIGMAPT", "SIGMApt",0.0);

    itsAttr[CORRX] = makeReal("CORRX", "CORRx",-0.5);
    itsAttr[CORRY] = makeReal("CORRY", "CORRy", 0.5);
    itsAttr[CORRT] = makeReal("CORRT", "CORRt", 0.0);

    itsAttr[TEMISSION] = makeReal("TEMISSION", "Time in seconds in which we have emission",  0.0);
    itsAttr[NBIN]      = makeReal("NBIN", "In case of emission how many bins should we use", 0.0);
    itsAttr[DEBIN]     = makeReal("DEBIN", "Energy band for a bin in keV, defines the rebinning", 1000000.0);

    itsAttr[RISETIME]    = makeReal("RISETIME", "Rise time for GUNGAUSSFLATTOP distribution type (m)", 1.0e-2);
    itsAttr[FALLTIME]    = makeReal("FALLTIME", "Fall time for GUNGAUSSFLATTOP distribution type (m)", 1.0e-2);
    itsAttr[FLATTOPTIME] = makeReal("FLATTOPTIME", "Flat top time for GUNGAUSSFLATTOP distribution type (m)", 1.0e-2);

    itsAttr[ELASER] = makeReal("ELASER", "Laser energy (eV)", 0.0);
    itsAttr[SIGLASER] = makeReal("SIGLASER", "Sigma of (uniform) laser spot size (m)", 0.0);
    itsAttr[W] = makeReal("W", "Workfunction of material (eV)", 0.0);
    itsAttr[FE] = makeReal("FE", "Fermi energy (eV)", 0.0);
    itsAttr[AG] = makeReal("AG", "Acceleration Gradient (MV/m", 0.0);

    // Set up default beam.
    Distribution *defDistribution = clone("UNNAMED_Distribution");
    defDistribution->builtin = true;

    try {
        defDistribution->update();
        OPAL.define(defDistribution);
    } catch (...) {
        delete defDistribution;
    }
    pbin_m = NULL;

}

Distribution::Distribution(const string &name, Distribution *parent):
    Definition(name, parent),
    reference(parent->reference)
{
    pbin_m = NULL;
}

Distribution::~Distribution()
{
    if (pbin_m) {
        delete pbin_m;
        pbin_m = NULL;
    }
}

#ifdef HAVE_ENVELOPE_SOLVER
void Distribution::createSlicedBunch(SLPartBunch &p) {
    int sl = (int) Attributes::getReal(itsAttr[NBIN]);
    INFOMSG("About to create sliced a bunch with n= " << sl << endl);
}
#endif


bool Distribution::canReplaceBy(Object *object)
{
    // Can replace only by another Distribution.
    return dynamic_cast<Distribution *>(object) != 0;
}


Distribution *Distribution::clone(const string &name)
{
    return new Distribution(name, this);
}


void Distribution::execute()
{
    update();
}

void Distribution::createBinom(Vector_t emit, Vector_t alpha, Vector_t beta, Vector_t gamma,
                               Vector_t bincoef, PartBunch &beam, size_t particles,
                               bool isBinned)
{
    const double two_pi = 2.0 * 4.0*atan(1.0);

    Vector_t M;
    Vector_t PM;
    Vector_t COSCHI;
    Vector_t SINCHI;
    Vector_t CHI;
    Vector_t AMI;
    Vector_t L;
    Vector_t PL;

    size_t pc = 0;
    size_t count = 0;

    for (int i=0;i<3;i++) {
        gamma[i] *= 4.0;
        beta[i] *= 4.0;
        M[i]       =  sqrt(emit[i]*beta[i]);
        PM[i]      =  sqrt(emit[i]*gamma[i]);
        COSCHI[i]  =  sqrt(1.0/(1.0+alpha[i]*alpha[i]));
        SINCHI[i]  = -alpha[i]*COSCHI[i];
        CHI[i]     =  atan2(SINCHI[i],COSCHI[i]);
        AMI[i]     =  1.0/bincoef[i];
        L[i]       =  sqrt((bincoef[i]+1.0)/2.0)*M[i];
        PL[i]      =  sqrt((bincoef[i]+1.0)/2.0)*PM[i];
    }

    Vector_t x;
    // now copy this over to the bunch
    Vector_t p;
    double betagamma_part;
    if (isBinned)
        {
            double ekin = Attributes::getReal(itsAttr[PT]);
            double mass = beam.getM();
            double gamma_part = 1. + ekin/mass;
            betagamma_part = sqrt(ekin*ekin/(mass*mass) + 2.*ekin/mass );
            pbin_m->setGamma(gamma_part);
            *gmsg << "* Gamma = " << gamma_part << "; Beta = " << betagamma_part / gamma_part  << endl;
        }

    for (size_t n=0; n < particles; ++n) {
        double S1,S2;
        double A,AL,U,V;
        for (int i=0;i<3;i++) {
            S1=IpplRandom();
            S2=IpplRandom();
            if (bincoef[i] <= 10000)
                {
                    A=sqrt(1.0 - pow(S1,AMI[i]));
                    AL=two_pi*S2;
                    U=A*cos(AL);
                    V=A*sin(AL);
                    x[i] = L[i]*U;
                    p[i] = PL[i]*(U*SINCHI[i]+V*COSCHI[i]);
                }
            else
                {
                    A=sqrt(2.0)/2.0 * sqrt(-log(S1));
                    AL=two_pi*S2;
                    U=A*cos(AL);
                    V=A*sin(AL);
                    x[i]=M[i]*U;
                    p[i]=PM[i]*(U*SINCHI[i]+V*COSCHI[i]);
                }
            p[2] += betagamma_part;
        }
        if (pc == Ippl::myNode()) {
            if (isBinned) {
                vector<double> tmp;
                tmp.push_back(x[0]);
                tmp.push_back(x[1]);
                tmp.push_back(x[2]);
                tmp.push_back(p[0]);
                tmp.push_back(p[1]);
                tmp.push_back(p[2]);
                tmp.push_back(0);
                pbin_m->fill(tmp);
            }
            else {
                beam.create(1);
                beam.R[count] = x;
                beam.P[count] = p;
                count++;
            }
        }

        if (!isBinned)
            {
                pc++;
                if (pc == Ippl::getNodes())
                    pc = 0;
            }
    }
    if (isBinned) {
        pbin_m->sortArray();
        // now copy this over to the bunch
        // so that we can emmit the particles
        beam.setPBins(pbin_m);
    }
}

void Distribution::createUniformTUniformL(Vector_t emit, Vector_t alpha, Vector_t beta, Vector_t gamma,
                                          PartBunch &beam, size_t particles, bool isBinned)
{
    const double &two_pi = Physics::two_pi;

    Vector_t M;
    Vector_t PM;
    Vector_t COSCHI;
    Vector_t SINCHI;
    Vector_t CHI;
    Vector_t AMI;
    Vector_t L;
    Vector_t PL;

    size_t pc = 0;
    size_t count = 0;

    for (int i=0;i<3;i++) {
        gamma[i]  *= 4.0;
        beta[i]   *= 4.0;
        M[i]       =  sqrt(emit[i]*beta[i]);
        PM[i]      =  sqrt(emit[i]*gamma[i]);
        COSCHI[i]  =  sqrt(1.0/(1.0+alpha[i]*alpha[i]));
        SINCHI[i]  = -alpha[i]*COSCHI[i];
        CHI[i]     =  atan2(SINCHI[i],COSCHI[i]);
    }

    Vector_t x;
    Vector_t p;
    double betagamma_part;
    if (isBinned)
        {
            double ekin = Attributes::getReal(itsAttr[PT]);
            double mass = beam.getM();
            double gamma_part = 1. + ekin/mass;
            betagamma_part = sqrt(ekin*ekin/(mass*mass) + 2.*ekin/mass );
            pbin_m->setGamma(gamma_part);
            *gmsg << "* Gamma = " << gamma_part << "; Beta = " << betagamma_part / gamma_part  << endl;
        }

    for (size_t n=0; n < particles; ++n) {
        double S1,S2,S3,S4,S5,S6;
        double A1,AL1,AC1,U1,V1,A2,AL2,AC2,U2,V2;
        double AB = 2.0;
        while (AB > 1.0)
            {
                S1=IpplRandom();
                S2=IpplRandom();
                S3=IpplRandom();
                S4=IpplRandom();
                S5=IpplRandom();
                S6=IpplRandom();
                AL1 = two_pi*S1;
                AL2 = two_pi*S2;
                A1=sqrt(1.0 - S3);
                A2=sqrt(1.0 - S4);
                AC1 = cos(AL1);
                U1=A1*AC1;
                AC2 = cos(AL2);
                U2=A2*AC2;
                AB = sqrt(U1*U1 + U2*U2);
            }

        V1=A1*sin(AL1);
        V2=A2*sin(AL2);
        x[0] = M[0]*U1;
        x[1] = M[1]*U2;
        p[0] = PM[0]*(U1*SINCHI[0]+V1*COSCHI[0]);
        p[1] = PM[1]*(U2*SINCHI[1]+V2*COSCHI[1]);


        S1=IpplRandom();
        S2=IpplRandom();

        A1=sqrt(1.0 - S5);
        AL1 = two_pi * S6;
        U1=A1*cos(AL1);
        // now copy this over to the bunch
        V1=A1*sin(AL1);
        x[2] = M[2] * U1;
        p[2] = betagamma_part + PM[2] * (U1 * SINCHI[2] + V1 * COSCHI[2]);


        if (pc == Ippl::myNode()) {
            if (isBinned) {
                vector<double> tmp;
                tmp.push_back(x[0]);
                tmp.push_back(x[1]);
                tmp.push_back(x[2]);
                tmp.push_back(p[0]);
                tmp.push_back(p[1]);
                tmp.push_back(p[2]);
                tmp.push_back(0);
                pbin_m->fill(tmp);
            }
            else {
                beam.create(1);
                beam.R[count] = x;
                beam.P[count] = p;
                count++;
            }
        }

        if (!isBinned)
            {
                pc++;
                if (pc == Ippl::getNodes())
                    pc = 0;
            }
    }
    if (isBinned) {
        pbin_m->sortArray();
        // now copy this over to the bunch
        // so that we can emmit the particles
        beam.setPBins(pbin_m);
    }
}

void Distribution::create(PartBunch &beam, size_t Np, bool scan) {
    if (beam.getTotalNum() != 0)
        scan_m = scan;
    else
        scan_m = false; // the first time we have to create particles
    create(beam, Np);
}

void Distribution::create(PartBunch &beam, size_t Np) {

    Inform msg("Distribution::create ");

    IpplTimings::startTimer(beam.distrCreate_m);

    int ebins = (int) Attributes::getReal(itsAttr[NBIN]);
    bool isBinned = (ebins > 0);

    if (isBinned) {
        if (pbin_m)
            delete pbin_m;
        pbin_m = new PartBins((int) Attributes::getReal(itsAttr[NBIN]));
        if (scan_m) {
            beam.destroy(beam.getLocalNum(),0);
            beam.update();
            INFOMSG("In scan mode: deleted BIN structure and all particles in the bunch" << endl;);
        }
    }
    else
        pbin_m = NULL;

    beam.setTEmission(Attributes::getReal(itsAttr[TEMISSION]));
    beam.setNumBunch(1);
    const string disttype = Attributes::getString(itsAttr[DISTRIBUTION]);
    if (disttype == "GUNGAUSS" || disttype == "GUNUNIFORM" || disttype == "GUNGAUSS3D" || disttype == "GUNGAUSSFLATTOP" || disttype == "GUNGAUSSFLATTOPTH")
        // Create a an initial beam bunch that is:
        // "GUNGAUSS": uniform in space transversely and with a Gaussian ("GUNGAUS") longitudinal profile
        // "GUNUNIFORM": uniform in space transversely and longitudinally.
        // "GUNGAUSS3D": Gaussian transversely and longitudinally.
        // "GUNGAUSSFLATTOP": uniform in space transversely, a Gaussian rise and fall time longitudinally with
        //                    a uniform flattop between.
        // "GUNGAUSSFLATTOPTH": uniform in space transversely, a Gaussian rise and fall time longitudinally with
        //                    a uniform flattop between, and a transvers thermal emittance

        {
            binnDistribution(beam, Np, disttype);
        }
    else if (disttype == "BINOMIAL")
        {

            Vector_t corr(Attributes::getReal(itsAttr[CORRX]),
                          Attributes::getReal(itsAttr[CORRY]),
                          Attributes::getReal(itsAttr[CORRT]));

            Vector_t sigX(Attributes::getReal(itsAttr[SIGMAX]),
                          Attributes::getReal(itsAttr[SIGMAY]),
                          Attributes::getReal(itsAttr[SIGMAT]));

            Vector_t sigPX(eVtoBetaGamma(Attributes::getReal(itsAttr[SIGMAPX]),beam.getM()),
                           eVtoBetaGamma(Attributes::getReal(itsAttr[SIGMAPY]),beam.getM()),
                           eVtoBetaGamma(Attributes::getReal(itsAttr[SIGMAPT]),beam.getM()));

            Vector_t bincoef(Attributes::getReal(itsAttr[MX]),
                             Attributes::getReal(itsAttr[MY]),
                             Attributes::getReal(itsAttr[MT]));

            Vector_t emit;
            Vector_t alpha;
            Vector_t beta;
            Vector_t gamma;

            for(int j=0; j<3; j++) {
                double chi = asin(corr[j]);
                emit[j] = sigX[j]*sigPX[j]*cos(chi);
            }
            for(int j=0; j<3; j++) {
                beta[j]  = sigX[j]*sigX[j]/emit[j];
                gamma[j] = sigPX[j]*sigPX[j]/emit[j];
                alpha[j] = -corr[j]*sqrt(beta[j]*abs(gamma[j]));
            }


            createBinom(emit, alpha, beta, gamma, bincoef, beam, Np, isBinned);
        }
    else if (disttype == "UNITUNIL")
        {

            Vector_t corr(Attributes::getReal(itsAttr[CORRX]),
                          Attributes::getReal(itsAttr[CORRY]),
                          Attributes::getReal(itsAttr[CORRT]));

            Vector_t sigX(Attributes::getReal(itsAttr[SIGMAX]),
                          Attributes::getReal(itsAttr[SIGMAY]),
                          Attributes::getReal(itsAttr[SIGMAT]));

            Vector_t sigPX(eVtoBetaGamma(Attributes::getReal(itsAttr[SIGMAPX]),beam.getM()),
                           eVtoBetaGamma(Attributes::getReal(itsAttr[SIGMAPY]),beam.getM()),
                           eVtoBetaGamma(Attributes::getReal(itsAttr[SIGMAPT]),beam.getM()));

            Vector_t emit;
            Vector_t alpha;
            Vector_t beta;
            Vector_t gamma;

            for(int j=0; j<3; j++) {
                double chi = asin(corr[j]);
                emit[j] = sigX[j]*sigPX[j]*cos(chi);
            }
            for(int j=0; j<3; j++) {
                beta[j]  = sigX[j]*sigX[j]/emit[j];
                gamma[j] = sigPX[j]*sigPX[j]/emit[j];
                alpha[j] = -corr[j]*sqrt(beta[j]*abs(gamma[j]));
            }

            createUniformTUniformL(emit, alpha, beta, gamma, beam, Np, isBinned);
        }
    else if (disttype == "GAUSS")
        {

            double corr[7];
            corr[0] = Attributes::getReal(itsAttr[CORRX]);
            corr[1] = Attributes::getReal(itsAttr[CORRY]);
            corr[2] = Attributes::getReal(itsAttr[CORRT]);
            corr[3] = Attributes::getReal(itsAttr[R61]);
            corr[4] = Attributes::getReal(itsAttr[R62]);
            corr[5] = Attributes::getReal(itsAttr[R51]);
            corr[6] = Attributes::getReal(itsAttr[R52]);
            double Hs2a = Attributes::getReal(itsAttr[SIGMAX]);
            double Hs2b = eVtoBetaGamma(Attributes::getReal(itsAttr[SIGMAPX]),beam.getM()); //in eV

            double Vs2a = Attributes::getReal(itsAttr[SIGMAY]);
            double Vs2b = eVtoBetaGamma(Attributes::getReal(itsAttr[SIGMAPY]),beam.getM()); //in eV

            double Ls2a = Attributes::getReal(itsAttr[SIGMAT]);
            double Ls2b = eVtoBetaGamma(Attributes::getReal(itsAttr[SIGMAPT]),beam.getM()); //in eV

            double avrgpt = eVtoBetaGamma(Attributes::getReal(itsAttr[PT]),beam.getM());
            double avrgt = Attributes::getReal(itsAttr[T]);

            RANLIB_class *rGen = new RANLIB_class(265314159,4);

            unsigned int pc = 0;
            unsigned int count = 0;

            if (beam.isZPeriodic())
                INFOMSG("Distribution: uniform in z (periodic BC) z = " << -0.5*beam.getGaBeLa() << " ... " << 0.5*beam.getGaBeLa() << endl);
            for(int i=0; i<Np; i++) {
                double x,y;       // generate independent Gaussians, then correlate and finaly scale

                x  = rGen->gauss(0.0,1.0);
                y  = rGen->gauss(0.0,1.0);
                double xx =x;
                double yy =y;
                double px0  = x*corr[0] + y*sqrt(1.0 - corr[0]*corr[0]);
                double x0   =  x*Hs2a;
                px0 *= Hs2b;

                x  = rGen->gauss(0.0,1.0);
                y  = rGen->gauss(0.0,1.0);
                double py0  = x*corr[1] + y*sqrt(1.0 - corr[1]*corr[1]);
                double y0   =  x*Vs2a;
                py0 *= Vs2b;

                double del0;
                double psi0;

                if (beam.isZPeriodic()) {
                    /*
                      create uniform distribution
                    */
                    //del0 = (IpplRandom()-0.5)*p[2];
                    //	psi0 = (IpplRandom()*beam.getGaBeLa()) - (0.5*beam.getGaBeLa());

                } else {
                    x  = rGen->gauss(0.0,1.0);
                    y  = rGen->gauss(0.0,1.0);
                    double l32=(corr[6]-corr[0]*corr[5])/sqrt(1.0 - corr[0]*corr[0]);
                    double l33=sqrt(1-corr[5]*corr[5]-l32*l32);
                    psi0 = xx*corr[5]+yy*l32+x*l33;
                    double l42=(corr[4]-corr[0]*corr[3])/sqrt(1.0 - corr[0]*corr[0]);
                    double l43=(corr[2]-corr[5]*corr[3]-l42*l32)/l33;
                    double l44=sqrt(1-corr[3]*corr[3]-l42*l42-l43*l43);
                    del0 = xx*corr[3]+yy*l42+x*l43+y*l44;
                    //      	  del0 = xx*corr[3]+yy*0.268109486023689+y*0.693221684242576;
                    psi0 = avrgt + psi0*Ls2a;
                    del0 = avrgpt + Ls2b*del0;
                }
                if (pc == Ippl::myNode()) {
                    beam.create(1);
                    beam.R[count] = Vector_t(x0,y0,psi0);
                    beam.P[count] = Vector_t(px0,py0,del0);
                    beam.Bin[count] = 0; // not initialized
                    //           dist.precision(8);
                    //           dist << x0 << "\t"
                    //                << y0 << "\t"
                    //                << psi0 << "\t"
                    //                << px0 << "\t"
                    //                << py0 << "\t";
                    //           dist.precision(15);
                    //           dist << del0 << endl;
                    count++;
                }
                pc++;
                if (pc == Ippl::getNodes())
                    pc = 0;
            }
        }
    else if (disttype == "UNIFORMXYZ")
        {

            double corr[3];

            corr[0] = Attributes::getReal(itsAttr[CORRX]);
            corr[1] = Attributes::getReal(itsAttr[CORRY]);
            corr[2] = Attributes::getReal(itsAttr[CORRT]);

            double Hs2a = Attributes::getReal(itsAttr[SIGMAX]);
            double Hs2b = Attributes::getReal(itsAttr[SIGMAPX]);

            double Vs2a = Attributes::getReal(itsAttr[SIGMAY]);
            double Vs2b = Attributes::getReal(itsAttr[SIGMAPY]);

            double Ls2a = Attributes::getReal(itsAttr[SIGMAT]);
            double Ls2b = Attributes::getReal(itsAttr[SIGMAPT]);

            double nBins = Attributes::getReal(itsAttr[NBIN]);

            RANLIB_class *rGen = new RANLIB_class(265314159,4);

            unsigned int pc = 0;
            unsigned int count = 0;

            double gamma = 1 + (Ls2b / beam.getM());
            double beta  = sqrt(1- (1/(gamma*gamma)));
            double bega  = beta*gamma;

            for(int i=0; i<Np; i++) {

                double x,y,z,px,py,pz;       // generate independent Gaussians, then correlate and finaly scale
                double R =3;

                while(R >1){

                    x = rGen->uniform(-1.0,1.0);
                    y = rGen->uniform(-1.0,1.0);
                    z = rGen->uniform(-1.0,1.0);

                    px = rGen->uniform(-1.0,1.0);
                    py = rGen->uniform(-1.0,1.0);
                    pz = rGen->uniform(-1.0,1.0);

                    R = sqrt(x*x + y*y + z*z);

                    // or can generate uniform distribution in 6D phase space by
                    // using following line instead of above one.
                    // R = sqrt(x*x + y*y + z*z + px*px + py*py +pz*pz);
                }

                px  = x*corr[0] + px*sqrt(1.0 - corr[0]*corr[0]);
                x   = x*Hs2a;
                px *= Hs2b;

                py  = y*corr[1] + py*sqrt(1.0 - corr[1]*corr[1]);
                y   = y*Vs2a;
                py *= Vs2b;

                pz  = z*corr[2] + pz*sqrt(1.0 - corr[2]*corr[2]);
                z   = z*Ls2a;
                pz *= Ls2b;

                if (pc == Ippl::myNode()) {
                    beam.create(1);
                    beam.R[count] = Vector_t(x,y,z);
                    beam.P[count] = Vector_t(px,py,pz );
                    beam.Bin[count] = 0; // not initialized
                    count++;
                }
                pc++;
                if (pc == Ippl::getNodes())
                    pc = 0;
            }

        }
    else if (disttype == "FROMFILE")
        {

            *gmsg <<"---------------------------------------------"<<endl;
            *gmsg <<"     READ ININITAL DISTRIBUTION FROM FILE    "<<endl;
            const string filename=Attributes::getString(itsAttr[FNAME]);

            double x0,px0,y0,py0,psi0,del0;

            std::ifstream fs;
            fs.open(filename.c_str());

            if(fs.fail()){
                throw OpalException("Distribution::Create()",
                                    "Open file operation failed, please check if \""
                                    + filename +  "\" really exists.");
            }

            fs >> Np;
            if ( Np <= 0 ){
                throw OpalException("Distribution::Create()",
                                    " The particle number should be bigger than zero! Please check the first line of file \""
                                    + filename +  "\".");
            }
            unsigned int pc = 0;
            unsigned int count = 0;
            // read in particles number from the first line of the file.

            for(unsigned int i=0; i<Np; i++) {

                if ( !fs.eof() ){
                    if (pc == Ippl::myNode()) {
                        // create 1 particle
                        beam.create(1);

                        fs>> x0 >> px0 >> y0 >> py0 >> psi0 >> del0;
                        // reasd in r and p
                        beam.R[count] = Vector_t(x0,y0,psi0);
                        beam.P[count] = Vector_t(px0,py0,del0);
                        beam.Bin[count] = 0; // not initialized
                        count++;
                    }
                }
                else {
                    throw OpalException("Distribution::Create()",
                                        "End of file reached before all particles imported, please check file \""
                                        + filename +  "\".");
                    return;
                }

            }
            if (count != Np) {

                throw OpalException("Distribution::Create()",
                                    "The number of actually read particles does NOT match Np (number of first line), please check file \""
                                    + filename +  "\".");
            }

            *gmsg <<"                     DONE !                  "<<endl;
            *gmsg <<"---------------------------------------------"<<endl;
            fs.close();

        }

    /*
      In the case of a binned distribution (gun)
      we have to do the boundp after emission.
    */
    if (!(isBinned)) {
        beam.boundp();
        beam.LastSection=0;
    }
    IpplTimings::stopTimer(beam.distrCreate_m);
}


void Distribution::binnDistribution(PartBunch &beam, size_t Np, string distType)
{
    const double &two_pi = Physics::two_pi;
    unsigned int pc = 0;

    double dEBins = Attributes::getReal(itsAttr[DEBIN]);
    pbin_m->setRebinEnergy(dEBins);

    double eneg= Attributes::getReal(itsAttr[PT]);
    double gamma = 1. + (eneg / beam.getM());
    double beta  = sqrt(1. - (1./(gamma*gamma)));
    double bega  = beta*gamma;

    pbin_m->setGamma(gamma); //has to be done by all processors!

    double corr[3] = {Attributes::getReal(itsAttr[CORRX]),
                      Attributes::getReal(itsAttr[CORRY]),
                      Attributes::getReal(itsAttr[CORRT])};

    double transvCutOff = Attributes::getReal(itsAttr[TRANSVCUTOFF]);
    double Hs2a = Attributes::getReal(itsAttr[SIGMAX]);
    double Hs2b = eVtoBetaGamma(Attributes::getReal(itsAttr[SIGMAPX]),beam.getM());

    double Vs2a = Attributes::getReal(itsAttr[SIGMAY]);
    double Vs2b = eVtoBetaGamma(Attributes::getReal(itsAttr[SIGMAPY]),beam.getM());

    double Ls2a = Attributes::getReal(itsAttr[SIGMAT]);
    double Ls2b = eVtoBetaGamma(Attributes::getReal(itsAttr[SIGMAPT]),beam.getM());

    double riseTime = Attributes::getReal(itsAttr[RISETIME]);
    double fallTime = Attributes::getReal(itsAttr[FALLTIME]);
    double flatTop = Attributes::getReal(itsAttr[FLATTOPTIME]);

    int particlesFront = Np * sqrt(2.0 * Physics::pi) * (riseTime / 2.0)
        / (flatTop + sqrt(2.0 * Physics::pi) * (riseTime + fallTime) / 2.0);
    int particlesTail = Np * sqrt(2.0 * Physics::pi) * (fallTime / 2.0)
        / (flatTop + sqrt(2.0 * Physics::pi) * (riseTime + fallTime) / 2.0);

    double nBins = Attributes::getReal(itsAttr[NBIN]);

    RANLIB_class *rGen = new RANLIB_class(265314159,4);

    /*
      prepare quantities for thermal emittance calculation

    */

    double workf = 0.0;         // eV
    double siglaser = 0.0;      // m
    double elaser = 0.0;        // eV
    double fe = 0.0;            // Fermi energy eV
    double ag = 0.0;            // Acceleration gradient eV/m
    double ekin = 0.0;          // eV
    double phimax = 0.0;        // rad
    double schottky = 0.0;      // eV
    double ptot = 0.0;          // beta gamma
    std::ofstream os;

    if (distType == "GUNGAUSSFLATTOPTH") {
        workf = Attributes::getReal(itsAttr[W]);
        siglaser = Attributes::getReal(itsAttr[SIGLASER]);
        elaser = Attributes::getReal(itsAttr[ELASER]);
        fe = Attributes::getReal(itsAttr[FE]);
        ag = Attributes::getReal(itsAttr[AG]) * 1000000;

        ekin = fe + elaser;

        schottky = sqrt(Physics::q_e*ag/4/Physics::pi/Physics::epsilon_0);

        workf -= schottky;

        phimax = acos(sqrt((fe+workf)/ekin));

        ptot = eVtoBetaGamma((workf+elaser-fe-schottky),beam.getM()),


            *gmsg << " ---------------------------------------------------------------" << endl;
        *gmsg << " ---------------------I N P U T --------------------------------" << endl;
        *gmsg << " GUNGAUSS FLAT TOP &  THERMAL EMITTANCE " << endl;
        *gmsg << " Laser energy = " << elaser << " [eV]  rrms= " << siglaser << endl;
        *gmsg << " Work function = " << Attributes::getReal(itsAttr[W]) << " [eV]  " << endl;
        *gmsg << " Fermi energy = " << fe << " [eV]  " << endl;
        *gmsg << " Shottky energy = " << schottky << " [eV]  " << endl;
        *gmsg << " ---------------------------------------------------------------" << endl;
        *gmsg << " -------D E R I V E D  Q U A N T I T I E S ---------------------" << endl;
        *gmsg << " Kinetic energy = " << ekin << " [eV]  " << endl;
        *gmsg << " Effective Work function = " << workf << " [eV]  " << endl;
        *gmsg << " Phi max = " << phimax*180/Physics::pi << " [deg]  " << endl;
        *gmsg << " Ptot = " << workf+elaser-fe-schottky << " [eV]  " << endl;

        *gmsg << " ---------------------------------------------------------------" << endl;
    }

    if (Ippl::getNodes() == 1) {
        *gmsg << " Write distribution to file dist.dat" << endl;
        string file("dist.dat");
        os.open(file.c_str());
        if (os.bad()) {
            *gmsg << "Unable to open output file " <<  file << endl;
        }
        os << "# x px y py z pz phi theta Ethx [eV] Ethy [eV] Ethz [eV]" << endl;
    }


    for(int i=0; i<Np; i++) {
        double x,y;       // generate independent Gaussians, then correlate and finally scale
        double u1, u2;
        double xy = 6;

        if (distType != "GUNGAUSS3D")
            {
                while (xy > 1){
                    x  = rGen->uniform(-1.0,1.0);
                    y  = rGen->uniform(-1.0,1.0);
                    xy = sqrt(x*x + y*y);
                }
            }
        else
            {
                while (xy > transvCutOff || isnan(xy))
                    {
                        u1 = rGen->uniform(0.0,1.0);
                        u2 = rGen->uniform(0.0,1.0);
                        x = sqrt(-2. * log(1. - u1)) * cos(two_pi * u2);
                        y = sqrt(-2. * log(1. - u1)) * sin(two_pi * u2);
                        xy = sqrt(x*x + y*y);
                    }
            }

        double x0   =  x*Hs2a;
        double y0   =  y*Vs2a;

        double del0 = 0.0;
        double psi0 = 0.0;

        if (distType == "GUNGAUSS" || distType == "GUNGAUSS3D") {
            x  = rGen->gauss(0.0,1.0);
            y  = rGen->gauss(0.0,1.0);
        } else if (distType == "GUNUNIFORM") {
            x = rGen->uniform(0.0, 1.0);
            y = rGen->uniform(0.0, 1.0);
        } else if (distType == "GUNGAUSSFLATTOP" || distType == "GUNGAUSSFLATTOPTH") {
            if (i < particlesFront) {
                // Fill rise.
                psi0 = -1.0;
                while (psi0 < 0.0)
                    psi0 = rGen->gauss(0.0, 1.0);
                psi0 = psi0 * riseTime + flatTop;
            } else if (i >= particlesFront && i < particlesFront + particlesTail) {
                // Fill fall.
                psi0 = 1.0;
                while (psi0 > 0.0)
                    psi0 = rGen->gauss(0.0, 1.0);
                psi0 = psi0 * fallTime;
            } else {
                // Fill flat top.
                psi0 = rGen->uniform(0.0, 1.0);
                psi0 *= flatTop;
            }
        }

        if (distType != "GUNGAUSSFLATTOP" && distType != "GUNGAUSSFLATTOPTH") {
            del0  = x*corr[2] + y*sqrt(1.0 - corr[2]*corr[2]);
            psi0  = x*Ls2a;
        }

        /*
          This is the longitudinal momenta in units of beta gamma.
        */

        del0  = bega + Ls2b*del0;

        double px0 = 0.0;
        double py0 = 0.0;
        double phi = 0.0;
        double theta = 0.0;

        if (distType == "GUNGAUSSFLATTOPTH") {

            /*
              Now calculate the thermal emittance
            */

            phi   = phimax * rGen->uniform(0.0, 1.0);
            theta = two_pi * rGen->uniform(0.0, 1.0);

            px0 = ptot * sin(phi) * cos(theta);
            py0 = ptot * sin(phi) * sin(theta);
            del0 = ptot * cos(phi);

        }

        if (pc == Ippl::myNode()) {
            vector<double> tmp;
            tmp.push_back(x0);
            tmp.push_back(y0);
            tmp.push_back(psi0);
            tmp.push_back(px0);
            tmp.push_back(py0);
            tmp.push_back(del0);
            tmp.push_back(0);
            pbin_m->fill(tmp);
        }
        pc++;
        if (pc == Ippl::getNodes())
            pc=0;

        if (Ippl::getNodes() == 1) {
            os << x0 << "\t " << px0    << "\t "
               << y0 << "\t " << py0    << "\t "
               << psi0 << "\t " << del0 << "\t "
               << phi*180./Physics::pi << "\t " << theta*180./Physics::pi << "\t "
               << betaGammatoeV(px0, beam.getM())  << "\t "
               << betaGammatoeV(py0, beam.getM())  << "\t "
               << betaGammatoeV(del0, beam.getM()) << "\t " << endl;
        }

    }
    if (Ippl::getNodes() == 1)
        os.close();

    pbin_m->sortArray();
    // now copy this over to the bunch
    // so that we can emit the particles
    beam.setPBins(pbin_m);
}

void Distribution::doRestart(PartBunch &beam, size_t Np, size_t restartStep)
{
    H5PartFile *H5file;
    string fn;

    IpplTimings::startTimer(beam.distrReload_m);

    if (OPAL.hasRestartFile()){
        fn=OPAL.getRestartFileName();
        *gmsg<<"Restart from a specified file:"<<fn<<endl;

    }
    else{
        fn= OPAL.getInputFn();
        int pos=fn.find(string("."),0);
        fn.erase(pos,fn.size()-pos);
        //        beam.setTEmission(Attributes::getReal(itsAttr[TEMISSION]));
        fn += string(".h5");
    }

#ifdef PARALLEL_IO
    H5file=H5PartOpenFileParallel(fn.c_str(),H5PART_READ,MPI_COMM_WORLD);
#else
    H5file=H5PartOpenFile(fn.c_str(),H5PART_READ);
#endif

    if(!H5file) {
        ERRORMSG("could not open file '" << fn << "';  exiting!" << endl);
        exit(0);
    }
    if (restartStep == -1) {
        restartStep = H5PartGetNumSteps(H5file) - 1;
        OPAL.setRestartStep(restartStep);
    } else {
        if (restartStep != H5PartGetNumSteps(H5file) - 1 && !OPAL.hasRestartFile()) {
            ERRORMSG("can't append to the file '" << fn << "' exiting!" << endl);
            exit(0);
        }
    }


    *gmsg << "restart step = " << restartStep << endl;
    H5PartSetStep(H5file,restartStep);
    int N=(int)H5PartGetNumParticles(H5file);

    h5part_int64_t totalSteps = H5PartGetNumSteps(H5file);
    *gmsg << "total number of particles = " << N << " total steps " << totalSteps << endl;

    //TODO: do a more sophisticated distribution of particles?
    //my guess is that the end range index is EXCLUSIVE!
    int numberOfParticlesPerNode = (int) floor((double) N / Ippl::getNodes());
    long long starti = Ippl::myNode() * numberOfParticlesPerNode;
    long long endi = 0;
    // ensure that we don't miss any particle in the end
    if(Ippl::myNode() == Ippl::getNodes() - 1)
        endi = -1;
    else
        endi = starti + numberOfParticlesPerNode;

    // cout << "Node " << Ippl::myNode() << " gets particles " << starti << " - " << endi << endl;

    H5PartSetView(H5file,starti,endi);
    N = (int)H5PartGetNumParticles(H5file);

    // cout << "Node " << Ippl::myNode() << " has " << N << " particles" << endl;

    double actualT;
    H5PartReadStepAttrib(H5file,"TIME",&actualT);

    beam.setT(actualT);

    /*
     * eventually read more attributes a la IMPACT-T
     */

    void *varray = malloc(N*sizeof(double));
    double *farray = (double*)varray;

    beam.create(N);

    H5PartReadDataFloat64(H5file,"x",farray);
    for (unsigned long int n=0; n < N; ++n) {
        beam.R[n](0)=farray[n];
        beam.Bin[n] = 0; // not initialized
    }
    H5PartReadDataFloat64(H5file,"y",farray);
    for (unsigned long int n=0; n < N; ++n)
        beam.R[n](1)=farray[n];

    H5PartReadDataFloat64(H5file,"z",farray);
    for (unsigned long int n=0; n < N; ++n)
        beam.R[n](2)=farray[n];

    H5PartReadDataFloat64(H5file,"px",farray);
    for (unsigned long int n=0; n < N; ++n)
        beam.P[n](0)=farray[n];

    H5PartReadDataFloat64(H5file,"py",farray);
    for (unsigned long int n=0; n < N; ++n)
        beam.P[n](1)=farray[n];

    H5PartReadDataFloat64(H5file,"pz",farray);
    for (unsigned long int n=0; n < N; ++n)
        beam.P[n](2)=farray[n];

    if(farray)
        free(farray);

    Ippl::Comm->barrier();
    H5PartCloseFile(H5file);
    beam.boundp();
    beam.LastSection=0;
    IpplTimings::stopTimer(beam.distrReload_m);
}

void Distribution::doRestart_cycl(PartBunch &beam, size_t Np, size_t restartStep, const int specifiedNumBunch)
{
    IpplTimings::startTimer(beam.distrReload_m);
    *gmsg<<"---------------- Start reading hdf5 file----------------"<<endl;
    H5PartFile *H5file;

    string fn;
    if (OPAL.hasRestartFile()){

        fn=OPAL.getRestartFileName();
        *gmsg<<"Restart from a specified file:"<<fn<<endl;
    }
    else{
        fn= OPAL.getInputFn();
        int pos=fn.find(string("."),0);
        fn.erase(pos,fn.size()-pos);

        beam.setTEmission(Attributes::getReal(itsAttr[TEMISSION]));

        fn += string(".h5");
    }

    *gmsg<<"Restart from hdf5 format file "<<fn<<", read phase space data of DumpStep "<<restartStep<<endl;

#ifdef PARALLEL_IO
    H5file=H5PartOpenFileParallel(fn.c_str(),H5PART_READ,MPI_COMM_WORLD);
#else
    H5file=H5PartOpenFile(fn.c_str(),H5PART_READ);
#endif

    if(!H5file) {
        ERRORMSG("File open failed:  exiting!" << endl);
        exit(0);
    }

    H5PartSetStep(H5file,restartStep);
    const int globalN=(int)H5PartGetNumParticles(H5file);

    h5part_int64_t totalSteps = H5PartGetNumSteps(H5file);
    *gmsg << "total number of particles = " <<globalN <<endl;

    int numberOfParticlesPerNode = (int) floor((double) globalN / Ippl::getNodes());
    long long starti = Ippl::myNode() * numberOfParticlesPerNode;
    long long endi = 0;

    if(Ippl::myNode() == Ippl::getNodes() - 1)
        endi = -1;
    else
        endi = starti + numberOfParticlesPerNode;

    H5PartSetView(H5file,starti,endi);
    const int localN = (int)H5PartGetNumParticles(H5file);

    // debug
    // Inform *gmsgAll;
    // gmsgAll = new  Inform("Message",INFORM_ALL_NODES);
    // *gmsgAll<< "total number of particles on this node = " <<localN <<endl;
    // end debug

    double actualT;
    H5PartReadStepAttrib(H5file,"TIME",&actualT);

    beam.setT(actualT);

    double lpath;
    H5PartReadStepAttrib(H5file,"LPATH",&lpath);
    beam.setLPath(lpath);

    h5part_int64_t tstep;
    H5PartReadStepAttrib(H5file,"TrackStep",&tstep);
    beam.setTrackStep((long long)tstep);

    h5part_int64_t numBunch;
    H5PartReadStepAttrib(H5file,"NumBunch",&numBunch);
    beam.setNumBunch((int)numBunch);
    *gmsg<<"There are "<<numBunch<<" Bunches(bins) exist in this file"<<endl;

    double gammaBin[numBunch];

    H5PartReadStepAttrib(H5file,"GammaBin",&gammaBin);

    void *varray = malloc(localN*sizeof(double));
    double *farray = (double*)varray;

    h5part_int64_t *larray = (h5part_int64_t *)varray;

    beam.create(localN);

    H5PartReadDataFloat64(H5file,"x",farray);
    for (unsigned long int n=0; n < localN; ++n)
        beam.R[n](0)=farray[n];

    H5PartReadDataFloat64(H5file,"y",farray);
    for (unsigned long int n=0; n < localN; ++n)
        beam.R[n](1)=farray[n];

    H5PartReadDataFloat64(H5file,"z",farray);
    for (unsigned long int n=0; n < localN; ++n)
        beam.R[n](2)=farray[n];

    H5PartReadDataFloat64(H5file,"px",farray);
    for (unsigned long int n=0; n < localN; ++n)
        beam.P[n](0)=farray[n];

    H5PartReadDataFloat64(H5file,"py",farray);
    for (unsigned long int n=0; n < localN; ++n)
        beam.P[n](1)=farray[n];

    H5PartReadDataFloat64(H5file,"pz",farray);
    for (unsigned long int n=0; n < localN; ++n)
        beam.P[n](2)=farray[n];

    H5PartReadDataInt64(H5file,"id",larray);
    for (unsigned long int n=0; n < localN; ++n)
        beam.ID[n]=larray[n];

    // only for multi-bunch mode
    if ( specifiedNumBunch > 1)
        {
            size_t partInBin[numBunch];
            for(int ii=0; ii<numBunch; ii++) partInBin[ii] = 0 ;

            // assign bin index for each particle and
            // calculate total particles number for each bin
            for (unsigned long int n=0; n < localN; ++n)
                {
                    double deltgamma[numBunch];
                    double gamma = sqrt(1.0 + dot(beam.P[n], beam.P[n]));
                    int index = 0;

                    for(int ii=0; ii<numBunch; ii++)
                        deltgamma[ii] = abs( gammaBin[ii] - gamma );

                    for(int ii=0; ii<numBunch; ii++)
                        if( *(deltgamma+index) > *(deltgamma+ii) )
                            index = ii;

                    beam.Bin[n]=index;
                    partInBin[index]++;
                }

            for(int ii=0; ii<numBunch; ii++)
                reduce(partInBin[ii],partInBin[ii],OpAddAssign());

            // instantiate PartBins class for restart run
            beam.setPBins( new PartBins(specifiedNumBunch, (int)numBunch, partInBin ));
        }


    if(farray) free(farray);

    Ippl::Comm->barrier();
    H5PartCloseFile(H5file);
    beam.boundp();
    *gmsg<<"----------------Finish reading hdf5 file----------------"<<endl;
    IpplTimings::stopTimer(beam.distrReload_m);
}


Distribution *Distribution::find(const string &name)
{
    Distribution *beam = dynamic_cast<Distribution*>(OPAL.find(name));

    if (beam == 0) {
        throw OpalException("Distribution::find()", "Distribution \"" + name + "\" not found.");
    }

    return beam;
}

/*
  double Distribution::getET() const
  {
  return Attributes::getReal(itsAttr[ET]);
  }

  void Distribution::setET(double value)
  {
  Attributes::setReal(itsAttr[ET], value);
  }

*/

const PartData &Distribution::getReference() const
{
    // Cast away const, to allow logically constant Distribution to update.
    const_cast<Distribution *>(this)->update();
    return reference;
}

void Distribution::update()
{

}


void Distribution::tfsDescriptors(std::ostream &os) const
{
    os << "@ Distribution     %s  " << getOpalName() << '\n' ;
}
// Emittances added to tfsDescriptors by JMJ 4/4/2000


Inform &Distribution::print(Inform &os) const {
    double transvCutOff = Attributes::getReal(itsAttr[TRANSVCUTOFF]);
    os << "* ************* D I S T R I B U T I O N ******************************************** " << endl;
    if (!OPAL.inRestartRun()) {
        os << "* Distribution:\t" << getOpalName() << endl;
        os << "* Distribution type:\t" << Attributes::getString(itsAttr[DISTRIBUTION]) << endl;
        os << "* sigmax=\t" << Attributes::getReal(itsAttr[SIGMAX]);
        os << "\t sigmay=\t" << Attributes::getReal(itsAttr[SIGMAY]);
        os << "\t sigmat=\t" << Attributes::getReal(itsAttr[SIGMAT]) << endl;
        os << "* transverse cut  " << transvCutOff << " [sigma] " << endl;
        os << "* sigmapx=\t" << Attributes::getReal(itsAttr[SIGMAPX]);
        os << "\t sigmapy=\t" << Attributes::getReal(itsAttr[SIGMAPY]);
        os << "\t pt +- sigmapt=\t" << Attributes::getReal(itsAttr[PT]) << "+-" << Attributes::getReal(itsAttr[SIGMAPT]) << endl;

        os << "* corr x-px=\t" << Attributes::getReal(itsAttr[CORRX]);
        os << "\t corr y=py=\t" << Attributes::getReal(itsAttr[CORRY]);
        os << "\t corr t-pt=\t" << Attributes::getReal(itsAttr[CORRT]) << endl;
        if(Attributes::getReal(itsAttr[TEMISSION]) > 0.0) {
            os << "* -------- G U N -------------------------------------------------------------------" << endl;
            os << "* temission [s] \t" << Attributes::getReal(itsAttr[TEMISSION]);
            os << "\t n bins \t" << (int) Attributes::getReal(itsAttr[NBIN]) << endl;
            os << "\t dEbin [keV] \t" << Attributes::getReal(itsAttr[DEBIN]) << endl;
            os << "* ----------------------------------------------------------------------------------" << endl;
        }
    }
    else
        os << "* Distribution from restart file" << endl;

    os << "* ********************************************************************************** " << endl;
}

/**
   For testing the implementation of particle bins do the following:

   Compile: g++ -DPartBinTest -o Bins Bins.cpp ranlib.o

   Run test: ./Bins 10 10000 1.0

*/


PartBins::PartBins(int bins) :
    bins_m(bins),
    nemittedBins_m(0),
    xmin_m(0.0),
    xmax_m(0.0),
    gamma_m(1.0)
{

    nBin_m = new size_t[bins_m];
    xbinmin_m = new double[bins_m];
    xbinmax_m = new double[bins_m];
    binsEmitted_m = new bool[bins_m];
    nDelBin_m = new size_t[bins_m];

    for (int i=0; i<bins_m; i++) {
        nDelBin_m[i] = nBin_m[i] = 0;
        xbinmax_m[i] = -DBL_MAX;
        xbinmin_m[i] = DBL_MAX;
        binsEmitted_m[i]=false;
    }
}

// constructer function for cyclotron
// for the time being, the total existing bins number equals the total existing bunches number
PartBins::PartBins(int specifiedNumBins, int bins, size_t  partInBin[]) :
    bins_m(specifiedNumBins),
    nemittedBins_m(bins)
{

    binsEmitted_m = new bool[bins_m];
    nBin_m = new size_t[bins_m];
    nDelBin_m = new size_t[bins_m];


    //dummy arrays, not used
    xbinmin_m = new double;
    xbinmax_m = new double;


    for (int i=0; i<bins_m; i++) {
        nDelBin_m[i] = 0;
        nBin_m[i] = 0;
        binsEmitted_m[i]=false;
    }

    for (int i=0; i<nemittedBins_m; i++) {
        nBin_m[i] = partInBin[i];
        *gmsg<<"Read in: Bin="<<i<<" Particles Num="<<nBin_m[i]<<endl;
        binsEmitted_m[i]=true;
    }

}

size_t PartBins::getTotalNum() {
    size_t s = 0;
    size_t sd = 0;
    size_t gs = 0;

    for (int i=0; i<getLastemittedBin(); i++) {
        s  += nBin_m[i];
        sd += nDelBin_m[i];
    }
    gs = s-sd;
    reduce(gs, gs, OpAddAssign());
    return gs;
}

size_t PartBins::getTotalNumPerBin(int b) {
    size_t s = 0;
    s  = nBin_m[b];
    reduce(s, s, OpAddAssign());
    return s;
}

void PartBins::updateStatus(const int bunchCount, const size_t partInBin){
    // array index of binsEmitted_m[] starts from 0
    // nemittedBins_m and bins_m index starts from 1
    binsEmitted_m[bunchCount-1] = true;
    size_t NpartInBin = partInBin;
    reduce(NpartInBin, NpartInBin, OpAddAssign());
    nBin_m[bunchCount-1] = NpartInBin;
    nemittedBins_m++;
}

void PartBins::updateDeletedPartsInBin(size_t countLost[])
{
    Inform m2all("updateDeletedPartsInBin ", INFORM_ALL_NODES);

    for (int ii=0; ii < getLastemittedBin(); ii++){
        bool flagNeedUpdate = false;
        flagNeedUpdate = (countLost[ii] > 0);

        reduce(flagNeedUpdate, flagNeedUpdate, OpOr());

        if(flagNeedUpdate){
            reduce(countLost[ii], countLost[ii], OpAddAssign());
            nDelBin_m[ii] = countLost[ii];
        }
        m2all << "In Bin: "<< ii <<", "<< nDelBin_m[ii] <<" particle(s) lost"<<endl;
    }
}

void PartBins::updatePartInBin(size_t countLost[]) {

    Inform msg0("updatePartInBin ");

    //  for (int ii=0; ii < bins_m; ii++){
    for (int ii=0; ii < nemittedBins_m; ii++){
        msg0 << "In Bin: "<< ii <<", "<< nBin_m[ii] <<" particles "<<endl;
    }

    for (int ii=0; ii < nemittedBins_m; ii++){
        bool flagNeedUpdate = false;
        if (countLost[ii] > 0) flagNeedUpdate = true;

        reduce(flagNeedUpdate, flagNeedUpdate, OpOr());

        if(flagNeedUpdate){
            reduce(countLost[ii], countLost[ii], OpAddAssign());
            nBin_m[ii] -= countLost[ii];
            msg0 << "In Bin: "<< ii <<", "<< countLost[ii] <<" particle(s) lost"<<endl;
        }
    }
}

void PartBins::resetPartInBin( size_t newPartNum[])
{
    for(int ii=0; ii<nemittedBins_m; ii++){
        reduce(newPartNum[ii], newPartNum[ii], OpAddAssign());
        nBin_m[ii] = newPartNum[ii];
        INFOMSG( "After reset Bin: "<< ii <<", particle(s): "<<newPartNum[ii]<<endl);
    }
}

PartBins::~PartBins()
{
    if (nBin_m) {
        delete nBin_m;
        delete xbinmax_m;
        delete xbinmin_m;
        delete binsEmitted_m;
    }
    tmppart_m.clear();
    isEmitted_m.clear();
}


bool PartBins::getPart(size_t n, int bin, vector<double> &p) {

    if (tmppart_m[n][6] == bin) {
        p = tmppart_m[n];
        return true;
    }
    else
        return false;
}

/** /brief There is only a local sort, no global yet */
void PartBins::sortArray() {
    extern Inform* gmsg2all;

    /** sort the vector of particles such that position of the particles decrease with increasing index.
        Then push the particles back by 1e-13 s * beta * c (approximately one step).
        In order that the method getBin(double x) works xmin_m has to be lowered a bit more.
    */

    double sshift = sqrt(1. - (1. / (gamma_m*gamma_m))) * Physics::c * 1e-13;
    std::sort(tmppart_m.begin(), tmppart_m.end(), DescendingLocationSort(2));
    xmax_m = tmppart_m[0][2];
    xmin_m = tmppart_m.back()[2];

    for (int n=0; n< tmppart_m.size(); n++)
        tmppart_m[n][2] -= xmax_m + sshift; /* push particles back */


    /**
       Because on CRAY we crash
       if writing like this

       #ifndef __LIBCATAMOUNT__

       ofstream dist("gun_distribution.txt");
       for (int n=0; n< tmppart_m.size(); n++)
       {
       dist.precision(8);
       dist << tmppart_m[n][0] << "\t"
       << tmppart_m[n][1] << "\t"
       << tmppart_m[n][2] << "\t"
       << tmppart_m[n][3] << "\t"
       << tmppart_m[n][4] << "\t";
       dist.precision(15);
       dist << tmppart_m[n][5] << endl;
       }

       dist.close();

       #else
       *gmsg << "On Cray can not dump gun_distribution.txt" << endl;
       #endif
    */
    xmin_m -= xmax_m + 0.0001*(xmax_m - xmin_m) + sshift; /* lower the limits */
    xmax_m = -sshift;

    reduce(xmin_m,xmin_m,OpMinAssign());
    reduce(xmax_m,xmax_m,OpMaxAssign());

    hBin_m = (fabs(xmax_m - xmin_m))/(bins_m);
    calcHBins();
    for (int n=0; n<bins_m; n++)
        if (nBin_m[n] == 0) setBinEmitted(n);
}

void PartBins::calcHBins() {

    for (int n=0; n < tmppart_m.size(); n++)
        tmppart_m[n][6] = getBin(tmppart_m[n][2]);
    calcExtrema();
}

size_t PartBins::getSum() {
    size_t s = 0;
    for (int n=0; n< bins_m; n++)
        s += nBin_m[n];
    return s;
}

void PartBins::calcGlobalExtrema(){
    xmin_m = DBL_MAX;
    xmax_m = -DBL_MAX;
    for (int n=0; n< tmppart_m.size(); n++) {
        if(tmppart_m[n][2]<=xmin_m)
            xmin_m = tmppart_m[n][2];
        if(tmppart_m[n][2]>=xmax_m)
            xmax_m = tmppart_m[n][2];
    }
    double xdiff = 0.01*(xmax_m - xmin_m);
    xmin_m -= xdiff;
    xmax_m += xdiff;
}

void PartBins::calcExtrema() {
    for (int n=0; n< tmppart_m.size(); n++) {
        if (xbinmin_m[(int)tmppart_m[n][6]] >= tmppart_m[n][2])
            xbinmin_m[(int)tmppart_m[n][6]] = tmppart_m[n][2];

        if (xbinmax_m[(int)tmppart_m[n][6]] <= tmppart_m[n][2])
            xbinmax_m[(int)tmppart_m[n][6]] = tmppart_m[n][2];
    }
}

Inform &PartBins::print(Inform &os) {

    os <<"-----------------------------------------"<<endl;
    os <<"     CREATE BINNED GAUSS DISTRIBUTION DONE        "<<endl;

    os << "Bins= " << bins_m << " hBin= " << hBin_m << endl;
    os << "minVal= " << xmin_m << " [m] maxVal= " << xmax_m << " [m] Particle vector length " << tmppart_m.size() << endl;

    for (int i=0; i<bins_m; i++)
        os << "Bin #" << i << " bin value " << getTotalNumPerBin(i) << "\t minx= " << xbinmin_m[i]  << "\t maxx= " << xbinmax_m[i] << endl;

    if (getLastemittedBin() >= 0)
        os << "Last emitted bin is " << getLastemittedBin() << endl;
    else
        os << "No bin is emitted !" << endl;
    return os;
}

int PartBins::getBin(double x) {
    /**
       returns the index of the bin to which the particle with z = 'x' belongs.
       If getBin returns b < 0 || b >= bins_m, then is x out of range!
    */
    int b = (int) floor(fabs(xmax_m - x)/hBin_m);
    nBin_m[b]++;
    return b;
}
