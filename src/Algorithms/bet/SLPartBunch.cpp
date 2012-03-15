#ifdef HAVE_ENVELOPE_SOLVER
// ------------------------------------------------------------------------
// $RCSfile: SLPartBunch.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class SLPartBunch
//   Interface to a particle bunch.
//   Can be used to avoid use of a template in user code.
//
// ------------------------------------------------------------------------
// Class category: Algorithms
// ------------------------------------------------------------------------
//
// $Date: 2008/07/19 18:57:53 $
// $Author: adelmann,hantschel $
//
// ------------------------------------------------------------------------

#include "Algorithms/bet/SLPartBunch.h"
#include <iostream>
#include <cfloat>
#include <fstream>
#include "Physics/Physics.h"
#include "AbstractObjects/OpalData.h"

using Physics::c;
using Physics::pi;

extern Inform* gmsg;

// Class SLPartBunch
// ------------------------------------------------------------------------

SLPartBunch::SLPartBunch(const PartData *ref):
    reference(ref),
    PartBunch(ref)
{
}


SLPartBunch::SLPartBunch(const SLPartBunch &rhs):
    reference(rhs.reference),
    PartBunch(rhs.reference)
{}


SLPartBunch::SLPartBunch(const std::vector<Particle> &rhs, const PartData *ref):
    reference(ref),
    PartBunch(ref)
{}


SLPartBunch::~SLPartBunch()
{

}

double SLPartBunch::calcGamma(double b) 
{
    return sqrt(1/(1+b*b));
}

/** 
 * main function for initialization of a BET bunch, 
 * obtained from BET file system.C 
 */
void SLPartBunch::iniBetBunch(int sli, double Q, double energy, double width, double frac, double current, double center, double bX, double bY, double mX, double mY, double Bz0)
{  
    EnvelopeBunch *tempB = new EnvelopeBunch(sli); 

    // Set charge Q and energy E0 
    tempB->setCharge(Q);
    tempB->setEnergy(energy);

    // Set slope parameter 
    wfraction=1; 

    /** Set the longitudinal and transversal shape of the bunch. 
     *  This determines the shape of all slices. 
     */

    // Shape of the bunch should be rectangular
    int i2 = 0; 

    /* set emmittance

       tempB->setEy(0);
       tempB->setEx(0);
    */

    tempB->setLShape(i2?bsGauss:bsRect,center,width,frac);
    tempB->setTShape(mX,mY,bX,bY,Bz0);

    *gmsg << "Bet bunch initialized " << endl << "Important values:" << endl;
    *gmsg << "Charge: Q=" << Q << "  No. of Slices: N=" << sli << "  Energy E0=" << energy << endl << endl;  

    // Set solver method 
    tempB->setSolver(12); // Solver = 12 in BET -> default

    BetBunch_m=tempB; 
}

int SLPartBunch::getN()
{
    return BetBunch_m->getN();
}

void SLPartBunch::setZ(int i, double zcoo)
{
    BetBunch_m->setZ(i,zcoo);
}

double SLPartBunch::getZ(int i) 
{
    return BetBunch_m->getZ(i);
}

double SLPartBunch::getX(int i) 
{
    return BetBunch_m->getX(i);
}

double SLPartBunch::getY(int i) 
{
    return BetBunch_m->getY(i);
}

double SLPartBunch::getX0(int i) 
{
    return BetBunch_m->getX0(i);
}

double SLPartBunch::getY0(int i) 
{
    return BetBunch_m->getY0(i);
}

double SLPartBunch::getPx(int i) 
{
    return BetBunch_m->getPx(i);
}

double SLPartBunch::getPy(int i) 
{
    return BetBunch_m->getPy(i);
}

double SLPartBunch::getPz(int i) 
{
    return BetBunch_m->getPz(i);
}

double SLPartBunch::getPx0(int i) 
{
    return BetBunch_m->getPx0(i);
}

double SLPartBunch::getPy0(int i) 
{
    return BetBunch_m->getPy0(i);
}

double SLPartBunch::getBeta(int i) 
{
    return BetBunch_m->getBeta(i);
}

double SLPartBunch::getGamma(int i) {
    return BetBunch_m->getGamma(i);
}

// functions to set value for K (needed for envelope equation)
void SLPartBunch::setKR(Vector_t value, int i)
{
    BetBunch_m->KR[i] = value;
}

void SLPartBunch::setKT(Vector_t value, int i)
{
    BetBunch_m->KT[i] = value;
}

// functions to set EF and BF 
void SLPartBunch::plotR()
{
    BetBunch_m->plotR();
}

void SLPartBunch::setEF(Vector_t value, int i)
{
    BetBunch_m->EF[i] = value;
}

void SLPartBunch::setBF(Vector_t value, int i)
{
    BetBunch_m->BF[i] = value;
}

// Function for one timestep, calls bunch->run() 
void SLPartBunch::tstep(double cat,double step, size_t nstep)
{
    Inform msg("tstep");

    // make sure, I-profile and spacecharge (SC) is calculated
    BetBunch_m->updateFields();  

    // timestep from BET
    BetBunch_m->run(cat,step);    

    // sets time for SLPartBunch
    setT(BetBunch_m->getT());

    // calculates new profile
    BetBunch_m->updateFields();   

    //msg << "zcoor: " << BetBunch_m->zAvg() << "\time SL: " << getT() << "\ttime BET: " << BetBunch_m->getT() << endl;

    msg << " Step " << nstep << " at " << BetBunch_m->zAvg() << " [m] t= " << BetBunch_m->getT() << " [s] E=" << BetBunch_m->Eavg()*1e-6 << " [MeV] " << endl;
}

// compare BET and OPAL time
void SLPartBunch::actT()
{
    setT(BetBunch_m->getT());
}

/** 
 * Complete set of output functions: Opens and closes file and writes output 
 *
 * 1st: writes output for slices 
*/
void SLPartBunch::writeSlicedBetBunch(const char *fName)
{
    // Define outputfile... different parameter as for Bunch->write()
    FILE *fSlice  = NULL;
    fSlice = fopen(fName,"w");
    BetBunch_m->writeSlice(fSlice,oFormat_sdds);
    fclose(fSlice);
}

// 2nd: Small output function, just uses the op-routine Bunch->write() 
void SLPartBunch::writeBetBunch(const char *fname)
{
    *gmsg << "Output for BetBunch:" << endl;
    BetBunch_m->write(fname);
    *gmsg << endl << "Finished BetBunch-Output" << endl; 
}

// 3rd: writeStatistics 
void SLPartBunch::writeBetStat(const char *fname)
{
    FILE *fStat = NULL;
    fStat = fopen(fname,"w");
    //*gmsg << "Opened STAT file" << endl;
    BetBunch_m->writeStats(fStat,wfraction,oFormat_sdds);
    fclose(fStat);
    //*gmsg << "Closed STAT file" << endl;
}

// 4th: MAIN output function of BET 
void SLPartBunch::BetOut(FILE* dat, FILE* sli)
{
    BetBunch_m->writeStats(dat,wfraction,oFormat_sdds); 
    BetBunch_m->writeSlice(sli,oFormat_sdds);
}

Inform &SLPartBunch::slprint(Inform &os)
{
    os << "* ************** S L B U N C H ********************************************************* " << endl;

    os << "* ********************************************************************************** " << endl;
    return os;
}

#endif
