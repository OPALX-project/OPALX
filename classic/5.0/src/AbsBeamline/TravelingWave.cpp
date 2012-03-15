// ------------------------------------------------------------------------
// $RCSfile: TravelingWave.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: TravelingWave
//   Defines the abstract interface for an accelerating structure.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/TravelingWave.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "Fields/Fieldmap.hh"
#include "Physics/Physics.h"
#include <iostream>
#include <fstream>

extern Inform *gmsg;

// Class TravelingWave
// ------------------------------------------------------------------------

TravelingWave::TravelingWave():
    Component(),
    NumCells_m(0),
    fast_m(false)
{}


TravelingWave::TravelingWave(const TravelingWave &right):
    Component(right),
    CoreFilename_m(right.CoreFilename_m),
    //   EntryFilename_m(right.EntryFilename_m),
    //   ExitFilename_m(right.ExitFilename_m),
    CoreFieldmap_m(right.CoreFieldmap_m),
    //   EntryFringeField_m(right.EntryFringeField_m),
    //   ExitFringeField_m(right.ExitFringeField_m),
    scale_m(right.scale_m),
    scaleCore_m(right.scaleCore_m),
    frequency_m(right.frequency_m),
    phase_m(right.phase_m),
    phaseCore1_m(right.phaseCore1_m),
    phaseCore2_m(right.phaseCore2_m),
    phaseExit_m(right.phaseExit_m),
    startField_m(right.startField_m),
    startCoreField_m(right.startCoreField_m),
    startExitField_m(right.startExitField_m),
    mappedStartExitField_m(right.mappedStartExitField_m),
    PeriodLength_m(right.PeriodLength_m),
    NumCells_m(right.NumCells_m),
    CellLength_m(right.CellLength_m),
    fast_m(right.fast_m),
    Mode_m(right.Mode_m)
{}


TravelingWave::TravelingWave(const string &name):
    Component(name)
{}


TravelingWave::~TravelingWave() {
    Fieldmap::deleteFieldmap(CoreFilename_m);
    //   Fieldmap::deleteFieldmap(EntryFilename_m);
    //   Fieldmap::deleteFieldmap(ExitFilename_m);
}


void TravelingWave::accept(BeamlineVisitor &visitor) const {
    visitor.visitTravelingWave(*this);
}

void TravelingWave::setFieldMapFN(string fn) {
    CoreFilename_m = fn;
}

// void TravelingWave::setEntryFieldMapFN(string fn)
// {
//   EntryFilename_m = fn;
// }

// void TravelingWave::setExitFieldMapFN(string fn)
// {
//   ExitFilename_m = fn;
// }

string TravelingWave::getFieldMapFN() const {
    return CoreFilename_m;
}

void TravelingWave::setAmplitudem(double vPeak) {
    scale_m = vPeak;
}

void TravelingWave::setFrequencym(double freq) {
    frequency_m = freq;
}

void TravelingWave::setPhasem(double phase) {
    using Physics::pi;

    phase_m = phase;
    phaseCore1_m = phase_m + pi * Mode_m / 2.0;
    phaseCore2_m = phase_m + pi * Mode_m * 1.5;
    phaseExit_m = phase_m + 2.0 * pi * (1.0 - NumCells_m * Mode_m);
}

double TravelingWave::getPhasem() const {
    return phase_m;
}

void TravelingWave::setNumCells(int NumCells) {
    NumCells_m = NumCells;
}

void TravelingWave::setFast(bool fast) {
    fast_m = fast;
}


bool TravelingWave::getFast() const {
    return fast_m;
}

/**
 * ENVELOPE COMPONENT for radial focussing of the beam
 * Calculates the transverse envelope component for the solenoid
 * element and adds it to the K vector
*/
void TravelingWave::addKR(int i, double t, Vector_t &K) {
    Inform msg("TravelingWave::addK()");

    // get field components:
    Vector_t tmpE(0.0, 0.0, 0.0);
    Vector_t tmpB(0.0, 0.0, 0.0);
    Vector_t tmpE_diff(0.0, 0.0, 0.0);
    Vector_t tmpB_diff(0.0, 0.0, 0.0);
    double px = RefPartBunch_m->getX(i) - dx_m;
    double py = RefPartBunch_m->getY(i) - dy_m;
    double pz = RefPartBunch_m->getZ(i) - startField_m - ds_m;
    const Vector_t tmpA0(px, py, pz);

    bool out_of_bounds;
    DiffDirection zdir(DZ);
    out_of_bounds = CoreFieldmap_m->getFieldstrength_fdiff(tmpA0, tmpE_diff, tmpB_diff, zdir);
    out_of_bounds = CoreFieldmap_m->getFieldstrength(tmpA0, tmpE, tmpB);

    double Ez = scale_m * tmpE(2);

    // calculate parameters
    double g = RefPartBunch_m->getGamma(i);
    double b = RefPartBunch_m->getBeta(i);
    double freq = CoreFieldmap_m->getFrequency();

    double phase_temp = getPhase();

    //BET: pi0x := Pi / x
    double pi02 = Physics::pi / 2;
    //BET, I used the phase "phase_temp", which is in rad
    double wtf = freq * t + phase_temp - pi02;
    double dfi = Physics::pi * phase_m / 180.0;
    double Emax = 1;

    //FIXME: the following equations were copied from the BET programme and NOT TESTED yet!
    //BET: "length of SW pattern in field mapping file"
    double zInc = PeriodLength_m;
    //FIXME: Andreas fragen!!!
    double d    = dfi * Physics::c / (2 * Physics::pi * freq);
    double kd   = freq * d / Physics::c;
    double L    = 2.0 * zInc + d * (NumCells_m - 1);
    double d3   = 3.0 * d;

    double z1 = 0.0, z2 = 0.0, k = 0.0;
    if(pz < zInc)
        k = Emax * (tmpE_diff(2) * sin(wtf) + b * freq * Ez * cos(wtf) / Physics::c);
    else if((L - pz) < zInc) {  // half end-cells
        //from BET
        z2 = pz - L + 2.0 * zInc + d3;
        const Vector_t tmpA1(px, py, z2);
        tmpE = (0.0, 0.0, 0.0);
        tmpB = (0.0, 0.0, 0.0);
        tmpE_diff = (0.0, 0.0, 0.0);
        tmpB_diff = (0.0, 0.0, 0.0);
        out_of_bounds = CoreFieldmap_m->getFieldstrength(tmpA1, tmpE, tmpB);
        out_of_bounds = CoreFieldmap_m->getFieldstrength_fdiff(tmpA1, tmpE_diff, tmpB_diff, zdir);

        wtf -= ((4 - NumCells_m) * Physics::pi / 3 + pi02);
        if((NumCells_m % 2) == 1) wtf += Physics::pi;

        k = Emax * (tmpE_diff(2) * cos(wtf) - b * freq * tmpE(2) * sin(wtf) / Physics::c);
    } else {
        wtf += Physics::pi / 6;
        z1 = zInc + fmod(pz - zInc, d3);
        z2 = zInc + fmod(pz + d - zInc, d3);

        // define new vectors for new z coordinates,
        //FIXME (ff): I guess that could be done easier and faster!
        const Vector_t tmpA1(px, py, z1);
        const Vector_t tmpA2(px, py, z2);

        Vector_t tmpE0(0.0, 0.0, 0.0);
        Vector_t tmpB0(0.0, 0.0, 0.0);
        Vector_t tmpE_diff0(0.0, 0.0, 0.0);
        Vector_t tmpB_diff0(0.0, 0.0, 0.0);

        Vector_t tmpE1(0.0, 0.0, 0.0);
        Vector_t tmpB1(0.0, 0.0, 0.0);
        Vector_t tmpE_diff1(0.0, 0.0, 0.0);
        Vector_t tmpB_diff1(0.0, 0.0, 0.0);

        out_of_bounds = CoreFieldmap_m->getFieldstrength(tmpA1, tmpE0, tmpB0);
        out_of_bounds = CoreFieldmap_m->getFieldstrength(tmpA2, tmpE1, tmpB1);
        out_of_bounds = CoreFieldmap_m->getFieldstrength_fdiff(tmpA1, tmpE_diff0, tmpB_diff0, zdir);
        out_of_bounds = CoreFieldmap_m->getFieldstrength_fdiff(tmpA2, tmpE_diff1, tmpB_diff1, zdir);
        k = Emax * (tmpE_diff0(2) * sin(wtf) - tmpE_diff1(2) * sin(wtf - kd) +
                    b * freq * (tmpE0(2) * cos(wtf) - tmpE_diff1(2) * cos(wtf - kd)) / Physics::c);
    }

    k *= (Physics::q_e / (2.0 * g * Physics::EMASS));
    /*
    if (Ndamp > 0) {
        Kloc *= fDamp(dz);
    }
    */
    K += Vector_t(k, k, 0.0);
}

/*
   BET code for radial K component

   //1) Code für die Berechnung von d3:
   dfi = pi*dfi_deg/180.0;
   d   = dfi*c/(pi2*f0);
   kd  = w0*d/c;
   L   = 2.0*zInc + d*(N-1);
   d3  = 3.0*d;

   //2) Code für die Berechung von K:

   if (dz < zInc) {
       k   =
           Emax*(dEdz->get(dz)*sin(wtf) + b*w0*Ez->get(dz)*cos(wtf)/c);
   } else if ((L - dz) < zInc) { // half end-cells
       z2  = dz - L + 2.0*zInc + d3;

       wtf -= ((4.0 - N)*pi/3.0 + pi02);
       if ((N % 2) == 1) wtf += pi;

       k = Emax*(dEdz->get(z2)*cos(wtf) - b*w0*Ez->get(z2)*sin(wtf)/c);
   } else {
       wtf += pi06;

       z1 = zInc + fmod(dz-zInc,d3);
       z2 = zInc + fmod(dz + d - zInc,d3);

       k = Emax*(dEdz->get(z1)*sin(wtf) - dEdz->get(z2)*sin(wtf - kd) +
                 b*w0*(Ez->get(z1)*cos(wtf) - Ez->get(z2)*cos(wtf - kd))/c);

   }
   k *= (e0/(2.0*g*m));
   Kloc  = Field(Cxy*k,Cxy*k,0.0);
   kLast = k;
   if (Ndamp > 0) {
     Kloc *= fDamp(dz);
   }
}
return Kloc;
*/

void TravelingWave::addKT(int i, double t, Vector_t &K) {
    // get K from radial function, solves the lastk problem from BET
    Vector_t tempK(0.0, 0.0, 0.0);
    addKR(i, t, tempK);

    // x and y component are identical and equal to k
    double oldK = tempK(0);

    //get coordinates
    double temp1 = RefPartBunch_m->getX(i) - dx_m;
    double temp2 = RefPartBunch_m->getY(i) - dy_m;

    double dx = RefPartBunch_m->getX0(i) - dx_m;
    double dy = RefPartBunch_m->getY0(i) - dy_m;

    K += Vector_t(oldK * dx, oldK * dy, 0.0);
}

/*
    Function to obtain z-k-component, copied from BET:
   Field Kloc;

   if (isActive(z)){
       double k,kx,ky;

       if (kLast == MAXFLOAT) {
           Kloc = getK(g,b,z,t,m);  // dummy call to set kLast
       }
       k   = kLast;
       kx  = 0.0; ky = 0.0;
       if (Cxy < 1.0) {
           Field Eloc = getE(z,t);
           double
               cf = -e0/(g*m);

           kx = cf*Eloc.x;
           ky = cf*Eloc.y;
       }
       Kloc = Field(k*Cxy*(x0-x) + kx,
               k*Cxy*(y0-y) + ky,
               0.0);
       if (Ndamp > 0) {
           Kloc *= fDamp(z - z0);
       }
   }

   return Kloc;
*/


bool TravelingWave::apply(const int &i, const double &t, double E[], double B[]) {
    Vector_t Ev(0, 0, 0), Bv(0, 0, 0);
    Vector_t Rt(RefPartBunch_m->getX(i), RefPartBunch_m->getY(i), RefPartBunch_m->getZ(i));
    if(apply(Rt, Vector_t(0.0), t, Ev, Bv)) return true;

    E[0] = Ev(0);
    E[1] = Ev(1);
    E[2] = Ev(2);
    B[0] = Bv(0);
    B[1] = Bv(1);
    B[2] = Bv(2);

    return false;
}

bool TravelingWave::apply(const int &i, const double &t, Vector_t &E, Vector_t &B) {
    double tmpcos, tmpsin;

    Vector_t tmpR(RefPartBunch_m->getX(i) - dx_m, RefPartBunch_m->getY(i) - dy_m , RefPartBunch_m->getZ(i) - startField_m - ds_m);
    Vector_t tmpE(0.0, 0.0, 0.0), tmpB(0.0, 0.0, 0.0);
    bool out_of_bounds = false;


    if(tmpR(2) < startCoreField_m) {
        tmpcos =  scale_m * cos(frequency_m * t + phase_m);
        tmpsin = -scale_m * sin(frequency_m * t + phase_m);

    } else if(tmpR(2) < startExitField_m) {
        Vector_t tmpE2(0.0, 0.0, 0.0), tmpB2(0.0, 0.0, 0.0);
        tmpR(2) -= startCoreField_m;
        const double z = tmpR(2);
        tmpR(2) = tmpR(2) - PeriodLength_m * floor(tmpR(2) / PeriodLength_m);
        tmpR(2) += startCoreField_m;

        tmpcos =  scaleCore_m * cos(frequency_m * t + phaseCore1_m);
        tmpsin = -scaleCore_m * sin(frequency_m * t + phaseCore1_m);
        out_of_bounds = CoreFieldmap_m->getFieldstrength(tmpR, tmpE, tmpB);
        E += tmpcos * tmpE;
        B += tmpsin * tmpB;

        tmpE = 0.0;
        tmpB = 0.0;

        tmpR(2) = z + CellLength_m;
        tmpR(2) = tmpR(2) - PeriodLength_m * floor(tmpR(2) / PeriodLength_m);
        tmpR(2) += startCoreField_m;

        tmpcos =  scaleCore_m * cos(frequency_m * t + phaseCore2_m);
        tmpsin = -scaleCore_m * sin(frequency_m * t + phaseCore2_m);

    } else {
        tmpcos =  scale_m * cos(frequency_m * t + phaseExit_m);
        tmpsin = -scale_m * sin(frequency_m * t + phaseExit_m);
        tmpR(2) -= mappedStartExitField_m;

    }

    out_of_bounds = out_of_bounds || CoreFieldmap_m->getFieldstrength(tmpR, tmpE, tmpB);
    E += tmpcos * tmpE;
    B += tmpsin * tmpB;

    return out_of_bounds;
}

bool TravelingWave::apply(const Vector_t &R, const Vector_t &centroid, const double &t, Vector_t &E, Vector_t &B) {
    double tmpcos, tmpsin;
    Vector_t tmpR(R(0) - dx_m, R(1) - dy_m , R(2) - startField_m - ds_m);
    Vector_t tmpE(0.0, 0.0, 0.0), tmpB(0.0, 0.0, 0.0);
    bool out_of_bounds = false;


    if(tmpR(2) < startCoreField_m) {
        tmpcos =  scale_m * cos(frequency_m * t + phase_m);
        tmpsin = -scale_m * sin(frequency_m * t + phase_m);

    } else if(tmpR(2) < startExitField_m) {
        Vector_t tmpE2(0.0, 0.0, 0.0), tmpB2(0.0, 0.0, 0.0);
        tmpR(2) -= startCoreField_m;
        const double z = tmpR(2);
        tmpR(2) = tmpR(2) - PeriodLength_m * floor(tmpR(2) / PeriodLength_m);
        tmpR(2) += startCoreField_m;

        tmpcos =  scaleCore_m * cos(frequency_m * t + phaseCore1_m);
        tmpsin = -scaleCore_m * sin(frequency_m * t + phaseCore1_m);
        out_of_bounds = CoreFieldmap_m->getFieldstrength(tmpR, tmpE, tmpB);
        E += tmpcos * tmpE;
        B += tmpsin * tmpB;

        tmpE = 0.0;
        tmpB = 0.0;

        tmpR(2) = z + CellLength_m;
        tmpR(2) = tmpR(2) - PeriodLength_m * floor(tmpR(2) / PeriodLength_m);
        tmpR(2) += startCoreField_m;

        tmpcos =  scaleCore_m * cos(frequency_m * t + phaseCore2_m);
        tmpsin = -scaleCore_m * sin(frequency_m * t + phaseCore2_m);

    } else {
        tmpcos =  scale_m * cos(frequency_m * t + phaseExit_m);
        tmpsin = -scale_m * sin(frequency_m * t + phaseExit_m);
        tmpR(2) -= mappedStartExitField_m;

    }

    out_of_bounds = out_of_bounds || CoreFieldmap_m->getFieldstrength(tmpR, tmpE, tmpB);
    E += tmpcos * tmpE;
    B += tmpsin * tmpB;

    return out_of_bounds;

}

void TravelingWave::initialise(PartBunch *bunch, double &startField, double &endField, const double &scaleFactor) {
    using Physics::pi;
    using Physics::two_pi;

    Inform msg("TravelingWave ");
    stringstream errormsg;
    double zbegin = 0.0;
    double zend = 0.0;
    double tmpDouble = 0.0;

    RefPartBunch_m = bunch;

    CoreFieldmap_m = Fieldmap::getFieldmap(CoreFilename_m, fast_m);
    if(CoreFieldmap_m != NULL) {
        double zBegin = 0.0, zEnd = 0.0, rBegin = 0.0, rEnd = 0.0;
        CoreFieldmap_m->getFieldDimensions(zBegin, zEnd, rBegin, rEnd);

        if(zEnd > zBegin) {
            msg << getName() << " using file ";
            CoreFieldmap_m->getInfo(&msg);
            if(fabs((frequency_m - CoreFieldmap_m->getFrequency()) / frequency_m) > 0.01) {
                errormsg << "FREQUENCY IN INPUT FILE DIFFERENT THAN IN FIELD MAP '" <<  CoreFilename_m + "';\n"
                         << frequency_m / two_pi * 1e-6 << " MHz <> "
                         << CoreFieldmap_m->getFrequency() / two_pi * 1e-6 << " MHz; TAKE ON THE LATTER\n";
                string errormsg_str = Fieldmap::typeset_msg(errormsg.str(), "warning");
                msg << errormsg_str << "\n"
                    << endl;
                if(Ippl::myNode() == 0) {
                    ofstream omsg("errormsg.txt", ios_base::app);
                    omsg << errormsg_str << endl;
                    omsg.close();
                }
                frequency_m = CoreFieldmap_m->getFrequency();
            }

            if(dx_m > 1e-10 || dy_m > 1e-10 || ds_m > 1e-10)
                msg << "misaligned by dx = " << dx_m << ", dy = " << dy_m << ", dz = " << ds_m << endl;

            if(hasAttribute("MODE")) {
                Mode_m = getAttribute("MODE");
            } else {
                Mode_m = 1. / 3.;
                errormsg.str("");
                errormsg  << "NO MODE GIVEN; 2\\pi/3 MODE ASSUMED.";
                string errormsg_str = Fieldmap::typeset_msg(errormsg.str(), "warning");
                msg << errormsg_str << "\n"
                    << endl;
                if(Ippl::myNode() == 0) {
                    ofstream omsg("errormsg.txt", ios_base::app);
                    omsg << errormsg_str << endl;
                    omsg.close();
                }
            }

            PeriodLength_m = (zEnd - zBegin) / 2.0;
            CellLength_m = PeriodLength_m * Mode_m;

            endField = startField + CellLength_m * NumCells_m + PeriodLength_m / 2.0;
            startField -= PeriodLength_m / 2.0;

            startField_m = startField;
            startCoreField_m = PeriodLength_m / 2.0;
            startExitField_m = startCoreField_m + NumCells_m * CellLength_m;
            mappedStartExitField_m = startExitField_m - 3.0 * PeriodLength_m / 2.0;

            scaleCore_m = scale_m / sin(2.0 * pi * Mode_m);
            phaseCore1_m = phase_m + pi * Mode_m / 2.0;
            phaseCore2_m = phase_m + pi * Mode_m * 1.5;
            phaseExit_m = phase_m + 2.0 * pi * (1.0 - NumCells_m * Mode_m);
        } else {
            endField = startField - 1e-3;
        }
    } else {
        endField = startField - 1e-3;
    }
}

void TravelingWave::finalise()
{}

bool TravelingWave::bends() const {
    return false;
}


void TravelingWave::goOnline() {
    Fieldmap::readMap(CoreFilename_m);
    online_m = true;
}

void TravelingWave::goOffline() {
    Fieldmap::freeMap(CoreFilename_m);
}

void TravelingWave::getDimensions(double &zBegin, double &zEnd) const {
    zBegin = startField_m;
    zEnd = startField_m + NumCells_m * CellLength_m;
}


const string &TravelingWave::getType() const {
    static const string type("TravelingWave");
    return type;
}

