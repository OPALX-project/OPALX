// ------------------------------------------------------------------------
// $RCSfile: RFCavity.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: RFCavity
//   Defines the abstract interface for an accelerating structure.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: andreas adelmann $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/RFCavity.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "Algorithms/PartBunch.h"
#include "Fields/Fieldmap.h"
#include "Utilities/GeneralClassicException.h"

#include "gsl/gsl_interp.h"
#include "gsl/gsl_spline.h"
#include <iostream>
#include <fstream>

#ifdef OPAL_NOCPLUSPLUS11_NULLPTR
#define nullptr NULL
#endif

extern Inform *gmsg;

using namespace std;

// Class RFCavity
// ------------------------------------------------------------------------

RFCavity::RFCavity():
    Component(),
    filename_m(""),
    scale_m(1.0),
    phase_m(0.0),
    frequency_m(0.0),
    ElementEdge_m(0.0),
    startField_m(0.0),
    endField_m(0.0),
    type_m(SW),
    fast_m(false),
    autophaseVeto_m(false),
    rmin_m(0.0),
    rmax_m(0.0),
    angle_m(0.0),
    sinAngle_m(0.0),
    cosAngle_m(0.0),
    pdis_m(0.0),
    gapwidth_m(0.0),
    phi0_m(0.0),
    RNormal_m(nullptr),
    VrNormal_m(nullptr),
    DvDr_m(nullptr),
    num_points_m(0),
    phase_td_m(nullptr),
    amplitude_td_m(nullptr),
    frequency_td_m(nullptr),
    phase_name_m(""),
    amplitude_name_m(""),
    frequency_name_m("")
{
    setElType(isRF);
}


RFCavity::RFCavity(const RFCavity &right):
    Component(right),
    filename_m(right.filename_m),
    scale_m(right.scale_m),
    phase_m(right.phase_m),
    frequency_m(right.frequency_m),
    ElementEdge_m(right.ElementEdge_m),
    startField_m(right.startField_m),
    endField_m(right.endField_m),
    type_m(right.type_m),
    fast_m(right.fast_m),
    autophaseVeto_m(right.autophaseVeto_m),
    rmin_m(right.rmin_m),
    rmax_m(right.rmax_m),
    angle_m(right.angle_m),
    sinAngle_m(right.sinAngle_m),
    cosAngle_m(right.cosAngle_m),
    pdis_m(right.pdis_m),
    gapwidth_m(right.gapwidth_m),
    phi0_m(right.phi0_m),
    RNormal_m(nullptr),
    VrNormal_m(nullptr),
    DvDr_m(nullptr),
    num_points_m(right.num_points_m),
    phase_td_m(right.phase_td_m),
    amplitude_td_m(right.amplitude_td_m),
    frequency_td_m(right.frequency_td_m),
    phase_name_m(right.phase_name_m),
    amplitude_name_m(right.amplitude_name_m),
    frequency_name_m(right.frequency_name_m) {

    setElType(isRF);

    std::vector<string>::const_iterator fname_it;
    for(fname_it = right.multiFilenames_m.begin(); fname_it != right.multiFilenames_m.end(); ++ fname_it) {
        multiFilenames_m.push_back(*fname_it);
    }
    std::vector<Fieldmap *>::const_iterator fmap_it;
    for(fmap_it = right.multiFieldmaps_m.begin(); fmap_it != right.multiFieldmaps_m.end(); ++ fmap_it) {
        multiFieldmaps_m.push_back(*fmap_it);
    }
    std::vector<double>::const_iterator scale_it;
    for(scale_it = right.multiScales_m.begin(); scale_it != right.multiScales_m.end(); ++ scale_it) {
        multiScales_m.push_back(*scale_it);
    }
    std::vector<double>::const_iterator phase_it;
    for(phase_it = right.multiPhases_m.begin(); phase_it != right.multiPhases_m.end(); ++ phase_it) {
        multiPhases_m.push_back(*phase_it);
    }
    std::vector<double>::const_iterator freq_it;
    for(freq_it = right.multiFrequencies_m.begin(); freq_it != right.multiFrequencies_m.end(); ++ freq_it) {
        multiFrequencies_m.push_back(*freq_it);
    }
}


RFCavity::RFCavity(const std::string &name):
    Component(name),
    filename_m(""),
    scale_m(1.0),
    phase_m(0.0),
    frequency_m(0.0),
    ElementEdge_m(0.0),
    startField_m(0.0),
    endField_m(0.0),
    type_m(SW),
    fast_m(false),
    autophaseVeto_m(false),
    rmin_m(0.0),
    rmax_m(0.0),
    angle_m(0.0),
    sinAngle_m(0.0),
    cosAngle_m(0.0),
    pdis_m(0.0),
    gapwidth_m(0.0),
    phi0_m(0.0),
    RNormal_m(nullptr),
    VrNormal_m(nullptr),
    DvDr_m(nullptr),
    phase_td_m(nullptr),
    amplitude_td_m(nullptr),
    frequency_td_m(nullptr),
    //     RNormal_m(std::nullptr_t(NULL)),
    //     VrNormal_m(std::nullptr_t(NULL)),
    //     DvDr_m(std::nullptr_t(NULL)),
    num_points_m(0) {
    setElType(isRF);
}


RFCavity::~RFCavity() {
    // FIXME: in deleteFielmak, a map find makes problems
    //       Fieldmap::deleteFieldmap(filename_m);
    //~ if(RNormal_m) {
    //~ delete[] RNormal_m;
    //~ delete[] VrNormal_m;
    //~ delete[] DvDr_m;
    //~ }
}


std::shared_ptr<AbstractTimeDependence> RFCavity::getAmplitudeModel() const {
  return amplitude_td_m;
}

std::shared_ptr<AbstractTimeDependence> RFCavity::getPhaseModel() const {
  return phase_td_m;
}

std::shared_ptr<AbstractTimeDependence> RFCavity::getFrequencyModel() const {
  return frequency_td_m;
}

void RFCavity::setAmplitudeModel(std::shared_ptr<AbstractTimeDependence> amplitude_td) {
  amplitude_td_m = amplitude_td;
}

void RFCavity::setPhaseModel(std::shared_ptr<AbstractTimeDependence> phase_td) {
  phase_td_m = phase_td;
}

void RFCavity::setFrequencyModel(std::shared_ptr<AbstractTimeDependence> frequency_td) {
  INFOMSG("frequency_td " << frequency_td << endl);
  frequency_td_m = frequency_td;
}


void RFCavity::setFrequencyModelName(std::string name) {
  frequency_name_m=name;
}




void RFCavity::accept(BeamlineVisitor &visitor) const {
    visitor.visitRFCavity(*this);
}

void RFCavity::dropFieldmaps() {
    std::vector<Fieldmap *>::iterator fmap_it;
    for(fmap_it = multiFieldmaps_m.begin(); fmap_it != multiFieldmaps_m.end(); ++ fmap_it) {
        *fmap_it = NULL;
    }
    multiFilenames_m.clear();
    multiFieldmaps_m.clear();
    multiScales_m.clear();
    multiPhases_m.clear();
    multiFrequencies_m.clear();
}

void RFCavity::setFieldMapFN(std::string fn) {
    multiFilenames_m.push_back(fn);
    filename_m = multiFilenames_m[0];
}

string RFCavity::getFieldMapFN() const {
    std::string filename("No_fieldmap");
    if(numFieldmaps() > 0) {
        filename = multiFilenames_m[0];
    }
    return filename;
}

void RFCavity::setAmplitudem(double vPeak) {
    multiScales_m.push_back(vPeak);
    scale_m = multiScales_m[0];
}

void RFCavity::setFrequency(double freq) {
    frequency_m = freq;
}

void RFCavity::setFrequencym(double freq) {
    multiFrequencies_m.push_back(freq);
    frequency_m = multiFrequencies_m[0];
}


double RFCavity::getFrequencym() const {
    return frequency_m;
}

void RFCavity::setPhasem(double phase) {
    multiPhases_m.push_back(phase);
    phase_m = multiPhases_m[0];
}

void RFCavity::updatePhasem(double phase) {
    if(multiPhases_m.size() == 0) {
        multiPhases_m.push_back(phase);
    } else {
        double diff = phase - multiPhases_m[0];
        multiPhases_m[0] = phase;
        for(size_t i = 1; i < numFieldmaps(); ++ i) {
            multiPhases_m[i] += diff;
        }
    }
    phase_m = multiPhases_m[0];
}

double RFCavity::getPhasem() const {
    return phase_m;
}

void RFCavity::setCavityType(std::string type) {

}

string RFCavity::getCavityType() const {
    return "SW";
}

void RFCavity::setFast(bool fast) {
    fast_m = fast;
}


bool RFCavity::getFast() const {
    return fast_m;
}



void RFCavity::setAutophaseVeto(bool veto) {
    autophaseVeto_m = veto;
}


bool RFCavity::getAutophaseVeto() const {
    return autophaseVeto_m;
}


/**
 * ENVELOPE COMPONENT for radial focussing of the beam
 * Calculates the transverse envelope component for the RF cavity
 * element and adds it to the K vector
*/
void RFCavity::addKR(int i, double t, Vector_t &K) {

    double pz = RefPartBunch_m->getZ(i) - startField_m - ds_m;
    const Vector_t tmpR(RefPartBunch_m->getX(i) - dx_m, RefPartBunch_m->getY(i) - dy_m, pz);
    double k = -Physics::q_e / (2.0 * RefPartBunch_m->getGamma(i) * Physics::EMASS);

    for(size_t j = 0; j < numFieldmaps(); ++ j) {
        Fieldmap *fieldmap = multiFieldmaps_m[j];
        double frequency = multiFrequencies_m[j];
        double scale = multiScales_m[j];
        double phase = multiPhases_m[j];

        Vector_t tmpE(0.0, 0.0, 0.0), tmpB(0.0, 0.0, 0.0);
        fieldmap->getFieldstrength(tmpR, tmpE, tmpB);
        double Ez = tmpE(2);

        tmpE = Vector_t(0.0);
        fieldmap->getFieldDerivative(tmpR, tmpE, tmpB, DZ);

        double wtf = frequency * t + phase;
        double kj = k * scale * (tmpE(2) * cos(wtf) - RefPartBunch_m->getBeta(i) * frequency * Ez * sin(wtf) / Physics::c);
        K += Vector_t(kj, kj, 0.0);
    }
}


/**
 * ENVELOPE COMPONENT for transverse kick (only has an impact if x0, y0 != 0)
 * Calculates the transverse kick component for the RF cavity element and adds it to
 * the K vector. Only important for off track tracking, otherwise KT = 0.
*/
void RFCavity::addKT(int i, double t, Vector_t &K) {

    RefPartBunch_m->actT();

    //XXX: BET parameter, default is 1.
    //If cxy != 1, then cxy = true
    bool cxy = false; // default
    double kx = 0.0, ky = 0.0;
    if(cxy) {
        for(size_t j = 0; j < numFieldmaps(); ++ j) {
            Fieldmap *fieldmap = multiFieldmaps_m[j];
            double scale = multiScales_m[j];
            double frequency = multiFrequencies_m[j];
            double phase = multiPhases_m[j];

            double pz = RefPartBunch_m->getZ(i) - startField_m - ds_m;
            const Vector_t tmpA(RefPartBunch_m->getX(i) - dx_m, RefPartBunch_m->getY(i) - dy_m, pz);

            Vector_t tmpE(0.0, 0.0, 0.0), tmpB(0.0, 0.0, 0.0);
            fieldmap->getFieldstrength(tmpA, tmpE, tmpB);

            double cwtf = cos(frequency * t + phase);
            double cf = -Physics::q_e / (RefPartBunch_m->getGamma(i) * Physics::m_e);
            kx += -cf * scale * tmpE(0) * cwtf;
            ky += -cf * scale * tmpE(1) * cwtf;
        }
    }

    double dx = RefPartBunch_m->getX0(i) - dx_m;
    double dy = RefPartBunch_m->getY0(i) - dy_m;

    Vector_t KR(0.0, 0.0, 0.0);
    addKR(i, t, KR);
    //FIXME ?? different in bet src
    K += Vector_t(KR(1) * dx + kx, KR(1) * dy + ky, 0.0);
    //
    //K += Vector_t(kx - KR(1) * dx, ky - KR(1) * dy, 0.0);
}


bool RFCavity::apply(const size_t &i, const double &t, double E[], double B[]) {
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

bool RFCavity::apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B) {
    bool out_of_bounds = true;
    const Vector_t tmpR(RefPartBunch_m->getX(i) - dx_m, RefPartBunch_m->getY(i) - dy_m , RefPartBunch_m->getZ(i) - startField_m - ds_m);

    if (tmpR(2) >= 0.0) {
        // inside the cavity
        for(size_t j = 0; j < numFieldmaps(); ++ j) {
            Fieldmap *fieldmap = multiFieldmaps_m[j];
            Vector_t tmpE(0.0, 0.0, 0.0), tmpB(0.0, 0.0, 0.0);
            const std::pair<double, double> & start_end = multi_start_end_field_m[j];

            if(tmpR(2) > start_end.first &&
               tmpR(2) < start_end.second &&
               !fieldmap->getFieldstrength(tmpR, tmpE, tmpB)) {
                const double phase = multiFrequencies_m[j] * t + multiPhases_m[j];
                const double &scale = multiScales_m[j];
                E += scale * cos(phase) * tmpE;
                B -= scale * sin(phase) * tmpB;
                out_of_bounds = false;
            }
        }
    }
    else {
        /*
           some of the bunch is still outside of the cavity
           so let them drift in
        */
        out_of_bounds = false;
    }
    return out_of_bounds;
}

bool RFCavity::apply(const Vector_t &R, const Vector_t &centroid, const double &t, Vector_t &E, Vector_t &B) {
    bool out_of_bounds = true;
    const Vector_t tmpR(R(0) - dx_m, R(1) - dy_m, R(2) - startField_m - ds_m);

    for(size_t i = 0; i < numFieldmaps(); ++ i) {
        Fieldmap *fieldmap = multiFieldmaps_m[i];
        Vector_t tmpE(0.0, 0.0, 0.0), tmpB(0.0, 0.0, 0.0);
        const std::pair<double, double> & start_end = multi_start_end_field_m[i];

        if(tmpR(2) > start_end.first &&
           tmpR(2) < start_end.second &&
           !fieldmap->getFieldstrength(tmpR, tmpE, tmpB)) {
            const double phase = multiFrequencies_m[i] * t + multiPhases_m[i];
            const double &scale = multiScales_m[i];
            E += scale * cos(phase) * tmpE;
            B -= scale * sin(phase) * tmpB;
            out_of_bounds = false;
        }
    }
    return out_of_bounds;
}

void RFCavity::initialise(PartBunch *bunch, double &startField, double &endField, const double &scaleFactor) {
    using Physics::two_pi;
    double zBegin = 0.0, zEnd = 0.0, rBegin = 0.0, rEnd = 0.0;
    Inform msg("RFCavity ", *gmsg);
    std::stringstream errormsg;
    RefPartBunch_m = bunch;

    for(size_t i = 0; i < numFieldmaps(); ++ i) {
        std::string fname = multiFilenames_m[i];
        Fieldmap *fmap = Fieldmap::getFieldmap(fname, fast_m);
        if(fmap != NULL) {
            multiFieldmaps_m.push_back(fmap);
        } else {
            for(size_t j = i; j < numFieldmaps() - 1; ++ j) {
                multiFieldmaps_m[j] = multiFieldmaps_m[j + 1];
            }
            multiFieldmaps_m.pop_back();
            -- i;
        }
    }
    if(multiScales_m.size() > 0) {
        while(multiScales_m.size() < numFieldmaps()) {
            multiScales_m.push_back(multiScales_m[0]);
        }
    } else {
        while(multiScales_m.size() < numFieldmaps()) {
            multiScales_m.push_back(1.0);
        }
    }
    if(multiFrequencies_m.size() > 0) {
        while(multiFrequencies_m.size() < numFieldmaps()) {
            multiFrequencies_m.push_back(multiFrequencies_m[0]);
        }
    } else {
        while(multiFrequencies_m.size() < numFieldmaps()) {
            multiFrequencies_m.push_back(0.0);
        }
    }
    if(multiPhases_m.size() > 0) {
        while(multiPhases_m.size() < numFieldmaps()) {
            multiPhases_m.push_back(multiPhases_m[0]);
        }
    } else {
        while(multiPhases_m.size() < multiFilenames_m.size()) {
            multiPhases_m.push_back(0.0);
        }
    }

    if(numFieldmaps() > 0) {
        double overall_zBegin = 999999999.99, overall_zEnd = -999999999.99;
        for(size_t i = 0; i < numFieldmaps(); ++ i) {
            Fieldmap *fmap = multiFieldmaps_m[i];
            double &frequency = multiFrequencies_m[i];
            const std::string &filename = multiFilenames_m[i];

            fmap->getFieldDimensions(zBegin, zEnd, rBegin, rEnd);
            if(zEnd > zBegin) {
                msg << level2 << getName() << " using file ";
                fmap->getInfo(&msg);
                if(fabs((frequency - fmap->getFrequency()) / frequency) > 0.01) {
                    errormsg << "FREQUENCY IN INPUT FILE DIFFERENT THAN IN FIELD MAP '" << filename << "';\n"
                             << frequency / two_pi * 1e-6 << " MHz <> "
                             << fmap->getFrequency() / two_pi * 1e-6 << " MHz; TAKE ON THE LATTER";
                    std::string errormsg_str = Fieldmap::typeset_msg(errormsg.str(), "warning");
                    ERRORMSG(errormsg_str << "\n" << endl);
                    if(Ippl::myNode() == 0) {
                        std::ofstream omsg("errormsg.txt", std::ios_base::app);
                        omsg << errormsg_str << std::endl;
                        omsg.close();
                    }
                    frequency = fmap->getFrequency();
                }
                multi_start_end_field_m.push_back(std::pair<double, double>(zBegin, zEnd));

                overall_zBegin = std::min(overall_zBegin, zBegin);
                overall_zEnd = std::max(overall_zEnd, zEnd);
            } else {
                for(size_t j = i; j < numFieldmaps() - 1; ++ j) {
                    multiFilenames_m[j] = multiFilenames_m[j + 1];
                    multiFieldmaps_m[j] = multiFieldmaps_m[j + 1];
                    multiScales_m[j] = multiScales_m[j + 1];
                    multiPhases_m[j] = multiPhases_m[j + 1];
                    multiFrequencies_m[j] = multiFrequencies_m[j + 1];
                }
                multiFilenames_m.pop_back();
                multiFieldmaps_m.back() = NULL;
                multiFieldmaps_m.pop_back();
                multiScales_m.pop_back();
                multiPhases_m.pop_back();
                multiFrequencies_m.pop_back();
                -- i;
            }
        }
        zBegin = overall_zBegin;
        zEnd = overall_zEnd;
        for(size_t i = 0; i < numFieldmaps(); ++ i) {
            multi_start_end_field_m[i].first -= zBegin;
            multi_start_end_field_m[i].second -= zBegin;
        }
    }

    if(zEnd > zBegin) {
        ElementEdge_m = startField;
        startField_m = startField = ElementEdge_m + zBegin;
        endField_m = endField = ElementEdge_m + zEnd;
    } else {
        endField = startField - 1e-3;
    }
}

// In current version ,this function reads in the cavity voltage profile data from file.
void RFCavity::initialise(PartBunch *bunch, const double &scaleFactor, std::shared_ptr<AbstractTimeDependence> freq_atd,
                          std::shared_ptr<AbstractTimeDependence> ampl_atd, std::shared_ptr<AbstractTimeDependence> phase_atd) {

    using Physics::pi;

    RefPartBunch_m = bunch;

    /// set the time dependent models
    setAmplitudeModel(ampl_atd);
    setPhaseModel(phase_atd);
    setFrequencyModel(freq_atd);

    ifstream in(filename_m.c_str());
    if(!in.good()) {
        throw GeneralClassicException("RFCavity::initialise",
                                      "failed to open file '" + filename_m + "', please check if it exists");
    }
    *gmsg << "* Read cavity voltage profile data" << endl;

    in >> num_points_m;

    RNormal_m  = std::unique_ptr<double[]>(new double[num_points_m]);
    VrNormal_m = std::unique_ptr<double[]>(new double[num_points_m]);
    DvDr_m     = std::unique_ptr<double[]>(new double[num_points_m]);

    for(int i = 0; i < num_points_m; i++) {
        if(in.eof()) {
            throw GeneralClassicException("RFCavity::initialise",
                                          "not enough data in file '" + filename_m + "', please check the data format");
        }
        in >> RNormal_m[i] >> VrNormal_m[i] >> DvDr_m[i];

        VrNormal_m[i] *= RefPartBunch_m->getQ();
        DvDr_m[i]     *= RefPartBunch_m->getQ();
    }
    sinAngle_m = sin(angle_m / 180.0 * pi);
    cosAngle_m = cos(angle_m / 180.0 * pi);

    if (frequency_td_m)
      *gmsg << "* timedependent frequency model " << frequency_name_m << endl;
    
    *gmsg << "* Cavity voltage data read successfully!" << endl;
}

void RFCavity::finalise()
{}

bool RFCavity::bends() const {
    return false;
}


void RFCavity::goOnline(const double &) {
    std::vector<string>::iterator fmap_it;
    for(fmap_it = multiFilenames_m.begin(); fmap_it != multiFilenames_m.end(); ++ fmap_it) {
        Fieldmap::readMap(*fmap_it);
    }
    online_m = true;
}

void RFCavity::goOffline() {
    std::vector<string>::iterator fmap_it;
    for(fmap_it = multiFilenames_m.begin(); fmap_it != multiFilenames_m.end(); ++ fmap_it) {
        Fieldmap::freeMap(*fmap_it);
    }
    online_m = false;
}

void  RFCavity::setRmin(double rmin) {
    rmin_m = rmin;
}

void  RFCavity::setRmax(double rmax) {
    rmax_m = rmax;
}

void  RFCavity::setAzimuth(double angle) {
    angle_m = angle;
}

void  RFCavity::setPerpenDistance(double pdis) {
    pdis_m = pdis;
}

void  RFCavity::setGapWidth(double gapwidth) {
    gapwidth_m = gapwidth;
}

void RFCavity::setPhi0(double phi0) {
    phi0_m = phi0;
}

double  RFCavity::getRmin() const {
    return rmin_m;
}

double  RFCavity::getRmax() const {
    return rmax_m;
}

double  RFCavity::getAzimuth() const {
    return angle_m;
}

double  RFCavity::getSinAzimuth() const {
    return sinAngle_m;
}

double  RFCavity::getCosAzimuth() const {
    return cosAngle_m;
}

double  RFCavity::getPerpenDistance() const {
    return pdis_m;
}

double  RFCavity::getGapWidth() const {
    return gapwidth_m;
}

double RFCavity::getPhi0() const {
    return phi0_m;
}

void RFCavity::setComponentType(std::string name) {
    if(name == "STANDING") {
        type_m = SW;
    } else if(name == "SINGLEGAP") {
        type_m = SGSW;
    } else {
        if(name != "") {
            std::stringstream errormsg;
            errormsg << "CAVITY TYPE " << name << " DOES NOT EXIST; \n"
                     << "CHANGING TO REGULAR STANDING WAVE";
            std::string errormsg_str = Fieldmap::typeset_msg(errormsg.str(), "warning");
            ERRORMSG(errormsg_str << "\n" << endl);
            if(Ippl::myNode() == 0) {
                ofstream omsg("errormsg.txt", ios_base::app);
                omsg << errormsg_str << endl;
                omsg.close();
            }
        }
        type_m = SW;
    }

}

string RFCavity::getComponentType()const {
    if(type_m == SGSW)
        return std::string("SINGLEGAP");
    else
        return std::string("STANDING");
}

double RFCavity::getCycFrequency()const {
    return  frequency_m;
}



/**
   \brief used in OPAL-cycl

   Is called from OPAL-cycl and can handle 
   time dependent frequency, amplitude and phase

   At the moment (test) only the frequence is time
   dependent

 */
void RFCavity::getMomentaKick(const double normalRadius, double momentum[], const double t, const double dtCorrt, const int PID, const double restMass, const int chargenumber) {
    using Physics::two_pi;
    using Physics::pi;
    using Physics::c;
    double derivate;
    double Voltage;

    double momentum2  = momentum[0] * momentum[0] + momentum[1] * momentum[1] + momentum[2] * momentum[2];
    double betgam = sqrt(momentum2);

    double gamma = sqrt(1.0 + momentum2);
    double beta = betgam / gamma;

    Voltage = spline(normalRadius, &derivate) * scale_m * 1.0e6; // V

    double transit_factor = 0.0;
    double Ufactor = 1.0;

    double frequency_save = frequency_m;
    frequency_m *= frequency_td_m->getValue(t);

    if(gapwidth_m > 0.0) {
        transit_factor = 0.5 * frequency_m * gapwidth_m * 1.0e-3 / (c * beta);
        Ufactor = sin(transit_factor) / transit_factor;
    }

    Voltage *= Ufactor;

    double dgam = 0.0;
    double nphase = (frequency_m * (t + dtCorrt) * 1.0e-9) - phi0_m / 180.0 * pi ; // rad/s, ns --> rad

    dgam = Voltage * cos(nphase) / (restMass);

    double tempdegree = fmod(nphase * 360.0 / two_pi, 360.0);
    if(tempdegree > 270.0) tempdegree -= 360.0;

    gamma += dgam;

    double newmomentum2 = pow(gamma, 2) - 1.0;

    double pr = momentum[0] * cosAngle_m + momentum[1] * sinAngle_m;
    double ptheta = sqrt(newmomentum2 - pow(pr, 2));
    double px = pr * cosAngle_m - ptheta * sinAngle_m ; // x
    double py = pr * sinAngle_m + ptheta * cosAngle_m; // y

    double rotate = -derivate * (scale_m * 1.0e6) / ((rmax_m - rmin_m) / 1000.0) * sin(nphase) / (frequency_m * two_pi) / (betgam * restMass / c / chargenumber); // radian

    /// B field effects
    momentum[0] =  cos(rotate) * px + sin(rotate) * py;
    momentum[1] = -sin(rotate) * px + cos(rotate) * py;

    if(PID == 0) {
      *gmsg << "* Cavity " << getName() << " Phase= " << tempdegree << " [deg] transit time factor=  " << Ufactor
            << " dE= " << dgam *restMass * 1.0e-6 << " [MeV]"
            << " E_kin= " << (gamma - 1.0)*restMass * 1.0e-6 << " [MeV] Time dep freq = " << frequency_td_m->getValue(t) << endl;
    }
    frequency_m = frequency_save;
}

/* cubic spline subrutine */
double RFCavity::spline(double z, double *za) {
    double splint;

    // domain-test and handling of case "1-support-point"
    if(num_points_m < 1) {
        throw GeneralClassicException("RFCavity::spline",
                                      "no support points!");
    }
    if(num_points_m == 1) {
        splint = RNormal_m[0];
        *za = 0.0;
        return splint;
    }

    // search the two support-points
    int il, ih;
    il = 0;
    ih = num_points_m - 1;
    while((ih - il) > 1) {
        int i = (int)((il + ih) / 2.0);
        if(z < RNormal_m[i]) {
            ih = i;
        } else if(z > RNormal_m[i]) {
            il = i;
        } else if(z == RNormal_m[i]) {
            il = i;
            ih = i + 1;
            break;
        }
    }

    double x1 =  RNormal_m[il];
    double x2 =  RNormal_m[ih];
    double y1 =  VrNormal_m[il];
    double y2 =  VrNormal_m[ih];
    double y1a = DvDr_m[il];
    double y2a = DvDr_m[ih];
    //
    // determination of the requested function-values and its derivatives
    //
    double dx  = x2 - x1;
    double dy  = y2 - y1;
    double u   = (z - x1) / dx;
    double u2  = u * u;
    double u3  = u2 * u;
    double dy2 = -2.0 * dy;
    double ya2 = y2a + 2.0 * y1a;
    double dy3 = 3.0 * dy;
    double ya3 = y2a + y1a;
    double yb2 = dy2 + dx * ya3;
    double yb4 = dy3 - dx * ya2;
    splint = y1  + u * dx * y1a +       u2 * yb4 +        u3 * yb2;
    *za    =            y1a + 2.0 * u / dx * yb4 + 3.0 * u2 / dx * yb2;
    // if(m>=1) za=y1a+2.0*u/dx*yb4+3.0*u2/dx*yb2;
    // if(m>=2) za[1]=2.0/dx2*yb4+6.0*u/dx2*yb2;
    // if(m>=3) za[2]=6.0/dx3*yb2;

    return splint;
}

void RFCavity::getDimensions(double &zBegin, double &zEnd) const {
    zBegin = startField_m;
    zEnd = endField_m;
}


ElementBase::ElementType RFCavity::getType() const {
    return RFCAVITY;
}

double RFCavity::getAutoPhaseEstimate(const double &E0, const double &t0, const double &q, const double &mass) {
    vector<double> t, E, t2, E2;
    std::vector< std::vector< double > > F(numFieldmaps());
    std::vector< std::pair< double, double > > G;
    std::vector< double > begin(numFieldmaps());
    std::vector< double > end(numFieldmaps());
    std::vector<gsl_spline *> onAxisInterpolants(numFieldmaps());
    std::vector<gsl_interp_accel *> onAxisAccel(numFieldmaps());

    unsigned int N;
    double A, B;
    double phi = 0.0, tmp_phi, dphi = 0.5 * Physics::pi / 180.;
    double min_dz = 1.0, max_length = 0.0, min_begin = 99999999.99;
    double max_frequency = 0.0;
    for(size_t i = 0; i < numFieldmaps(); ++ i) {
        multiFieldmaps_m[i]->getOnaxisEz(G);
        begin[i] = (G.front()).first;
        end[i] = (G.back()).first;
        max_frequency = std::max(max_frequency, multiFrequencies_m[i]);
        std::unique_ptr<double[]> zvals(new double[G.size()]);
        std::unique_ptr<double[]> onAxisField(new double[G.size()]);

        for(size_t j = 0; j < G.size(); ++ j) {
            zvals[j] = G[j].first;
            onAxisField[j] = G[j].second;
        }
        onAxisInterpolants[i] = gsl_spline_alloc(gsl_interp_cspline, G.size());
        onAxisAccel[i] = gsl_interp_accel_alloc();
        gsl_spline_init(onAxisInterpolants[i], zvals.get(), onAxisField.get(), G.size());

        double length = end[i] - begin[i];
        min_dz = std::min(min_dz, length / G.size());
        max_length = std::max(max_length, length);
        min_begin = std::min(min_begin, begin[i]);

        G.clear();
        //~ delete[] zvals;
        //~ delete[] onAxisField;
    }

    N = (int)floor(max_length / min_dz + 1);
    min_dz = max_length / N;

    for(size_t i = 0; i < numFieldmaps(); ++ i) {
        F[i].resize(N);
        double z = min_begin;
        for(size_t j = 0; j < N; ++ j, z += min_dz) {
            if(z >= begin[i] && z <= end[i]) {
                F[i][j] = gsl_spline_eval(onAxisInterpolants[i], z, onAxisAccel[i]);
            } else {
                F[i][j] = 0.0;
            }
        }
        gsl_spline_free(onAxisInterpolants[i]);
        gsl_interp_accel_free(onAxisAccel[i]);
    }

    //~ delete[] onAxisInterpolants;
    //~ delete[] onAxisAccel;

    if(N > 0) {

        t.resize(N, t0);
        t2.resize(N, t0);
        E.resize(N, E0);
        E2.resize(N, E0);

        double z = min_begin + min_dz;
        for(unsigned int i = 1; i < N; ++ i, z += min_dz) {
            E[i] = E[i - 1] + min_dz * scale_m / mass;
            E2[i] = E[i];
        }

        for(int iter = 0; iter < 10; ++ iter) {
            A = B = 0.0;
            for(unsigned int i = 1; i < N; ++ i) {
                t[i] = t[i - 1] + getdT(i, E, min_dz, mass);
                t2[i] = t2[i - 1] + getdT(i, E2, min_dz, mass);
                for(size_t j = 0; j < numFieldmaps(); ++ j) {
                    const double &frequency = multiFrequencies_m[j];
                    const double &scale = multiScales_m[j];
                    const double Dphi = multiPhases_m[j] - multiPhases_m[0];
                    A += scale * (1. + frequency * (t2[i] - t[i]) / dphi) * getdA(i, t, min_dz, Dphi, frequency, F[j]);
                    B += scale * (1. + frequency * (t2[i] - t[i]) / dphi) * getdB(i, t, min_dz, Dphi, frequency, F[j]);
                }
            }

            if(fabs(B) > 0.0000001) {
                tmp_phi = atan(A / B);
            } else {
                tmp_phi = Physics::pi / 2;
            }
            if(q * (A * sin(tmp_phi) + B * cos(tmp_phi)) < 0) {
                tmp_phi += Physics::pi;
            }

            if(fabs(phi - tmp_phi) < max_frequency * (t[N - 1] - t[0]) / (10 * N)) {
                for(unsigned int i = 1; i < N; ++ i) {
                    E[i] = E[i - 1];
                    for(size_t j = 0; j < numFieldmaps(); ++ j) {
                        const double &frequency = multiFrequencies_m[j];
                        const double &scale = multiScales_m[j];
                        const double Dphi = multiPhases_m[j] - multiPhases_m[0];
                        E[i] += q * scale * getdE(i, t, min_dz, phi + Dphi, frequency, F[j]) ;
                    }
                }
                const int prevPrecision = Ippl::Info->precision(8);
                INFOMSG(level2 << "estimated phase= " << tmp_phi << " rad, "
                        << "Ekin= " << E[N - 1] << " MeV" << setprecision(prevPrecision) << endl);

                return tmp_phi;
            }
            phi = tmp_phi - floor(tmp_phi / Physics::two_pi + 0.5) * Physics::two_pi;


            for(unsigned int i = 1; i < N; ++ i) {
                E[i] = E[i - 1];
                E2[i] = E2[i - 1];
                for(size_t j = 0; j < numFieldmaps(); ++ j) {
                    const double &frequency = multiFrequencies_m[j];
                    const double &scale = multiScales_m[j];
                    const double Dphi = multiPhases_m[j] - multiPhases_m[0];
                    E[i] += q * scale * getdE(i, t, min_dz, phi + Dphi, frequency, F[j]) ;
                    E2[i] += q * scale * getdE(i, t2, min_dz, phi + dphi + Dphi, frequency, F[j]);
                }
                t[i] = t[i - 1] + getdT(i, E, min_dz, mass);
                t2[i] = t2[i - 1] + getdT(i, E2, min_dz, mass);

                E[i] = E[i - 1];
                E2[i] = E2[i - 1];
                for(size_t j = 0; j < numFieldmaps(); ++ j) {
                    const double &frequency = multiFrequencies_m[j];
                    const double &scale = multiScales_m[j];
                    const double Dphi = multiPhases_m[j] - multiPhases_m[0];
                    E[i] += q * scale * getdE(i, t, min_dz, phi + Dphi, frequency, F[j]) ;
                    E2[i] += q * scale * getdE(i, t2, min_dz, phi + dphi + Dphi, frequency, F[j]);
                }
            }

            double totalEz0 = 0.0, cosine_part = 0.0, sine_part = 0.0;
            double p0 = sqrt((E0 / mass + 1) * (E0 / mass + 1) - 1);
            for(size_t j = 0; j < numFieldmaps(); ++ j) {
                const double &frequency = multiFrequencies_m[j];
                const double &scale = multiScales_m[j];
                const double Dphi = multiPhases_m[j] - multiPhases_m[0];
                cosine_part += scale * cos(frequency * t0 + Dphi) * F[j][0];
                sine_part += scale * sin(frequency * t0 + Dphi) * F[j][0];
            }
            totalEz0 = cos(phi) * cosine_part - sin(phi) * sine_part;

            if(p0 + q * totalEz0 * (t[1] - t[0]) * Physics::c / mass < 0) {
                // make totalEz0 = 0
                tmp_phi = atan(cosine_part / sine_part);
                if(abs(tmp_phi - phi) > Physics::pi) {
                    phi = tmp_phi + Physics::pi;
                } else {
                    phi = tmp_phi;
                }
            }
        }

        const int prevPrecision = Ippl::Info->precision(8);
        INFOMSG(level2 << "estimated phase= " << tmp_phi << " rad, "
                << "Ekin= " << E[N - 1] << " MeV" << setprecision(prevPrecision) << endl);

        return phi;
    } else {
        return 0.0;
    }
}

pair<double, double> RFCavity::trackOnAxisParticle(const double &p0,
        const double &t0,
        const double &dt,
        const double &q,
        const double &mass) {
    double p = p0;
    double t = t0;
    double cdt = Physics::c * dt;
    vector<pair<double, double> > F;
    vector<size_t> length_fieldmaps(numFieldmaps());
    for(size_t i = 0; i < numFieldmaps(); ++ i) {
        vector<pair<double, double> > G;
        multiFieldmaps_m[i]->getOnaxisEz(G);
        length_fieldmaps[i] = G.size();
        F.insert(F.end(), G.begin(), G.end());
    }

    std::unique_ptr<double[]> zvals(new double[F.size()]);
    std::unique_ptr<double[]> onAxisField(new double[F.size()]);
    for(unsigned int i = 0; i < F.size(); ++i) {
        zvals[i] = F[i].first;
        onAxisField[i] = F[i].second;
    }
    double zbegin = zvals[0];
    double zend = zvals[F.size() - 1];
    std::vector<std::pair<double, double> > begin_end(numFieldmaps());
    std::vector<gsl_spline *> onAxisInterpolants(numFieldmaps());
    std::vector<gsl_interp_accel *> onAxisAccel(numFieldmaps());
    size_t start = 0;
    for(size_t i = 0; i < numFieldmaps(); ++ i) {
        begin_end[i].first = zvals[start];
        onAxisInterpolants[i] = gsl_spline_alloc(gsl_interp_cspline, length_fieldmaps[i]);
        onAxisAccel[i] = gsl_interp_accel_alloc();
        gsl_spline_init(onAxisInterpolants[i], &(zvals[start]), &(onAxisField[start]), length_fieldmaps[i]);
        start += length_fieldmaps[i];
        begin_end[i].second = zvals[start - 1];

        zbegin = std::min(zbegin, begin_end[i].first);
        zend = std::max(zend, begin_end[i].second);
    }

    //~ delete[] zvals;
    //~ delete[] onAxisField;

    double z = zbegin;
    double dz = 0.5 * p / sqrt(1.0 + p * p) * cdt;
    while(z + dz < zend && z + dz > zbegin) {
        z += dz;
        for(size_t i = 0; i < numFieldmaps(); ++ i) {
            if(z >= begin_end[i].first && z <= begin_end[i].second) {
                double phase = multiFrequencies_m[i] * t + multiPhases_m[i];
                double scale = multiScales_m[i];
                double ez = scale * gsl_spline_eval(onAxisInterpolants[i], z, onAxisAccel[i]);
                p += cos(phase) * ez * q * cdt / mass;
            }
        }
        dz = 0.5 * p / sqrt(1.0 + p * p) * cdt;
        z += dz;
        t += dt;
    }

    for(size_t i = 0; i < numFieldmaps(); ++ i) {
        gsl_spline_free(onAxisInterpolants[i]);
        gsl_interp_accel_free(onAxisAccel[i]);
    }
    //~ delete[] onAxisInterpolants;
    //~ delete[] onAxisAccel;

    const double beta = sqrt(1. - 1 / (p * p + 1.));
    const double tErr  = (z - zend) / (Physics::c * beta);

    return pair<double, double>(p, t - tErr);
}
