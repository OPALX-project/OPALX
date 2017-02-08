#include "Utilities/OpalSection.h"
#include "Structure/OpalWake.h"
#include "Structure/SurfacePhysics.h"
#include "Solvers/WakeFunction.hh"
#include "Solvers/SurfacePhysicsHandler.hh"
#include "Structure/BoundaryGeometry.h"

extern Inform *gmsg;

OpalSection::OpalSection(const CompVec &elements, const double &start, const double &end):
    elements_m(elements.begin(), elements.end()),
    start_m(start),
    end_m(end),
    bends_m(false),
    has_wake_m(false),
    has_boundarygeometry_m(false),
    has_surface_physics_m(false),
    is_live_m(false),
    wakefunction_m(NULL),
    sphys_handler_m(NULL),
    boundarygeometry_m(NULL),
    orientation_m(0.0),
    exit_face_angle_m(0.0),
    previous_is_glued_m(false),
    glued_to_m(NULL),
    StartCache_m(),
    EndCache_m()
{
    for(CompVec::const_iterator clit = elements_m.begin(); clit != elements_m.end(); ++ clit) {
        if((*clit)->bends()) {
            bends_m = true;
            // (*clit)->getOrientation(orientation_m, exit_face_angle_m);
        }
        if((*clit)->hasWake()) {
            if(has_wake_m && wakefunction_m != (*clit)->getWake()) {
            /*--------- Modified by Xiaoying Pang 04/22/2014 ---------------
             * If two overlapping elements in one section both have wake functions,
             * the first wake function is applied for now. Later, we will implement
             * mulitple wake functions.*/
            //    *gmsg << "more than one wake function in one section! dismiss all." << endl;
            //    wakefunction_m = NULL;
                *gmsg << " more than one wake function in one section! use the wake function from the first element." << endl;

            } else {
                wakefunction_m = (*clit)->getWake();
                wakeFunctionOwner_m = (*clit);
            }
            has_wake_m = true;
        }
        if((*clit)->hasSurfacePhysics()) {
            if(has_surface_physics_m && sphys_handler_m != (*clit)->getSurfacePhysics()) {
                *gmsg << "more than one surface physics handler in one section! dismiss all." << endl;
                sphys_handler_m = NULL;
            } else {
                sphys_handler_m = (*clit)->getSurfacePhysics();
            }
            has_surface_physics_m = true;
        } else
            has_surface_physics_m = false;

        if((*clit)->hasBoundaryGeometry()) {
            /**
               we maybe want to have a boundary geometry handler
            */
            boundarygeometry_m = (*clit)->getBoundaryGeometry();
            has_boundarygeometry_m = true;
        }
    }
    if(has_wake_m && !wakefunction_m) {
        // we then dismissed them all or there is
        // an other error
        has_wake_m = false;
    }

    if(has_surface_physics_m && !sphys_handler_m) {
        has_surface_physics_m = false;
    }

    updateStartCache();
    updateEndCache();
}

OpalSection::~OpalSection() {
    for(CompVec::iterator clit = elements_m.begin(); clit != elements_m.end(); ++ clit) {
        *clit = NULL;
    }
}

void OpalSection::setOrientation(const Vector_t &angle) {
    orientation_m = angle;
    if(glued_to_m) {
        glued_to_m->setOrientation(angle);
    }
    updateStartCache();
    updateEndCache();
}

bool OpalSection::find(std::shared_ptr<const Component> check) const {
    bool found = false;
    for(size_t index = 0; index < elements_m.size(); ++ index) {
        if(elements_m[index] == check) {
            found = true;
            break;
        }
    }
    return found;
}

void OpalSection::print(Inform &msg) const {
    std::stringstream mymsg;
    static std::string closure("------------------------------------------------------------------------------------\n");
    if(glued_to_m) {
        mymsg << "--- "
              << start_m << " m -- "
              << end_m << " m -- (glued to next) ";
        if(boundarygeometry_m)
            mymsg << " has boundary geometry start at " << boundarygeometry_m->getS() ;
        msg << mymsg.str() << closure.substr(mymsg.str().length());
    } else {
        mymsg << "--- "
              << start_m << " m -- "
              << end_m << " m -- ";
        if(boundarygeometry_m)
            mymsg  << " has boundary geometry ";

        if(hasSurfacePhysics())
            mymsg  << " has surface physics ";
        msg << mymsg.str() << closure.substr(mymsg.str().length());
    }
    for(CompVec::const_iterator clit = elements_m.begin(); clit != elements_m.end(); ++ clit) {
        msg << (*clit)->getName() << '\n';
    }
}

void OpalSection::updateStartCache() {
    if(previous_is_glued_m || cos(orientation_m(0)) < 1.e-8) {
        // Setting the factors to zero will cause start_m to be returned
        StartCache_m.u_factor = 0.0;
        StartCache_m.v_factor = 0.0;
    } else if(!doesBend()) {
        double const cosa = cos(orientation_m(0));
        double const tana = tan(orientation_m(0));
        double const tanb = tan(orientation_m(1));
        double const cosc = cos(orientation_m(2));
        double const sinc = sin(orientation_m(2));
        StartCache_m.u_factor = tanb * sinc / cosa - tana * cosc;
        StartCache_m.v_factor = -tanb * cosc / cosa - tana * sinc;
    }
}

void OpalSection::updateEndCache() {
    if(glued_to_m || cos(orientation_m(0) - exit_face_angle_m) < 1.e-8) {
        // Setting the factors to zero will cause end_m to be returned
        EndCache_m.u_factor = 0.0;
        EndCache_m.v_factor = 0.0;
    } else if(!doesBend()) {
        double const cosa = cos(orientation_m(0) - exit_face_angle_m);
        double const tana = tan(orientation_m(0) - exit_face_angle_m);
        double const tanb = tan(orientation_m(1));
        double const cosc = cos(orientation_m(2));
        double const sinc = sin(orientation_m(2));
        EndCache_m.u_factor = tanb * sinc / cosa - tana * cosc;
        EndCache_m.v_factor = -tanb * cosc / cosa - tana * sinc;
    }
}

bool OpalSection::doDipoleFieldsOverlap() const {
    if (!bends_m) return false;

    unsigned int numFieldContributions = 0;
    for (auto it = elements_m.begin(); it != elements_m.end(); ++ it) {
        switch((*it)->getType()) {
        case ElementBase::CORRECTOR:
        case ElementBase::MULTIPOLE:
        case ElementBase::RBEND:
        case ElementBase::RFCAVITY:
        case ElementBase::RFQUADRUPOLE:
        case ElementBase::SBEND:
        case ElementBase::SOLENOID:
        case ElementBase::TRAVELINGWAVE:
            ++ numFieldContributions;
            break;
        default:
            break;
        }
    }
    return (numFieldContributions > 1);
}
