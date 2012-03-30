#include "Utilities/OpalSection.h"
#include "Structure/OpalWake.h"
#include "Structure/SurfacePhysics.h"
#include "Solvers/WakeFunction.hh"
#include "Solvers/SurfacePhysicsHandler.hh"
#include "Structure/BoundaryGeometry.h"

extern Inform *gmsg;

OpalSection::OpalSection(CompVec &elements, const double &start, const double &end):
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
    glued_to_m(NULL) {
    for(CompVec::const_iterator clit = elements_m.begin(); clit != elements_m.end(); ++ clit) {
        if((*clit)->bends()) {
            bends_m = true;
            (*clit)->getOrientation(orientation_m, exit_face_angle_m);
        }
        if((*clit)->hasWake()) {
            if(has_wake_m && wakefunction_m != (*clit)->getWake()) {
                *gmsg << "more than one wake function in one section! dismiss all." << endl;
                wakefunction_m = NULL;
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
        }

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
}

OpalSection::~OpalSection() {
    for(CompVec::iterator clit = elements_m.begin(); clit != elements_m.end(); ++ clit) {
        *clit = NULL;
    }
}

double OpalSection::getStart(const double &u, const double &v) const {
    if(previous_is_glued_m || cos(orientation_m(0)) < 1.e-8) {
        return start_m;
    } else {
        //        return start_m - (sin(orientation_m(0)) * u + tan(orientation_m(1)) * v) / cos(orientation_m(0));

        double cosa = cos(orientation_m(0));
        double tana = tan(orientation_m(0));
        double tanb = tan(orientation_m(1));
        double cosc = cos(orientation_m(2));
        double sinc = sin(orientation_m(2));

        return start_m - tanb * (-u * sinc + v * cosc) / cosa - tana * (u * cosc + v * sinc);

    }
}

double OpalSection::getEnd(const double &u, const double &v) const {
    if(glued_to_m || cos(orientation_m(0) - exit_face_angle_m) < 1.e-8) {
        return end_m;
    } else {
        //        return end_m - (sin(orientation_m(0) - exit_face_angle_m) * u + tan(orientation_m(1)) * v) / cos(orientation_m(0) - exit_face_angle_m);

        double cosa = cos(orientation_m(0) - exit_face_angle_m);
        double tana = tan(orientation_m(0) - exit_face_angle_m);
        double tanb = tan(orientation_m(1));
        double cosc = cos(orientation_m(2));
        double sinc = sin(orientation_m(2));

        return end_m - tanb * (-u * sinc + v * cosc) / cosa - tana * (u * cosc + v * sinc);
    }
}

void OpalSection::setOrientation(const Vector_t &angle) {
    orientation_m = angle;
    if(glued_to_m) {
        glued_to_m->setOrientation(angle);
    }
}

bool OpalSection::find(const Component *check) const {
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
    static string closure("------------------------------------------------------------------------------------\n");
    if(glued_to_m) {
        mymsg << "--- "
              << start_m << " m -- "
              << end_m << " m -- "
              << "alpha = " << orientation_m(0) << " -- "
              << "beta = "  << orientation_m(1) << " (glued to next) ";
        if(boundarygeometry_m)
            mymsg << " has boundary geometry start at " << boundarygeometry_m->getS() ;
        msg << mymsg.str() << closure.substr(mymsg.str().length());
    } else {
        mymsg << "--- "
              << start_m << " m -- "
              << end_m << " m -- "
              << "alpha = " << orientation_m(0) << " -- "
              << "beta = "  << orientation_m(1) << " ";
        if(boundarygeometry_m)
            mymsg  << " has boundary geometry ";
        msg << mymsg.str() << closure.substr(mymsg.str().length());

    }
    for(CompVec::const_iterator clit = elements_m.begin(); clit != elements_m.end(); ++ clit) {
        msg << (*clit)->getName() << '\n';
    }
}
