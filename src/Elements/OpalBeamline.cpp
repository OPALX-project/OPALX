#include "Elements/OpalBeamline.h"
#include "Utilities/OpalException.h"

using namespace std;

extern Inform *gmsg;

OpalBeamline::OpalBeamline():
    sections_m(),
    elements_m(),
    prepared_m(false),
    online_secs_m(NULL)
{
    online_sections_m = new int[5 * Ippl::getNodes()];
    // This is needed to communicate across the nodes which sections are online. We should allocate memory to store
    // (# of nodes) * (max # of sections which are online) integers.
    // While we know the first number the second depends on many things. It is *ASSUMED* that the bunch isn't ever
    // in more than 5 sections at the same time. Though it might be that this assumption is wrong. Then the
    // the communication between the nodes has to be adjusted to a scheme which works always or by using any higher
    // number.
}

OpalBeamline::~OpalBeamline() {
    elements_m.clear();
    sections_m.clear();
    delete[] online_sections_m;
    if(online_secs_m)
        delete[] online_secs_m;
}

CompVec OpalBeamline::dummy_list_m = CompVec();
OpalSection OpalBeamline::dummy_section_m = OpalSection(dummy_list_m, 0., 0.);

CompVec &OpalBeamline::getPredecessors(std::shared_ptr<const Component> element) {
    size_t index;
    bool found = false;
    for(index = 0; index < sections_m.size(); ++ index) {
        if(sections_m[index].find(element)) {
            found = true;
            break;
        }
    }

    if(found && index > 0) {
        return sections_m[index - 1].getElements();
    } else {
        if(dummy_list_m.size() != 0) {
            dummy_list_m.erase(dummy_list_m.begin(), dummy_list_m.end());
        }
        return dummy_list_m;
    }

}

CompVec &OpalBeamline::getSuccessors(std::shared_ptr<const Component> element) {
    size_t index;
    bool found = false;

    if (sections_m.size() == 0)
        return dummy_list_m;

    if (element == NULL)
        return sections_m[0].getElements();

    for(index = 0; index < sections_m.size(); ++ index) {
        if(sections_m[index].find(element)) {
            found = true;
            break;
        }
    }

    if(found && index + 1 < sections_m.size()) {
        return sections_m[index + 1].getElements();
    } else {
        if(dummy_list_m.size() != 0) {
            dummy_list_m.erase(dummy_list_m.begin(), dummy_list_m.end());
        }
        return dummy_list_m;
    }
}

OpalSection &OpalBeamline::getSectionAt(const Vector_t &pos, long &sectionIndex) {
    if(sectionIndex >= (long) sections_m.size()) sectionIndex = static_cast<long>(sections_m.size()) - 1;

    if(prepared_m) {
        if(pos(2) < sections_m[0].getStart(pos(0), pos(1))) return dummy_section_m;
        if(pos(2) > sections_m.back().getEnd(pos(0), pos(1))) {
            sectionIndex = BEAMLINE_EOL;
            return dummy_section_m;
        }
        if(sectionIndex == -1) sectionIndex = sections_m.size() - 1;

        while(pos(2) < sections_m[sectionIndex].getStart(pos(0), pos(1))) --sectionIndex;
        while(pos(2) > sections_m[sectionIndex].getEnd(pos(0), pos(1))) ++sectionIndex;

        if(pos(2) >= sections_m[sectionIndex].getStart(pos(0), pos(1)) &&
           pos(2) <= sections_m[sectionIndex].getEnd(pos(0), pos(1))) return sections_m[sectionIndex];

        return dummy_section_m;
    }
    *gmsg << "** WARNING *************************************\n"
          << "in OpalBeamline::getSectionAt(): section list not prepared\n"
          << "************************************************" << endl;
    return dummy_section_m;
}

OpalSection &OpalBeamline::getSection(const unsigned int &index) {
    if(index < sections_m.size()) {
        return sections_m[index];
    }

    return dummy_section_m;
}

void OpalBeamline::getSectionIndexAt(const Vector_t &pos, long &initial_guess) const {
    if(initial_guess > -1 && (size_t) initial_guess >= sections_m.size()) initial_guess = static_cast<long>(sections_m.size()) - 1;

    if(prepared_m) {
        if (sections_m.size() == 0) return;

        if(initial_guess == -1 && pos(2) > sections_m[0].getStart(pos(0), pos(1)))
            initial_guess = 0;

        if(initial_guess > -1) {
            if(pos(2) < sections_m[initial_guess].getStart(pos(0), pos(1))) {
                do {
                    -- initial_guess;
                } while(initial_guess > -1 && pos(2) < sections_m[initial_guess].getEnd(pos(0), pos(1)));
                ++ initial_guess;
            } else {
                while(initial_guess < (long) sections_m.size() && pos(2) > sections_m[initial_guess].getEnd(pos(0), pos(1))) {
                    ++ initial_guess;
                }
                if(pos(2) > sections_m.back().getEnd(pos(0), pos(1))) {
                    initial_guess = BEAMLINE_EOL;
                } else {
                    if(pos(2) < sections_m[initial_guess].getStart(pos(0), pos(1))) {
                        -- initial_guess;
                    }
                }
            }
        }
        return;
    }
    *gmsg << "** WARNING *************************************\n"
          << "in OpalBeamline::getSectionIndexAt(): section list not prepared\n"
          << "************************************************" << endl;
}

void OpalBeamline::getKFactors(const unsigned int &index, const Vector_t &pos, const long &sindex, const double &t, Vector_t &KR, Vector_t &KT) {

    if(sindex > -1 && (size_t) sindex < sections_m.size()) {
        OpalSection &section = getSection(sindex);
        setStatus(sindex, true);
        if(pos(2) >= section.getStart(pos(0), pos(1)) && pos(2) < section.getEnd(pos(0), pos(1))) {
            const CompVec &elements = section.getElements();
            for(CompVec::const_iterator elit = elements.begin(); elit != elements.end(); elit++) {
                (*elit)->addKR(index, t, KR);
                (*elit)->addKT(index, t, KT);
            }
        }
    }
}

unsigned long OpalBeamline::getFieldAt(const unsigned int &index, const Vector_t &pos, const long &sindex, const double &t, Vector_t &E, Vector_t &B) {

    /*
      If one uses static rtv has
      after the first particles enters the Collimator
      always 16384

      static unsigned int rtv = 0x00;
     */

    if(sindex >= (long) sections_m.size()) return BEAMLINE_EOL;
    if(sindex < 0) return 0x00;

    unsigned long rtv = 0x00;

    B = Vector_t(0.0);
    E = Vector_t(0.0);

    OpalSection &section = getSection(sindex);
    setStatus(sindex, true);
    if(pos(2) < section.getStart(pos(0), pos(1)) || pos(2) > section.getEnd(pos(0), pos(1))) return 0x00;

    const CompVec &elements = section.getElements();
    for(CompVec::const_iterator elit = elements.begin(); elit != elements.end(); ++ elit) {
        if((*elit)->apply(index, t, E, B)) {
            return BEAMLINE_OOB;
        }
    }

    const Vector_t &ori = section.getOrientation();
    if(fabs(ori(0)) > 1.e-10 || fabs(ori(1)) > 1.e-10 || fabs(ori(2)) > 1.e-10) {

        // Rotate field out of the element's local coordinate system and back to lab frame.
        //
        // 1) Rotate about the x axis by angle -ori(1).
        // 2) Rotate about the y axis by angle ori(0).
        // 3) Rotate about the z axis by angle ori(3).

        // Tenzor<double, 3> rota(cos(ori(0)), 0, sin(ori(0)), 0, 1, 0, -sin(ori(0)), 0, cos(ori(0)));
        // Tenzor<double, 3> rotb(1, 0, 0, 0, cos(ori(1)), sin(ori(1)), 0, -sin(ori(1)), cos(ori(1)));
        // Tenzor<double, 3> rotc(cos(ori(2)), -sin(ori(2)), 0, sin(ori(2)), cos(ori(2)), 0, 0, 0, 1);

        // Tenzor<double, 3> orientation = dot(rotc, dot(rota, rotb));
        // E = dot(orientation, E);
        // B = dot(orientation, B);

        const  double sina = sin(ori(0)),
                      cosa = cos(ori(0)),
                      sinb = sin(ori(1)),
                      cosb = cos(ori(1)),
                      sinc = sin(ori(2)),
                      cosc = cos(ori(2));

        Vector_t temp = E;

        E(0) =  cosa * cosc * temp(0) + (-sina * sinb * cosc - cosb * sinc) * temp(1) + (sina * cosb * cosc - sinb * sinc) * temp(2);
        E(1) =  cosa * sinc * temp(0) + (-sina * sinb * sinc + cosb * cosc) * temp(1) + (sina * cosb * sinc + sinb * cosc) * temp(2);
        E(2) = -sina *        temp(0) + (-cosa * sinb) * temp(1) + (cosa * cosb) * temp(2);
        temp = B;
        B(0) =  cosa * cosc * temp(0) + (-sina * sinb * cosc - cosb * sinc) * temp(1) + (sina * cosb * cosc - sinb * sinc) * temp(2);
        B(1) =  cosa * sinc * temp(0) + (-sina * sinb * sinc + cosb * cosc) * temp(1) + (sina * cosb * sinc + sinb * cosc) * temp(2);
        B(2) = -sina *        temp(0) + (-cosa * sinb) * temp(1) + (cosa * cosb) * temp(2);

    }

    if(section.doesBend()) {
        rtv |= BEAMLINE_BEND;
    }
    if(section.hasWake()) {
        rtv |= BEAMLINE_WAKE;
    }
    if(section.hasSurfacePhysics()) {
        rtv |= BEAMLINE_SURFACEPHYSICS;
    }

    return rtv;
}

unsigned long OpalBeamline::getFieldAt(const Vector_t &pos, const Vector_t &centroid, const double &t, Vector_t &E, Vector_t &B) {
    unsigned long rtv = 0x00;

    long initial_guess = -1;
    B = Vector_t(0.0);
    E = Vector_t(0.0);

    OpalSection &section = getSectionAt(pos, initial_guess);
    if(!(initial_guess & BEAMLINE_EOL || initial_guess < 0)) {
        unsigned int ini_guess = initial_guess;
        setStatus(ini_guess, true);

        if(pos(2) >= section.getStart(pos(0), pos(1)) && pos(2) <= section.getEnd(pos(0), pos(1))) {
            const CompVec &elements = section.getElements();
            for(CompVec::const_iterator elit = elements.begin(); elit != elements.end(); ++ elit) {
                if((*elit)->apply(pos, centroid, t, E, B)) {

                    rtv |= BEAMLINE_OOB;
                }
            }
            if(! rtv & BEAMLINE_OOB) {

                const Vector_t &ori = section.getOrientation();
                if(fabs(ori(0)) > 1.e-10 || fabs(ori(1)) > 1.e-10 || fabs(ori(2)) > 1.e-10) {

                    // Rotate field out of the element's local coordinate system and back to lab frame.
                    //
                    // 1) Rotate about the x axis by angle ori(1).
                    // 2) Rotate about the y axis by angle ori(0).
                    // 3) Rotate about the z axis by angle negative ori(3).

                    const  double sina = sin(ori(0)),
                                  cosa = cos(ori(0)),
                                  sinb = sin(ori(1)),
                                  cosb = cos(ori(1)),
                                  sinc = sin(ori(2)),
                                  cosc = cos(ori(2));

                    Vector_t temp = E;

                    E(0) =  cosa * cosc * temp(0) + (-sina * sinb * cosc - cosb * sinc) * temp(1) + (sina * cosb * cosc - sinb * sinc) * temp(2);
                    E(1) =  cosa * sinc * temp(0) + (-sina * sinb * sinc + cosb * cosc) * temp(1) + (sina * cosb * sinc + sinb * cosc) * temp(2);
                    E(2) = -sina *        temp(0) + (-cosa * sinb) * temp(1) + (cosa * cosb) * temp(2);
                    temp = B;
                    B(0) =  cosa * cosc * temp(0) + (-sina * sinb * cosc - cosb * sinc) * temp(1) + (sina * cosb * cosc - sinb * sinc) * temp(2);
                    B(1) =  cosa * sinc * temp(0) + (-sina * sinb * sinc + cosb * cosc) * temp(1) + (sina * cosb * sinc + sinb * cosc) * temp(2);
                    B(2) = -sina *        temp(0) + (-cosa * sinb) * temp(1) + (cosa * cosb) * temp(2);
                }
            }
            if(section.doesBend()) {
                rtv |= BEAMLINE_BEND;
            }
            if(section.hasWake()) {
                rtv |= BEAMLINE_WAKE;
            }
            if(section.hasSurfacePhysics()) {
                rtv |= BEAMLINE_SURFACEPHYSICS;
            }

            return rtv;
        }
        return 0x00;
    } else {
        if(initial_guess > -1) {
            return BEAMLINE_EOL;
        }
        return 0x00;
    }
}

void OpalBeamline::switchElements(const double &min, const double &max, const double &kineticEnergy, const bool &nomonitors) {
    for(FieldList::iterator flit = elements_m.begin(); flit != elements_m.end(); ++ flit) {
        // don't set online monitors if the centroid of the bunch is allready inside monitor
        // or if explicitly not desired (eg during auto phasing)
        if(flit->getElement()->getType() == ElementBase::MONITOR) {
            double spos = (max + min) / 2.;
            if(!nomonitors && spos < (*flit).getStart()) {
                if(!(*flit).isOn() && max > (*flit).getStart()) {
                    (*flit).setOn(kineticEnergy);
                }
            }

        } else {
            if(!(*flit).isOn() && max > (*flit).getStart() && min < (*flit).getEnd()) {
                (*flit).setOn(kineticEnergy);
            }
        }
        /////////////////////////
        // does not work like that
        //         if (min > (*flit).getEnd()) {
        //             (*flit).setOff();
        //         }
        /////////////////////////

    }
}


void OpalBeamline::switchAllElements() {
    for(FieldList::iterator flit = elements_m.begin(); flit != elements_m.end(); ++ flit) {
        if(!(*flit).isOn()) {
            (*flit).setOn(-1.0);
        }
    }

}

void OpalBeamline::switchElementsOff(const double &min, ElementBase::ElementType eltype) {
    if(eltype == ElementBase::ANY) {
        for(FieldList::iterator flit = elements_m.begin(); flit != elements_m.end(); ++ flit) {
            if((*flit).isOn() && min >= (*flit).getEnd()) {
                (*flit).setOff();
            }

        }
    } else {
        for(FieldList::iterator flit = elements_m.begin(); flit != elements_m.end(); ++ flit) {
            if((*flit).isOn() && min >= (*flit).getEnd() && (*flit).getElement()->getType() == eltype) {
                (*flit).setOff();
            }
        }
    }
}

void OpalBeamline::switchElementsOff() {
    for(FieldList::iterator flit = elements_m.begin(); flit != elements_m.end(); ++ flit)
        (*flit).setOff();
}

void OpalBeamline::prepareSections() {
    list<double> start_end;
    CompVec tmp;
    FieldList::iterator flit;
    list<double>::iterator pos_it, next_it, last_it;
    const double tolerance = 1.e-10;

    if (elements_m.size() == 0) {
        prepared_m = true;
        return;
    }

    /* there might be elements with length zero or extremely short ones.
       we delete them such that they don't appear in the simulation
     */
    elements_m.sort(ClassicField::SortAsc);
    for(flit = elements_m.begin();  flit != elements_m.end(); ++ flit) {
        while(flit != elements_m.end() && (*flit).getLength() < tolerance) {
            FieldList::iterator temp = flit;
            ++temp;
            elements_m.erase(flit);
            flit = temp;
        }
        if(flit != elements_m.end()) {
            start_end.push_back((*flit).getStart());
            start_end.push_back((*flit).getEnd());
        }
    }
    start_end.sort();

    if(start_end.size() == 0) {
        throw OpalException("OpalBeamline::prepareSections",
                            "no valid elements found");
    }

    next_it = start_end.begin();
    ++next_it;
    for(pos_it = start_end.begin(); next_it != start_end.end(); ++ pos_it, ++ next_it) {
        if(*next_it - *pos_it < tolerance) {
            *next_it = *pos_it;
        }
    }

    start_end.unique();  // remove duplicate entries

    next_it = start_end.begin();
    ++next_it;
    for(pos_it = start_end.begin(); next_it != start_end.end(); ++ pos_it, ++ next_it) {
        tmp.clear();
        for(flit = elements_m.begin(); flit != elements_m.end(); ++ flit) {
            if((*flit).getStart() <= *pos_it  + tolerance && (*flit).getEnd()   >= *next_it - tolerance) {
                tmp.push_back((*flit).getElement());
            } else {
                if((*flit).getStart() >= *next_it) {
                    break;
                }
            }
        }
        if(tmp.size() > 0) {
            sections_m.push_back(OpalSection(tmp, *pos_it, *next_it));
        }
    }

    for(int i = 0; i < -1 + static_cast<long>(sections_m.size()); ++ i) {
        if(sections_m[i].doesBend()) {
            const CompVec &els = sections_m[i].getElements();
            int num_bending_elements = 0;
            for(CompVec::const_iterator elit = els.begin(); elit != els.end(); ++ elit) {
                if((*elit)->bends()) {
                    ++ num_bending_elements;
                }
            }
            if(num_bending_elements > 1) {
                if(els.size() == 2) {
                    // TODO:
                    // set end of section i - 1 to (End_{i-1} + End_{i})/2 and
                    // set start of section i + 1 to (Start_{i+1} + Start_{i})/2
                    //
                    // mark section i for deletion
                }
            }
        }
        if (sections_m[i].doDipoleFieldsOverlap()) {
            const CompVec &elements = sections_m[i].getElements();
            auto it = elements.begin();
            std::stringstream elementNamesStream;
            std::string elementNames;
            for (; it != elements.end(); ++ it) {
                double start, end;
                (*it)->getDimensions(start, end);
                elementNamesStream << (*it)->getName() << ": "
                                   << "start= " << start << " m, "
                                   << "end= " << end << " m,\n";
            }
            elementNames = elementNamesStream.str();
            elementNames.erase(elementNames.length() - 2, 1);
            throw OpalException("OpalBeamline::prepareSections",
                                "Fields overlap with dipole fields; not supported yet;\n*** affected elements:\n" + elementNames +
                                "***");
        }
    }


    for(int i = 0; i < static_cast<long>(sections_m.size()) - 1; ++ i) {
        if(sections_m[i].doesBend() && sections_m[i + 1].doesBend()) {

            if(fabs(sections_m[i].getEnd(0.0, 0.0) - sections_m[i + 1].getStart(0.0, 0.0)) < 1.e-12) {
                // TODO:
                // check that elements realy overlap
                const CompVec &els = sections_m[i].getElements();
                for(CompVec::const_iterator elit = els.begin(); elit != els.end(); ++ elit) {
                    if((*elit)->bends()) {
                        if(sections_m[i + 1].find(*elit)) {
                            sections_m[i].glue_to(&sections_m[i + 1]);
                            sections_m[i + 1].previous_is_glued();
                        }
                    }
                }
            }
        }
    }

    online_secs_m = new bool[sections_m.size()];
    for(size_t i = 0; i < sections_m.size(); ++ i) {
        online_secs_m[i] = false;
    }
    prepared_m = true;
}

void OpalBeamline::print(Inform &msg) const {
    SectionList::const_iterator sec_it;

    msg << level1 << "\n--- BEGIN FIELD LIST ---------------------------------------------------------------\n\n";
    for(sec_it = sections_m.begin(); sec_it != sections_m.end(); ++ sec_it) {
        (*sec_it).print(msg);
    }
    msg << level1 << "\n--- END FIELD LIST -----------------------------------------------------------------\n" << endl;
}


double OpalBeamline::calcBeamlineLenght() {
    SectionList::const_iterator sec_it;
    double l = -1000000000.0;
    for(sec_it = sections_m.begin(); sec_it != sections_m.end(); ++ sec_it) {
        double tmp = (*sec_it).getEnd();
        if(tmp > l)
            l = tmp;
    }
    return l;
}


FieldList OpalBeamline::getElementByType(ElementBase::ElementType type) {
    FieldList elements_of_requested_type;
    for(FieldList::iterator fit = elements_m.begin(); fit != elements_m.end(); ++ fit) {
        if((*fit).getElement()->getType() == type) {
            elements_of_requested_type.push_back((*fit));
        }
    }
    return elements_of_requested_type;
}