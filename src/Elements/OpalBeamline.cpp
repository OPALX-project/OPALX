#include "Elements/OpalBeamline.h"

using namespace std;

OpalBeamline::OpalBeamline():
    elements_m(),
    sections_m(),
    prepared_m(false)
{
    
}

OpalBeamline::~OpalBeamline()
{
    elements_m.clear();
    sections_m.clear();
}

CompVec OpalBeamline::dummy_list_m = CompVec();
OpalSection OpalBeamline::dummy_section_m = OpalSection(dummy_list_m, 0., 0.);

CompVec& OpalBeamline::getPredecessors(const Component* element)
{
    int index;
    bool found = false;
    for (index = 0; index < sections_m.size(); ++ index) {
        if (sections_m[index].find(element)) {
            found = true;
            break;
        }
    }
    
    if (found && index > 0) {
        return sections_m[index - 1].getElements();
    } else {
        if (dummy_list_m.size() != 0) {
            dummy_list_m.erase(dummy_list_m.begin(), dummy_list_m.end());
        }
        return dummy_list_m;
    }
    
}

CompVec& OpalBeamline::getSuccessors(const Component* element)
{
    int index;
    bool found = false;

    for (index = 0; index < sections_m.size(); ++ index) {
        if (sections_m[index].find(element)) {
            found = true;
            break;
        }
    }
    
    if (found && index < sections_m.size() - 1) {
        return sections_m[index + 1].getElements();
    } else {
        if (dummy_list_m.size() != 0) {
            dummy_list_m.erase(dummy_list_m.begin(), dummy_list_m.end());
        }
        return dummy_list_m;
    }    
}

OpalSection& OpalBeamline::getSectionAt(const Vector_t& pos, unsigned int& initial_guess)
{
    if (initial_guess >= sections_m.size()) initial_guess = sections_m.size() - 1;

    if (prepared_m) {
        if (pos(2) < sections_m[initial_guess].getStart(0., 0.)) {
            do {
                -- initial_guess;
            } while (initial_guess > -1 && pos(2) < sections_m[initial_guess].getEnd(0., 0.));
            ++ initial_guess;
            return sections_m[initial_guess];
        } else {
            while (initial_guess < sections_m.size() && pos(2) > sections_m[initial_guess].getStart(0., 0.)) {
                ++ initial_guess;
            }
            if (pos(2) > sections_m.back().getEnd(pos(0),pos(1))) {
                initial_guess = BEAMLINE_EOL;
                if (dummy_list_m.size() != 0) {
                    dummy_list_m.erase(dummy_list_m.begin(), dummy_list_m.end());
                }
                return dummy_section_m;
            } else {
                -- initial_guess;
                return sections_m[initial_guess];
            }
        }
    } else {
        *gmsg << "** WARNING *************************************\n" 
              << "in OpalBeamline::getSectionAt(): section list no prepared\n"
              << "************************************************" << endl;
        return dummy_section_m;
    }    
}

OpalSection& OpalBeamline::getSection(const unsigned int& index)
{
    if (index < sections_m.size()) {
        return sections_m[index];
    } else {
        return dummy_section_m;
    }
}

void OpalBeamline::getSectionIndexAt(const Vector_t& pos, unsigned int& initial_guess) const
{
    if (initial_guess >= sections_m.size()) initial_guess = sections_m.size() - 1;

    if (prepared_m) {
        if (pos(2) < sections_m[initial_guess].getStart(pos(0), pos(1))) {
            do {
                -- initial_guess;
            } while (initial_guess > -1 && pos(2) < sections_m[initial_guess].getEnd(pos(0), pos(1)));
            ++ initial_guess;
        } else {
            while (initial_guess < sections_m.size() && pos(2) > sections_m[initial_guess].getStart(pos(0), pos(1))) {
                ++ initial_guess;
            }
            if (pos(2) > sections_m.back().getEnd(pos(0),pos(1))) {
                initial_guess = BEAMLINE_EOL;
            } else {
                -- initial_guess;
            }
        }
    } else {
        *gmsg << "** WARNING *************************************\n" 
              << "in OpalBeamline::getSectionIndexAt(): section list no prepared\n"
              << "************************************************" << endl;
    }
}

const unsigned int& OpalBeamline::getFieldAt(const unsigned int& index, const Vector_t& pos, const unsigned int& sindex, const double& t, Vector_t& E, Vector_t& B)
{
    
    /*
      If one uses static rtv has
      after the first particles enteres the Collimator 
      always 16384 

      static unsigned int rtv = 0x00;
    */

    unsigned int rtv = 0x00;

    B = Vector_t(0.0);
    E = Vector_t(0.0);
    if (!(sindex & BEAMLINE_EOL)) {
        OpalSection& section = getSection(sindex);
        setStatus(sindex, true);
        if (pos(2) >= section.getStart(pos(0),pos(1)) && pos(2) <= section.getEnd(pos(0),pos(1))) {
            const CompVec& elements = section.getElements();
            if (section.doesBend()) {
                rtv |= BEAMLINE_BEND;
            }
            if (section.hasWake()) {
                rtv |= BEAMLINE_WAKE;
            }
            for (CompVec::const_iterator elit = elements.begin(); elit != elements.end(); ++ elit) {
                if ((*elit)->apply(index, t, E, B)) {
                    rtv |= BEAMLINE_OOB;
                } else {
                    const Vector_t& ori = section.getOrientation();
                    if (ori(0) > 1.e-6 || ori(1) > 1.e-6) {
                        const  double sina = sin(ori(0)),
                                      cosa = cos(ori(0)),
                                      sinb = sin(ori(1)),
                                      cosb = cos(ori(1));
                        static Vector_t temp = E;
                        E(0) = cosa * temp(0) - sina * sinb * temp(1) + sina *cosb * temp(2);
                        E(1) = cosb * temp(1) + sinb * temp(2);
                        E(2) = -sina * temp(0) + cosa * sinb * temp(1) + cosa * cosb * temp(2);
                        temp = B;
                        B(0) = cosa * temp(0) - sina * sinb * temp(1) + sina *cosb * temp(2);
                        B(1) = cosb * temp(1) + sinb * temp(2);
                        B(2) = -sina * temp(0) + cosa * sinb * temp(1) + cosa * cosb * temp(2);
                    }
                }
                /*
                  This speeds up the loop in case of an lost particle
                */
                if (rtv & BEAMLINE_OOB)
                    break;
            }
        }
    } else {
        rtv |= BEAMLINE_EOL;
    }
    return rtv;
}

const unsigned int& OpalBeamline::getFieldAt(const Vector_t& pos, const double& t, Vector_t& E, Vector_t& B)
{
    static unsigned int rtv;
    static unsigned int initial_guess;
    rtv = 0x00;
    initial_guess = 0;
    B = Vector_t(0.0);
    E = Vector_t(0.0);
    OpalSection& section = getSectionAt(pos, initial_guess);
    if (!(initial_guess & BEAMLINE_EOL)) {
        setStatus(initial_guess, true);
        if (pos(2) >= section.getStart(pos(0),pos(1)) && pos(2) <= section.getEnd(pos(0),pos(1))) {
            const CompVec& elements = section.getElements();
            if (section.doesBend()) {
                rtv |= BEAMLINE_BEND;
            }
            if (section.hasWake()) {
                rtv |= BEAMLINE_WAKE;
            }
            for (CompVec::const_iterator elit = elements.begin(); elit != elements.end(); ++ elit) {
                if ((*elit)->apply(pos, t, E, B)) {
                    rtv |= BEAMLINE_OOB;
                } else {
                    const Vector_t& ori = section.getOrientation();
                    if (ori(0) > 1.e-6 || ori(1) > 1.e-6) {
                        double sina = sin(ori(0)),
                            cosa = cos(ori(0)),
                            sinb = sin(ori(1)),
                            cosb = cos(ori(1));
                        Vector_t temp = E;
                        E(0) = cosa * temp(0) - sina * sinb * temp(1) + sina *cosb * temp(2);
                        E(1) = cosb * temp(1) + sinb * temp(2);
                        E(2) = -sina * temp(0) + cosa * sinb * temp(1) + cosa * cosb * temp(2);
                        temp = B;
                        B(0) = cosa * temp(0) - sina * sinb * temp(1) + sina *cosb * temp(2);
                        B(1) = cosb * temp(1) + sinb * temp(2);
                        B(2) = -sina * temp(0) + cosa * sinb * temp(1) + cosa * cosb * temp(2);
                    }
                }
            }
        }
    } else {
        rtv |= BEAMLINE_EOL;
    }
    return rtv;
}

void OpalBeamline::switchElements(const double& min, const double& max) 
{
    for (FieldList::iterator flit = elements_m.begin(); flit != elements_m.end(); ++ flit) {
        if (!(*flit).isOn() && max > (*flit).getStart() && min < (*flit).getEnd()) {
            (*flit).setOn();
        }
/////////////////////////
// does not work like that        
//         if (min > (*flit).getEnd()) {
//             (*flit).setOff();
//         }
/////////////////////////

    }
}

void OpalBeamline::switchElementsOff(const double& min, const double& max)
{
    for (FieldList::iterator flit = elements_m.begin(); flit != elements_m.end(); ++ flit) {
        if ((*flit).isOn() && max >= (*flit).getEnd() && min <= (*flit).getStart()) {
            (*flit).setOff();
        }
    }
}

void OpalBeamline::switchElementsOff()
{
    for (FieldList::iterator flit = elements_m.begin(); flit != elements_m.end(); ++ flit) 
	(*flit).setOff();
}

void OpalBeamline::prepareSections()
{
    list<double> start_end;
    CompVec tmp;
    FieldList::iterator flit;
    list<double>::iterator pos_it, next_it, last_it;
    const double tolerance = 1.e-4;

    /* there might be elements with length zero or extremely short ones.
       we delete them such that they don't appear in the simulation
    */
    elements_m.sort(OpalField::SortAsc);

    for (flit = elements_m.begin();  flit != elements_m.end(); ++ flit) {
        while ((*flit).getLength() < tolerance && flit != elements_m.end()) {
            FieldList::iterator temp = flit;
            ++temp;
            elements_m.erase(flit);
            flit = temp;
        }
        if (flit != elements_m.end()) {
            start_end.push_back((*flit).getStart());
            start_end.push_back((*flit).getEnd());
        }
    }
    start_end.sort();

    next_it = start_end.begin(); 
    ++next_it;
    for (pos_it = start_end.begin(); next_it != start_end.end(); ++ pos_it, ++ next_it) {
        if (*next_it - *pos_it < tolerance) {
            *next_it = *pos_it;
        }
    }

    start_end.unique();  // remove duplicate entries

    next_it = start_end.begin(); 
    ++next_it;
    for (pos_it = start_end.begin(); next_it != start_end.end(); ++ pos_it, ++ next_it){
        tmp.clear();
        for (flit = elements_m.begin(); flit != elements_m.end(); ++ flit) {
            if ((*flit).getStart() <= *pos_it  + tolerance && (*flit).getEnd()   >= *next_it - tolerance) {
                tmp.push_back((*flit).getElement());
            } else {
                if ((*flit).getStart() >= *next_it) {
                    break;
                }
            }
        }
        if (tmp.size() > 0) {
            sections_m.push_back(OpalSection(tmp,*pos_it,*next_it));
        }
    }
    prepared_m = true;
}

void OpalBeamline::print(Inform& msg) const
{
    SectionList::const_iterator sec_it;

    msg << "\n--- BEGIN FIELD LIST ----------------------------------------------------------------------\n" << endl;
    for (sec_it = sections_m.begin(); sec_it != sections_m.end(); ++ sec_it) {
        (*sec_it).print(msg);
    }
    msg << "\n--- END   FIELD LIST ----------------------------------------------------------------------\n" << endl;
}


