// ------------------------------------------------------------------------
// $RCSfile: PlaceRep.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: PlaceRep
//   A class used to represent a place specification.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:34 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/PlaceRep.h"
#include "AbstractObjects/Element.h"
#include "AbstractObjects/OpalData.h"
#include "Utilities/OpalException.h"

// Class PlaceRep
// ------------------------------------------------------------------------

PlaceRep::PlaceRep() : data(), is_selected(false) { initialize(); }

PlaceRep::PlaceRep(const PlaceRep& rhs) : data(rhs.data), is_selected(rhs.is_selected) {
    initialize();
}

PlaceRep::PlaceRep(const std::string& def) : data(), is_selected(def == "SELECTED") {
    append(def, 1);
    initialize();
}

PlaceRep::~PlaceRep() {}

const PlaceRep& PlaceRep::operator=(const PlaceRep& rhs) {
    data = rhs.data;
    initialize();
    return *this;
}

void PlaceRep::append(const std::string& name, int occur) {
    data.push_back(std::make_pair(name, occur));
}

void PlaceRep::initialize() {
    status = false;
    seen   = 0;
}

bool PlaceRep::isActive() const { return status; }

bool PlaceRep::isSelected() const { return is_selected; }

void PlaceRep::print(std::ostream& os) const {
    if (data.empty()) {
        os << "#S";
    } else {
        os << data[0].first;
        if (data[0].second > 0) os << '[' << data[0].second << ']';

        for (Data::size_type i = 1; i < data.size(); ++i) {
            os << "::" << data[i].first;
            if (data[i].second) os << '[' << data[i].second << ']';
        }
    }

    return;
}
