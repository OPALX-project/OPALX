#include "Utilities/OpalField.h"

OpalField::OpalField(Component* element, const double& start, const double& end):
    element_m(element),
    start_m(start),
    end_m(end),
    is_on_m(false)
{ }

OpalField::~OpalField()
{
    element_m = NULL;
}
