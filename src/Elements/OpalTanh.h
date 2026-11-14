
#ifndef OPAL_OPALTANH_H
#define OPAL_OPALTANH_H

#include "Elements/OpalElement.h"

/** OpalTanh provides user interface information for the Tanh end field model
 */
class OpalTanh : public OpalElement {
public:
    /** enum maps string to integer value for UI definitions */
    enum {
        LENGTH = COMMON,
        FRINGE_LENGTH,
        SIZE  // size of the enum
    };

    /** Default constructor initialises UI parameters. */
    OpalTanh();

    /** Destructor does nothing */
    virtual ~OpalTanh() {}

    /** Inherited copy constructor */

    /** Update the ScalingFFA with new parameters from UI parser */
    virtual void update();

private:
    // Not implemented.
    OpalTanh(const OpalTanh&);
    void operator=(const OpalTanh&);
};
#endif  // OPAL_OPALENGE_H
