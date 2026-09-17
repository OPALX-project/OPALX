
#ifndef OPAL_OPALTANH_H
#define OPAL_OPALTANH_H

#include "Elements/OpalElement.h"

/** OpalTanh provides user interface information for the Tanh end field model
 */
class OpalTanh : public OpalElement {
  public:
    /** enum maps string to integer value for UI definitions */
    enum {
        X0 = COMMON,
        LAMBDA,
        SIZE // size of the enum
    };

    /** Default constructor initialises UI parameters. */
    OpalTanh();

    /** Destructor does nothing */
    virtual ~OpalTanh() {}

    /** Inherited copy constructor */
    virtual OpalTanh *clone(const std::string &name);

    /** Update the ScalingFFA with new parameters from UI parser */
    virtual void update();

  private:
    // Not implemented.
    OpalTanh(const OpalTanh &);
    void operator=(const OpalTanh &);

    // Clone constructor.
    OpalTanh(const std::string &name, OpalTanh *parent);
};
#endif // OPAL_OPALENGE_H
