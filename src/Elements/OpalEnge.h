
#ifndef OPAL_OPALENGE_H
#define OPAL_OPALENGE_H

#include "Elements/OpalElement.h"

/** OpalEnge provides user interface information for the Enge and AsymmetricEnge objects
 */
class OpalEnge : public OpalElement {
  public:
    /** enum maps string to integer value for UI definitions */
    enum {
        X0 = COMMON,
        LAMBDA,
        COEFFICIENTS,
        SIZE // size of the enum
    };

    /** Default constructor initialises UI parameters. */
    OpalEnge();

    /** Destructor does nothing */
    virtual ~OpalEnge() {}

    /** Inherited copy constructor */
    virtual OpalEnge *clone(const std::string &name);

    /** Update the ScalingFFA with new parameters from UI parser */
    virtual void update();

  private:
    // Not implemented.
    OpalEnge(const OpalEnge &);
    void operator=(const OpalEnge &);

    // Clone constructor.
    OpalEnge(const std::string &name, OpalEnge *parent);
};
#endif // OPAL_OPALENGE_H
