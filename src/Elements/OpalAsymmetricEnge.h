
#ifndef OPAL_OPALASYMMETRICENGE_H
#define OPAL_OPALASYMMETRICENGE_H

#include "Elements/OpalElement.h"

/** OpalEnge provides user interface information for the Enge and AsymmetricEnge objects
 */
class OpalAsymmetricEnge : public OpalElement {
  public:
    /** enum maps string to integer value for UI definitions */
    enum {
        X0_START = COMMON,
        LAMBDA_START,
        COEFFICIENTS_START,
        X0_END,
        LAMBDA_END,
        COEFFICIENTS_END,
        SIZE // size of the enum
    };

    /** Default constructor initialises UI parameters. */
    OpalAsymmetricEnge();

    /** Destructor does nothing */
    virtual ~OpalAsymmetricEnge();

    /** Inherited copy constructor */
    virtual OpalAsymmetricEnge *clone(const std::string &name);

    /** Update the ScalingFFA with new parameters from UI parser */
    virtual void update();

  private:
    // Not implemented.
    OpalAsymmetricEnge(const OpalAsymmetricEnge &);
    OpalAsymmetricEnge& operator=(const OpalAsymmetricEnge &);

    // Clone constructor.
    OpalAsymmetricEnge(const std::string &name, OpalAsymmetricEnge *parent);
};
#endif // OPAL_OPALASYMMETRICENGE_H
