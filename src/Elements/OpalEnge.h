
#ifndef OPAL_OPALENGE_H
#define OPAL_OPALENGE_H

#include "Elements/OpalElement.h"

/** OpalEnge provides user interface information for the Enge and AsymmetricEnge objects
 */
class ElementValues;

class OpalEnge : public OpalElement {
public:
    /** enum maps string to integer value for UI definitions */
    enum {
        LENGTH = COMMON,
        FRINGE_LENGTH,
        COEFFICIENTS,
        SIZE  // size of the enum
    };

    /** Default constructor initialises UI parameters. */
    OpalEnge();

    /** Destructor does nothing */
    virtual ~OpalEnge() {}

    /** Inherited copy constructor */

    /**
     * @brief Make the model of input object @p values and register it under its name, where
     *        the elements look it up. Called by buildElements() at RUN.
     */
    static void registerModel(const ElementValues& values);

private:
    // Not implemented.
    OpalEnge(const OpalEnge&);
    void operator=(const OpalEnge&);
};
#endif  // OPAL_OPALENGE_H
