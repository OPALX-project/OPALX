/* 
 *  Copyright (c) 2017, Chris Rogers
 *  All rights reserved.
 *  Redistribution and use in source and binary forms, with or without 
 *  modification, are permitted provided that the following conditions are met: 
 *  1. Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer. 
 *  2. Redistributions in binary form must reproduce the above copyright notice, 
 *     this list of conditions and the following disclaimer in the documentation 
 *     and/or other materials provided with the distribution.
 *  3. Neither the name of STFC nor the names of its contributors may be used to 
 *     endorse or promote products derived from this software without specific 
 *     prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
 *  POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef OPAL_OPALSCALINGFFAGMAGNET_H
#define OPAL_OPALSCALINGFFAGMAGNET_H

#include "Elements/OpalBend.h"

/** OpalScalingFFAGMagnet provides user interface information for the SCALINGFFAG object
 *
 *  Defines three parameters - field map name, units for field, length for field
 */
class OpalScalingFFAGMagnet : public OpalElement {
  public:
    /** enum maps string to integer value for UI definitions */
    enum {
        B0 = COMMON,
        R0,
        FIELD_INDEX,
        TAN_DELTA,
        MAX_Y_POWER,
        END_LENGTH,
        CENTRE_LENGTH,
        RADIAL_NEG_EXTENT,
        RADIAL_POS_EXTENT,
        PHI_START,
        PHI_END,
        AZIMUTHAL_EXTENT,
        SIZE // size of the enum
    };

    /** Default constructor initialises UI parameters. */
    OpalScalingFFAGMagnet();

    /** Destructor does nothing */
    virtual ~OpalScalingFFAGMagnet();

    /** Inherited copy constructor */
    virtual OpalScalingFFAGMagnet *clone(const std::string &name);

    /** Fill in all registered attributes
     *
     *  Just calls fillRegisteredAttributes on the base class
     */
    virtual void fillRegisteredAttributes(const ElementBase &, ValueFlag);

    /** Update the ScalingFFAG with new parameters from UI parser */
    virtual void update();

  private:
    // Not implemented.
    OpalScalingFFAGMagnet(const OpalScalingFFAGMagnet &);
    void operator=(const OpalScalingFFAGMagnet &);

    // Clone constructor.
    OpalScalingFFAGMagnet(const std::string &name, OpalScalingFFAGMagnet *parent);
};

#endif // OPAL_OPALSCALINGFFAGMAGNET_H

