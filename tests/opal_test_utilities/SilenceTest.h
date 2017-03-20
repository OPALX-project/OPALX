/*
 *  Copyright (c) 2012, Chris Rogers
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

#ifndef OPALTESTUTILITIES_SILENCETEST_H_

#include <sstream>
#include <iostream>

namespace OpalTestUtilities {
/** Shutup test output
 *
 *  If more than one is called, will shutup output on any alloc if it is loud
 *  and will make loud on any dealloc if it is quiet.
 */
class SilenceTest {
  public:
    SilenceTest(bool willSilence) {
        if (willSilence && _defaultCout == NULL ) {
            _defaultCout = std::cout.rdbuf(_debugOutput.rdbuf());
            _defaultCerr = std::cerr.rdbuf(_debugOutput.rdbuf());
        }
    }

    ~SilenceTest() { // return buffer to normal on delete
        if (_defaultCout != NULL) {
            std::cout.rdbuf(_defaultCout);
            std::cerr.rdbuf(_defaultCerr);
            _defaultCout = NULL;
            _defaultCerr = NULL;
        }
    }

  private:
    SilenceTest(); // disable default ctor
    SilenceTest(const SilenceTest& test); // disable default copy ctor

    std::ostringstream _debugOutput;
    static std::streambuf *_defaultCout;
    static std::streambuf *_defaultCerr;
};
}

#endif //OPALTESTUTILITIES_SILENCETEST_H_
