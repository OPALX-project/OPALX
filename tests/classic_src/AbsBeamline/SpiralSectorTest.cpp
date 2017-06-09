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

#include <cmath>
#include <fstream>

#include "gtest/gtest.h"


#include "Classic/AbsBeamline/EndFieldModel/Tanh.h"
#include "Classic/AbsBeamline/SpiralSector.h"

class SpiralSectorTest : public ::testing::Test { 
public: 
    SpiralSectorTest() : sector_m(NULL), fout_m("spiralsectortest.out") { 
    }
 
    void SetUp( ) { 
        sector_m = new SpiralSector("test");
        // characteristic length is R*dphi => 0.6545 m
        endfieldmodel::Tanh* tanh = new endfieldmodel::Tanh(M_PI/12., M_PI/96., 10);
        sector_m->setEndField(tanh);
        sector_m->setTanDelta(0.); //tan(M_PI/3.));
        sector_m->setR0(10.);
        sector_m->setDipoleConstant(1.);
        sector_m->setFieldIndex(6);
        sector_m->setMaxOrder(4);
        sector_m->calculateDfCoefficients();
    }
 
    void TearDown( ) { 
        delete sector_m;
        sector_m = NULL;
    }

    ~SpiralSectorTest() {
    }

    SpiralSector* sector_m;
    std::ofstream fout_m;

    double get_div_B(Vector_t pos, Vector_t delta) {
        Vector_t mom(0., 0., 0.);
        Vector_t E(0., 0., 0.);
        Vector_t B(0., 0., 0.);
        double t = 0;
        Vector_t div_elements(0., 0., 0.);
        // dB_x/dx ~ (B_x(x+dx) - B_x(x-dx))/2/dx
        for (size_t i = 0; i < 3; ++i) {
            pos[i] += delta[i];
            sector_m->apply(pos, mom, t, E, B);
            div_elements[i] += B[i];
            pos[i] -= 2.*delta[i];
            sector_m->apply(pos, mom, t, E, B);
            div_elements[i] -= B[i];
            pos[i] += delta[i];
            div_elements[i] /= delta[i]*2.;
        }
        std::cerr << "\n\n";
        sector_m->apply(pos, mom, t, E, B);
        std::cerr << "Div elements " << div_elements << std::endl;
        return div_elements[0] + div_elements[1] + div_elements[2];
    }

    double get_div_B_cylindrical(Vector_t posCyl, Vector_t delta) {
        Vector_t mom(0., 0., 0.);
        Vector_t E(0., 0., 0.);
        Vector_t B(0., 0., 0.);
        double t = 0;
        Vector_t div_elements(0., 0., 0.);
        // dB_x/dx ~ (B_x(x+dx) - B_x(x-dx))/2/dx
        for (size_t i = 0; i < 3; ++i) {
            posCyl[i] += delta[i];
            sector_m->applyCylindrical(posCyl, mom, t, E, B);
            div_elements[i] += B[i];
            posCyl[i] -= 2.*delta[i];
            sector_m->applyCylindrical(posCyl, mom, t, E, B);
            div_elements[i] -= B[i];
            posCyl[i] += delta[i];
            div_elements[i] /= delta[i]*2.;
        }
        sector_m->applyCylindrical(posCyl, mom, t, E, B);
        div_elements[0] += B[0]/posCyl[0];
        div_elements[1] /= posCyl[0];
        return div_elements[0] + div_elements[1] + div_elements[2];
    }

private:

};

Vector_t get_curl_B(Vector_t pos, Vector_t delta) {
    return Vector_t();
}

TEST_F(SpiralSectorTest, ApplyTest) {
    std::vector<std::vector<double> > coeff = sector_m->getDfCoefficients();
    std::cerr << "Coefficients" << std::endl;
    for (size_t n = 0; n < coeff.size(); ++n) {
        for (size_t i = 0; i < coeff[n].size(); ++i) {
            std::cerr << coeff[n][i] << " ";
        }
        std::cerr << std::endl;
    }
    std::cerr << std::endl;
    for (double y = 0.01; y < 0.055; y += 0.1 /*0.01*/) {
        for (double r = 10.; r < 12.1; r += 10. /*0.5*/) {
                for (double phi = M_PI/12./*-M_PI/3.*/; phi < M_PI/3.; phi += M_PI/*/600.*/) {
                double x = r*cos(phi);
                double z = r*sin(phi);
                Vector_t pos(r, y, phi);
                Vector_t posCart(x, y, z);
                Vector_t mom(0., 0., 0.);
                Vector_t E(0., 0., 0.);
                Vector_t B(0., 0., 0.);
                Vector_t BCart(0., 0., 0.);
                double t = 0;
                
                sector_m->applyCylindrical(pos, mom, t, E, B);
                //sector_m->apply(posCart, mom, t, E, BCart);
                double div_B_cyl = 0.; //get_div_B_cylindrical(pos, Vector_t(1e-3, 1e-3, 1e-3));
                double div_B_cart = get_div_B(posCart, Vector_t(1e-3, 1e-3, 1e-3));
                fout_m << x << "    " << y << "    " << z << "    " << r << "    " << phi << "    " << B[0] << "    " << B[1] << "    "
                       << B[2] << "    " << div_B_cyl << "    " << BCart[0] << "    " << BCart[1] << "    "
                       << BCart[2] << "     " << div_B_cart << std::endl;
            }
        }
        fout_m << std::endl;
    }
}


