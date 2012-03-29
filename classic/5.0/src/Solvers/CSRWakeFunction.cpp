#include "Solvers/CSRWakeFunction.hh"
#include "Physics/Physics.h"
#include "AbsBeamline/RBend.h"
#include "AbsBeamline/SBend.h"

using namespace std;

extern Inform *gmsg;

CSRWakeFunction::CSRWakeFunction(const std::string &name, ElementBase *element, vector<Filter *> filters, const unsigned int &N):
    WakeFunction(name, element),
    filters_m(filters.begin(), filters.end()),
    lineDensity_m(),
    dlineDensitydz_m(),
    d2lineDensitydz2_m(),
    R_m(0.0),
    Phi_m(0.0),
    NBin_m(N)
{ }

void CSRWakeFunction::apply(PartBunch &bunch) {
#ifdef CSRDEBUG
    static string bendName;
    string bendNameOld = bendName;
#endif
    if(R_m < 1e-8) {
        double End;
        if(dynamic_cast<RBend *>(element_ref_m)) {
            RBend *bend = dynamic_cast<RBend *>(element_ref_m);
            R_m = bend->getR();
            bend->getDimensions(Begin_m, End);
            Length_m = bend->getEffectiveLength();
            FieldBegin_m = Begin_m + bend->getEffectiveCenter() - Length_m / 2.0;
            Phi_m = bend->getBendAngle();

#ifdef CSRDEBUG
            bendName = bend->getName();
#endif //CSRDEBUG
        } else if(dynamic_cast<SBend *>(element_ref_m)) {
            SBend *bend = dynamic_cast<SBend *>(element_ref_m);
            R_m = bend->getR();
            bend->getDimensions(Begin_m, End);
            Length_m = bend->getEffectiveLength();
            FieldBegin_m = Begin_m + bend->getEffectiveCenter() - Length_m / 2.0;
            Phi_m = bend->getBendAngle();

#ifdef CSRDEBUG
            bendName = bend->getName();
#endif //CSRDEBUG
        }
    }

    double prefactor = 1. / (4. * Physics::pi * Physics::epsilon_0);
    double pref1 = -2. * prefactor / pow(3. * R_m * R_m, 1. / 3.);
    double pref2 = -4. * prefactor / R_m;

    double Phi, SlippageLength;
    double min_relative_s, Ds_max, Ds_max2;
    double frac, x, dx1, dx2, dx3, dx4;
    double scaleFactor = 1.0; // Physics::c * bunch.getdT();
    Vector_t smin, smax;
    const double hz = bunch.getMesh().get_meshSpacing(2) * scaleFactor;
    const double hzsup = pow(hz, -1. / 3.);
    const Vector_t origin = bunch.getMesh().get_origin() * scaleFactor;

#ifdef CSRDEBUG
    vector<double> temp_hist;
    static unsigned long counter = 0;
    if(bendNameOld != bendName) counter = 0;
    stringstream filename_str;
    const int every = 1;
    bool print_criterion = (counter + 1) % every == 0 && Ippl::myNode() == 0;
#endif

    bunch.get_bounds(smin, smax);
    bunch.calcLineDensity();
    bunch.getLineDensity(lineDensity_m);

    const unsigned int N = lineDensity_m.size();

#ifdef CSRDEBUG
    if(print_criterion)
        temp_hist.assign(lineDensity_m.begin(), lineDensity_m.end());
#endif

    if(Ez_m.size() != N) {
        Ez_m.resize(N, 0.0);
    }

    if (Psi_m.size() != N) {
        Psi_m.resize(N, 0.0);
    }

    for(vector<Filter *>::const_iterator fit = filters_m.begin(); fit != filters_m.end(); ++ fit) {
        (*fit)->apply(lineDensity_m);
    }
    dlineDensitydz_m.assign(lineDensity_m.begin(), lineDensity_m.end());
    filters_m.back()->calc_derivative(dlineDensitydz_m, hz);

    Ez_m[0] = 0.0;

    min_relative_s = smin(2) * scaleFactor - FieldBegin_m;
    for(unsigned int i = 1; i < N; ++i) {

        // Compute current bend angle of particle.
        // Particle does not bend until it enters the start
        // position of the effective dipole.
        if(min_relative_s + i * hz <= 0.0)
            Phi = 0.0;
        else
            Phi = (min_relative_s + i * hz) / R_m;

        SlippageLength = Phi * Phi * Phi * R_m / 24.;
        Ez_m[i] = 0.0;
        if(Phi > 0.0) {
            if(Phi <= Phi_m) {
                if(SlippageLength > i * hz) {

                    // Break integral into sum of integrals between grid points, then use linear interpolation
                    // between each grid point.

                    dx1 = pow(i, 2. / 3.);
                    dx2 = pow(i, 5. / 3.);
                    dx3 = pow(i - 1., 5. / 3.);
                    Ez_m[i] += 0.3 * hzsup * dlineDensitydz_m[0] * (5. * dx1 - 3. * dx2 + 3. * dx3);
                    for(unsigned int j = 1; j < i; ++ j) {
                        dx1 = dx2;
                        dx2 = dx3;
                        dx3 = pow(i - j - 1., 5. / 3.);
                        Ez_m[i] += 0.9 * hzsup * dlineDensitydz_m[j] * (dx1 - 2.* dx2 + dx3);
                    }
                    Ez_m[i] += 0.9 * hzsup * dlineDensitydz_m[i];

                } else if(SlippageLength < hz) {

                    // First do transient term.
                    if (4.0 * SlippageLength <= hz) {

                        Ez_m[i] += 3.0 * pow(SlippageLength, 2.0 / 3.0) * (lineDensity_m[i] - lineDensity_m[i - 1]) / hz;
                        frac = SlippageLength / hz;

                    } else {

                        if (4.0 * SlippageLength < i * hz) {

                            int j = i - static_cast<int>(floor(4.0 * SlippageLength / hz));
                            frac = 4.0 * SlippageLength / hz - (i - j);
                            Ez_m[i] -= (frac * lineDensity_m[j - 1] + (1. - frac) * lineDensity_m[j]) / pow(SlippageLength, 1. / 3.);

                        }

                        frac = SlippageLength / hz;
                        Ez_m[i] += (frac * lineDensity_m[i - 1] + (1. - frac) * lineDensity_m[i]) / pow(SlippageLength, 1. / 3.);

                    }

                    // Now do steady state term.
                    Ez_m[i] += (0.3 / hz) * pow(SlippageLength, 2. / 3.) * (5. * dlineDensitydz_m[i] - 2. * frac * (dlineDensitydz_m[i] - dlineDensitydz_m[i - 1]));

                } else {

                    if(4. * SlippageLength < i * hz) {

                        int j = i - static_cast<int>(floor(4. * SlippageLength / hz));
                        frac = 4. * SlippageLength / hz - (i - j);
                        Ez_m[i] -= (frac * lineDensity_m[j - 1] + (1. - frac) * lineDensity_m[j]) / pow(SlippageLength, 1. / 3.);

                    }

                    int j = i - static_cast<int>(floor(SlippageLength / hz));
                    frac = SlippageLength / hz - (i - j);
                    Ez_m[i] += (frac * lineDensity_m[j - 1] + (1. - frac) * lineDensity_m[j]) / pow(SlippageLength, 1. / 3.);

                    dx1 = pow(i - j + frac, 2. / 3.);
                    dx2 = pow(i - j, 2. / 3.);
                    dx3 = pow(i - j + frac, 5. / 3.);
                    dx4 = pow(i - j, 5. / 3.);

                    Ez_m[i] += 1.5 * hzsup * dlineDensitydz_m[j - 1] * (dx1 - dx2);
                    Ez_m[i] += 0.3 * hzsup * (dlineDensitydz_m[j] - dlineDensitydz_m[j - 1]) * (5.*(dx1 - dx2) + 3.*(dx3 - dx4));

                    dx1 = dx2;
                    dx2 = dx4;
                    dx3 = pow(i - j - 1., 5. / 3.);
                    Ez_m[i] += 0.3 * hzsup * dlineDensitydz_m[j] * (5.*dx1 - 3.*dx2 + 3.*dx3);
                    for(unsigned int k = j + 1; k < i; ++ k) {
                        dx1 = dx2;
                        dx2 = dx3;
                        dx3 = pow(i - k - 1., 5. / 3.);
                        Ez_m[i] += 0.9 * hzsup * dlineDensitydz_m[k] * (dx1 - 2.*dx2 + dx3);
                    }
                    Ez_m[i] += 0.9 * hzsup * dlineDensitydz_m[i];
                }
                Ez_m[i] *= pref1;

            } else {
                x = Phi - Phi_m;
                Ds_max = R_m * Phi_m * Phi_m * Phi_m / 24. * (4. - 3.* Phi_m / Phi);

                // First do contribution from particles whose retarded position is
                // prior to the bend.
                Ds_max2 = R_m * Phi_m * Phi_m / 6. * (3. * Phi - 2. * Phi_m);
                int j = 0;
                if(Ds_max2 < i * hz) {
                    j = i - static_cast<int>(floor(Ds_max2 / hz));
                    frac = Ds_max2 / hz - (i - j);
                    Ez_m[i] -= (frac * lineDensity_m[j - 1] + (1. - frac) * lineDensity_m[j]) / (2. * Phi - Phi_m);
                }

                // Now do delta function contribution for particles whose retarded position
                // is in the bend.
                if (Ds_max < i * hz) {
                    j = i - static_cast<int>(floor(Ds_max / hz));
                    frac = Ds_max / hz - (i - j);
                    Ez_m[i] += (frac * lineDensity_m[j - 1] + (1.0 - frac) * lineDensity_m[j]) / (2. * Phi - Phi_m);
                }

                // Now do integral contribution for particles whose retarded position is in
                // the bend.

                // First term.
                int k = i;
                if (Ds_max < i * hz) {
                    k = j;
                    Psi_m[k] = calcPsi(Psi_m[k], x, hz * k + frac * hz);
                    if (Psi_m[k] > 0 && Psi_m[k] < Phi_m)
                        Ez_m[i] += 0.5 * (frac * dlineDensitydz_m[i - k - 1] + (1.0 - frac) * dlineDensitydz_m[i - k]) / (Psi_m[k] + 2.0 * x);
                } else {
                    Psi_m[0] = calcPsi(Psi_m[0], x, hz * i);
                    if (Psi_m[0] > 0 && Psi_m[0] < Phi_m)
                        Ez_m[i] += 0.5 * dlineDensitydz_m[0] / (Psi_m[0] + 2.0 * x);
                }

                // Do rest of integral.
                for(unsigned int l = i - k + 1; l < i; ++ l) {
                    Psi_m[l] = calcPsi(Psi_m[l], x, hz * (i - l));
                    if (Psi_m[l] > 0 && Psi_m[l] < Phi_m)
                        Ez_m[i] += dlineDensitydz_m[l] / (Psi_m[l] + 2.0 * x);
                }

                // We don't go right to the end as there is a singularity in the numerical integral that we don't quite know
                // how to deal with properly yet. This introduces a very slight error in the calculation (fractions of a percent).
                Psi_m[i] = calcPsi(Psi_m[i], x, hz / 4.0);
                if (Psi_m[i] > 0 && Psi_m[i] < Phi_m)
                    Ez_m[i] += 0.5 * dlineDensitydz_m[i] / (Psi_m[i] + 2.0 * x);

                Ez_m[i] *= pref2;
            }
        }
    }

    for(unsigned int i = 0; i < bunch.getLocalNum(); ++i) {
        Vector_t R = bunch.R[i] * scaleFactor;
        int indexz = (int)floor((R(2) - origin(2)) / hz);
        double leverz = (R(2) - origin(2)) / hz - indexz;
        bunch.Ef[i](2) += (1. - leverz) * Ez_m[indexz] + leverz * Ez_m[indexz + 1];
    }


#ifdef CSRDEBUG
    double zPositionOfBunch = bunch.get_sPos() * scaleFactor;
    if (zPositionOfBunch > 1.759) {
        *gmsg << "Phi: " << Phi << " x: " << x << " Ds_max: " << Ds_max << " hz: " << hz << endl;
    }
    if(print_criterion) {
        static unsigned int file_number = 0;
        if(bendNameOld != bendName) file_number = 0;
        ++ file_number;
        filename_str << bendName << "-CSRWake" << file_number << ".txt";
        ofstream csr(filename_str.str().c_str());
        csr << zPositionOfBunch << endl;
        for(int i = 0; i < N; ++ i) {
            csr << i *hz << "\t"
                << Ez_m[i] << "\t"
                << temp_hist[i] << "\t"
                << lineDensity_m[i] << "\t"
                << dlineDensitydz_m[i] << endl;
        }
        csr.close();
        *gmsg << "** wrote " << filename_str.str() << endl;
    }
    ++ counter;
#endif

}

double CSRWakeFunction::calcPsi(const double &psiInitial, const double &x, const double &Ds) const {
    /** solve the equation
     *  \f[
     *  \Delta s = \frac{R \Psi^3}{24} \frac{\Psi + 4x}{\Psi + x}
     *  \f]
     *  for \f$\Psi\f$ using Newtons method.
     */

    const int Nmax = 100;
    const double eps = 1e-10;
    double residual = 0.0;
    double psi = pow(24. * Ds / R_m, 1. / 3.);
    if (psiInitial != 0.0) psi = psiInitial;

    for(int i = 0; i < Nmax; ++i) {
        residual = R_m * psi * psi * psi * psi + 4. * x * R_m * psi * psi * psi - 24. * Ds * psi - 24. * Ds * x;
        if(std::abs(residual) < eps)
            return psi;
        psi -= residual / (4. * R_m * psi * psi * psi + 12. * x * R_m * psi * psi - 24. * Ds);
    }
    cerr << "In CSRWakeFunction::calcPsi(): exceed maximum number of iterations!" << endl;
    return psi;
}

const string CSRWakeFunction::getType() const {
    return "CSRWakeFunction";
}


