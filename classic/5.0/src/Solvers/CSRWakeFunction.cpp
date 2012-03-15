#include "Solvers/CSRWakeFunction.hh"
#include "Physics/Physics.h"
#include "AbsBeamline/RBend.h"
#include "AbsBeamline/SBend.h"

extern Inform *gmsg;

CSRWakeFunction::CSRWakeFunction(const std::string &name, ElementBase* element, vector<Filter*> filters, const unsigned int &N):
    WakeFunction(name, element),
    filters_m(filters.begin(), filters.end()),
    lineDensity_m(),
    dlineDensitydz_m(),
    d2lineDensitydz2_m(),
    NBin_m(N),
    R_m(0.0)
{ }

void CSRWakeFunction::apply(PartBunch &bunch)
{
    if (R_m < 1e-8) {
        double End;
        if (dynamic_cast<RBend*>(element_ref_m)) {
            RBend* bend = dynamic_cast<RBend*>(element_ref_m);
            R_m = bend->getR();
            bend->getDimensions(Begin_m, End);
        } else if (dynamic_cast<SBend*>(element_ref_m)) {
            SBend* bend = dynamic_cast<SBend*>(element_ref_m);
            R_m = bend->getR();
            bend->getDimensions(Begin_m, End);
        }
        Length_m = End - Begin_m;        
    }
    double prefactor = 1. / (4. * Physics::pi * Physics::epsilon_0);
    double pref1 = -2. * prefactor / pow(3. * R_m*R_m, 1./3.);
    double pref2 = -4. * prefactor / R_m;

    double Phi, Phi_m, Psi, SlippageLength;
    double min_relative_s, Ds_max, Ds_max2;
    double frac, x, dx1, dx2, dx3, dx4;
    double scaleFactor = Physics::c * bunch.getdT();
    Vector_t smin, smax;
    const double hz = bunch.getMesh().get_meshSpacing(2) * scaleFactor;
    const double hzsup = pow(hz,2./3.);
    const Vector_t origin = bunch.getMesh().get_origin() * scaleFactor;

#ifdef CSRDEBUG
    vector<double> temp_hist;
    static unsigned long counter = 0;
    stringstream filename_str;
    const int every = 50;
    bool print_criterion = (counter + 1) % every == 0 && Ippl::myNode() == 0;
#endif

    bunch.get_bounds(smin, smax);

    bunch.calcLineDensity();
    bunch.getLineDensity(lineDensity_m);

    const unsigned int N = lineDensity_m.size();

#ifdef CSRDEBUG
    if (print_criterion)
        temp_hist.assign(lineDensity_m.begin(), lineDensity_m.end());
#endif

    if (Ez_m.size() != N) {
        Ez_m.resize(N, 0.0);
    }

    for (vector<Filter*>::const_iterator fit = filters_m.begin(); fit != filters_m.end(); ++ fit) {
        (*fit)->apply(lineDensity_m);
    }
    dlineDensitydz_m.assign(lineDensity_m.begin(), lineDensity_m.end());
    filters_m.back()->calc_derivative(dlineDensitydz_m, hz);

    Ez_m[0] = 0.0;
    Phi_m = Length_m / R_m;

    min_relative_s = smin(2) * scaleFactor - Begin_m;
    for (int i = 1; i < N; ++i)
        {
            Phi = (min_relative_s + i * hz) / R_m;
            SlippageLength = Phi * Phi * Phi * R_m / 24.;
            Ez_m[i] = 0.0;
            if (Phi > 0.0) {
                if (Phi < Phi_m) {
                    if (SlippageLength > i * hz) {
                        dx1 = pow(i, 2./3.);
                        dx2 = pow(i, 5./3.);
                        dx3 = pow(i - 1., 5./3.);
                        Ez_m[i] += 0.3 * hzsup * dlineDensitydz_m[0] * (5. * dx1 - 3. * dx2 + 3. * dx3);
                        for (int j = 1; j < i; ++ j) {
                            dx1 = dx2;
                            dx2 = dx3;
                            dx3 = pow(i - j - 1., 5./3.);
                            Ez_m[i] += 0.3 * hzsup * dlineDensitydz_m[j] * (dx1 - 2.* dx2 + dx3);
                        }
                        Ez_m[i] += 0.9 * hzsup * dlineDensitydz_m[i];

                    } else if (SlippageLength < hz) {
                        if (4. * SlippageLength < i * hz) {
                            int j = i - static_cast<int>(floor(4. * SlippageLength / hz));
                            frac = 4. * SlippageLength / hz - (i - j);
                            Ez_m[i] -= (frac * lineDensity_m[j-1] + (1. - frac) * lineDensity_m[j]) / pow(SlippageLength, 1./3.);
                        }

                        frac = SlippageLength / hz;
                        Ez_m[i] += (frac * lineDensity_m[i-1] + (1. - frac) * lineDensity_m[i]) / pow(SlippageLength, 1./3.);

                        Ez_m[i] += 0.3 * pow(SlippageLength, 2./3.) * (5. * dlineDensitydz_m[i] - 2. * frac * (dlineDensitydz_m[i] - dlineDensitydz_m[i-1]));

                    } else {
                        if (4. * SlippageLength < i * hz) {
                            int j = i - static_cast<int>(floor(4. * SlippageLength / hz));
                            frac = 4. * SlippageLength / hz - (i - j);
                            Ez_m[i] -= (frac * lineDensity_m[j-1] + (1. - frac) * lineDensity_m[j]) / pow(SlippageLength, 1./3.);
                        }

                        int j = i - static_cast<int>(floor(SlippageLength / hz));
                        frac = SlippageLength / hz - (i - j);
                        Ez_m[i] += (frac * lineDensity_m[j-1] + (1. - frac) * lineDensity_m[j]) / pow(SlippageLength, 1./3.);

                        dx1 = pow(i - j + frac, 2./3.);
                        dx2 = pow(i - j, 2./3.);
                        dx3 = pow(i - j + frac, 5./3.);
                        dx4 = pow(i - j, 5./3.);

                        Ez_m[i] += 1.5 * hzsup * dlineDensitydz_m[j-1] * (dx1 - dx2);
                        Ez_m[i] += 0.3 * hzsup * (dlineDensitydz_m[j] - dlineDensitydz_m[j-1]) * (5.*(dx1 - dx2) + 3.*(dx3 - dx4));

                        dx1 = dx2;
                        dx2 = dx4;
                        dx3 = pow(i - j - 1., 5./3.);
                        Ez_m[i] += 0.3 * hzsup * dlineDensitydz_m[j] * (5.*dx1 - 3.*dx2 + 3.*dx3);
                        for (int k = j + 1; k < i; ++ k) {
                            dx1 = dx2;
                            dx2 = dx3;
                            dx3 = pow(i - k - 1., 5./3.);
                            Ez_m[i] += 0.9 * hzsup * dlineDensitydz_m[k] * (dx1 - 2.*dx2 + dx3);
                        }
                        Ez_m[i] += 0.9 * hzsup * dlineDensitydz_m[i];
                    }
                    Ez_m[i] *= pref1;
                } else {
                    x = Phi - Phi_m;
                    Ds_max = R_m * Phi_m*Phi_m*Phi_m / 24. * (4. - 3.* Phi_m / Phi);

                    if (Ds_max > hz * i) {
                        Psi = calcPsi(x, hz * i);
                        Ez_m[i] += 0.5 * dlineDensitydz_m[0] / (Psi + 2.*x) * hz;
                        for (int j = 1; j < i; ++ j) {
                            Psi = calcPsi(x, hz * (i - j));
                            Ez_m[i] += dlineDensitydz_m[j] / (Psi + 2.*x) * hz;
                        }
                        Ez_m[i] += 0.25 * dlineDensitydz_m[i] / x * hz;

                    } else if (Ds_max < hz) {
                        Ds_max2 = R_m * Phi_m*Phi_m / 6. * (3.* Phi - 2. * Phi_m);
                        if (Ds_max2 < i * hz) {
                            int j = i - static_cast<int>(floor(Ds_max2/hz));
                            frac = Ds_max2 / hz - (i - j);
                            Ez_m[i] += (frac * lineDensity_m[j-1] + (1. - frac) * lineDensity_m[j]) / (2. * Phi - Phi_m);
                        }
                        int j = i - static_cast<int>(floor(Ds_max/hz));
                        frac = Ds_max / hz - (i - j);
                        Ez_m[i] += (frac * lineDensity_m[j-1] + (1.-frac) * lineDensity_m[j]) / (2. * Phi - Phi_m);

                    } else {
                        Ds_max2 = R_m * Phi_m*Phi_m / 6. * (3. * Phi - 2. * Phi_m);
                        if (Ds_max2 < i * hz) {
                            int j = i - static_cast<int>(floor(Ds_max2/hz));
                            frac = Ds_max2/hz - (i - j);
                            Ez_m[i] -= (frac * lineDensity_m[j-1] + (1. - frac) * lineDensity_m[j]) / (2. * Phi - Phi_m);
                        }

                        int j = i - static_cast<int>(floor(Ds_max/hz));
                        Ez_m[i] += lineDensity_m[j-1]/(2. * Phi - Phi_m);

                        Psi = calcPsi(x, hz * (i - j + 1.));
                        Ez_m[i] += 0.5 * dlineDensitydz_m[j-1] / (Psi + 2.*x) * hz;
                        for (int k = j; k < i - 1; ++ k) {
                            Psi = calcPsi(x, hz * (i - k));
                            Ez_m[i] += dlineDensitydz_m[k] / (Psi + 2.*x) * hz;
                        }
                        Psi = calcPsi(x, hz);
                        Ez_m[i] += 0.5 * dlineDensitydz_m[i-1] / (Psi + 2.*x) * hz;

                        frac = hz / Psi;
                        Ez_m[i] += dlineDensitydz_m[i-1] * log(Psi/(2.*x) + 1.) * frac
                            + (dlineDensitydz_m[i] - dlineDensitydz_m[i-1]) * ((1. + 2.*x / Psi) * log(Psi/(2. * x) + 1.) - 1) * frac;
                    }
                    Ez_m[i] *= pref2;
                }
            }
        }

    for (int i = 0; i < bunch.getLocalNum(); ++i) {
        Vector_t R = bunch.R[i] * scaleFactor;
        int indexz = (int)floor((R(2) - origin(2))/hz);
        double leverz = (R(2) - origin(2))/hz - indexz;
        bunch.Ef[i](2) += (1.-leverz) * Ez_m[indexz] + leverz * Ez_m[indexz + 1];
    }

#ifdef CSRDEBUG
    if (print_criterion) {
        static unsigned int file_number = 0;
        ++ file_number;
        filename_str << "CSRWake" << file_number << ".txt";
        ofstream csr(filename_str.str().c_str());
        for (int i = 0; i < N; ++ i) {
            csr << i * hz << "\t"
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

double CSRWakeFunction::calcPsi(const double& x, const double& Ds) const
{
    /** solve the equation
     *  \f[
     *  \Delta s = \frac{R \Psi^3}{24} \frac{\Psi + 4x}{\Psi + x}
     *  \f]
     *  for \f$\Psi\f$ using Newtons method.
     */

    const int Nmax = 100;
    const double eps = 1e-10;
    double residual = 0.0;
    double psi = pow(24. * Ds / R_m, 1./3.);

    for (int i = 0; i < Nmax; ++i)
        {
            residual = R_m * psi*psi*psi*psi + 4. * x * R_m * psi*psi*psi - 24. * Ds * psi - 24. * Ds * x;
            if (fabs(residual) < eps)
                return psi;
            psi -= residual / (4. * R_m * psi*psi*psi + 12. * x * R_m * psi*psi - 24. * Ds);
        }
    cerr << "In CSRWakeFunction::calcPsi(): exceed maximum number of iterations!" << endl;
    return psi;
}

const string CSRWakeFunction::getType() const
{
    return "CSRWakeFunction";
}


