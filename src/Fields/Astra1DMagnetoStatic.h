#ifndef OPALX_AstraFIELDMAP1DMAGNETOSTATIC_HH
#define OPALX_AstraFIELDMAP1DMAGNETOSTATIC_HH

#include "Fields/Fieldmap.h"
#include "Physics/Physics.h"

#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>

class Astra1DMagnetoStatic : public Fieldmap {
public:
    bool getFieldstrength(
            const Vector_t<double, 3>& R, Vector_t<double, 3>& E,
            Vector_t<double, 3>& B) const override;

    void getFieldDimensions(double& zBegin, double& zEnd) const override;

    void getFieldDimensions(
            double& xIni, double& xFinal, double& yIni, double& yFinal, double& zIni,
            double& zFinal) const override;

    bool getFieldDerivative(
            const Vector_t<double, 3>& R, Vector_t<double, 3>& E, Vector_t<double, 3>& B,
            const DiffDirection& dir) const override;

    void swap() override;

    void getInfo(Inform* msg) override;

    double getFrequency() const override;

    void setFrequency(double freq) override;

    bool isInside(const Vector_t<double, 3>& r) const override {
        return r(2) >= zbegin_m && r(2) < zend_m;
    }

    template <class ViewType>
    KOKKOS_INLINE_FUNCTION static void computeField(
            const Vector_t<double, 3>& R, Vector_t<double, 3>& /*E*/, Vector_t<double, 3>& B,
            const ViewType& FourCoefs, double zbegin, double length, int accuracy) {
        const double RR2 = R(0) * R(0) + R(1) * R(1);

        const double two_pi = Physics::two_pi;
        const double pi     = Physics::pi;
        const double dk     = Physics::two_pi / length;
        const double kz     = two_pi * (R(2) - zbegin) / length + pi;

        double bz    = FourCoefs(0);
        double bzp   = 0.0;
        double bzpp  = 0.0;
        double bzppp = 0.0;

        int n = 1;
        for (int l = 1; l < accuracy; ++l, n += 2) {
            const double base = dk * l;

            const double coskzl = Kokkos::cos(kz * l);
            const double sinkzl = Kokkos::sin(kz * l);

            bz += FourCoefs(n) * coskzl - FourCoefs(n + 1) * sinkzl;
            bzp += base * (-FourCoefs(n) * sinkzl - FourCoefs(n + 1) * coskzl);
            bzpp += base * base * (-FourCoefs(n) * coskzl + FourCoefs(n + 1) * sinkzl);
            bzppp += base * base * base * (FourCoefs(n) * sinkzl + FourCoefs(n + 1) * coskzl);
        }

        const double BfieldR = -bzp / 2.0 + bzppp * RR2 / 16.0;

        B(0) += BfieldR * R(0);
        B(1) += BfieldR * R(1);
        B(2) += bz - bzpp * RR2 / 4.0;
    }

    void applyField(std::shared_ptr<ParticleContainer_t> pc, double scale) override;

private:
    Astra1DMagnetoStatic(const std::string& filename);
    ~Astra1DMagnetoStatic();

    void readMap() override;
    void freeMap() override;

    Kokkos::DualView<double*> FourCoefs_m;

    double zbegin_m;
    double zend_m;
    double length_m;

    int accuracy_m;
    int num_gridpz_m;

    friend class Fieldmap;
};

#endif