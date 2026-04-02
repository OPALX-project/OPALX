#ifndef CLASSIC_Laser_HH
#define CLASSIC_Laser_HH

#include "AbsBeamline/Component.h"

class Laser : public Component {
public:
    explicit Laser(const std::string& name);
    Laser();
    Laser(const Laser& right);
    ~Laser() override;

    void accept(BeamlineVisitor&) const override;

    void initialise(PartBunch_t* bunch, double& startField, double& endField) override;
    void finalise() override;
    bool bends() const override;
    void getDimensions(double& zBegin, double& zEnd) const override;
    ElementType getType() const override;
    int getRequiredNumberOfTimeSteps() const override;

    void setWavelength(double wavelength);
    double getWavelength() const;

    void setPulseEnergy(double pulseEnergy);
    double getPulseEnergy() const;

    void setPulseLength(double pulseLength);
    double getPulseLength() const;

    void setWaistX(double waistX);
    double getWaistX() const;

    void setWaistY(double waistY);
    double getWaistY() const;

    void setDirection(const Vector_t<double, 3>& direction);
    const Vector_t<double, 3>& getDirection() const;

    void setStokes(const Vector_t<double, 3>& stokes);
    const Vector_t<double, 3>& getStokes() const;

private:
    double startField_m;
    double wavelength_m;
    double pulseEnergy_m;
    double pulseLength_m;
    double waistX_m;
    double waistY_m;
    Vector_t<double, 3> direction_m;
    Vector_t<double, 3> stokes_m;
};

inline int Laser::getRequiredNumberOfTimeSteps() const {
    return 1;
}

#endif  // CLASSIC_Laser_HH
