#ifndef TRIM_COIL_H
#define TRIM_COIL_H

/// Abstract TrimCoil class

class TrimCoil {

public:

    TrimCoil(double bmax, double rmin, double rmax);
    /// Apply the trim coil at position r and z to Bfields br and bz
    /// Calls virtual doApplyField
    void applyField(const double r, const double z, double *br, double *bz);

    virtual ~TrimCoil() { };

protected:

    /// Maximum B field (kG)
    double bmax_m;
    /// Minimum radius (m)
    double rmin_m;
    /// Maximum radius (m)
    double rmax_m;

private:

    /// virtual implementation of applyField
    virtual void doApplyField(const double r, const double z, double *br, double *bz) = 0;
};

#endif //TRIM_COIL_H
