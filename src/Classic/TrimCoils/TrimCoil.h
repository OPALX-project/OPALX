#ifndef TRIM_COIL_H
#define TRIM_COIL_H

/// Abstract TrimCoil class

class TrimCoil {

public:

    TrimCoil();
    /// Apply the trim coil at position r and z to Bfields br and bz
    /// Calls virtual doApplyField
    void applyField(const double r, const double z, double *br, double *bz);

    virtual ~TrimCoil() { };

private:

    /// virtual implementation of applyField
    virtual void doApplyField(const double r, const double z, double *br, double *bz) = 0;
};

#endif //TRIM_COIL_H
