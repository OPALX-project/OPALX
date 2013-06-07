#ifndef __LOMB__
#define __LOMB__
/*****************************************************************************/
/*                                                                           */
/* Class for Lomb Periodograms & Co.                                         */
/* =================================                                         */
/*                                                                           */
/*****************************************************************************/
#include <vector>
#include <cmath>

typedef struct {
    double x, y;
} LOMB_TYPE;

typedef std::vector<LOMB_TYPE>::const_iterator CI_lt;
typedef std::vector<double>::const_iterator CI_vd;
typedef std::vector<float >::const_iterator CI_vf;

class Lomb_eq : public std::unary_function<LOMB_TYPE, bool> {
    double b;
public:
    explicit Lomb_eq(const double &a) : b(a) {}
    bool operator()(const LOMB_TYPE &c) const {return c.y == b;}
};

class Compare_XPos {
public:
    bool operator()(const LOMB_TYPE &a, const LOMB_TYPE &b)
    { return a.y < b.y; }

};

class Countlomb {
public:
    bool operator()(const LOMB_TYPE &a, const double &b)
    { return a.y == b; }

};

class LOMB_class {
private:

    double TWOPID;

public:

    LOMB_class(int);          // constructor
    virtual ~LOMB_class(void); //destructor

    int period(std::vector<LOMB_TYPE> *indata, std::vector<LOMB_TYPE> *outdata,
               double ofac, double hifac, int *nout, int *jmax, double *prob,
               int amp);

    int avevar(std::vector<LOMB_TYPE> *data, double *ave, double *var);
    double Signi(double *peak, int *nout, double *ofac);

    int moment(std::vector<LOMB_TYPE> *indata, double *ave, double *adev,
               double *sdev, double *var, double *skew, double *curt);


};

#endif
