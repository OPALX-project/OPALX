#include "Fields/BMultipoleField.h"
#include "Physics/Physics.h"

#include "FixedAlgebra/FTps.h"
#include "FixedAlgebra/FVps.h"
#include "FixedAlgebra/FTpsMath.h"

using namespace Physics;

enum { X, PX, Y, PY, T, PT };       // !!!!! used by sbendmap.hh and quadmap.h

typedef FTps<double, 6> Series;

class RbendMap {


public:
    RbendMap() { };

    FVps<double, 6> getBodyMap(const BMultipoleField &field, double length, double beta,
                               double scale, double p, double mass,
                               FVps<double, 6> map);

    Series getHamiltonian(const BMultipoleField &field, double beta,
                          double scale, double p, double mass);

    FVps<double, 6> getEntranceFringeMap(double angle, double curve,
                                         const BMultipoleField &field,
                                         double scale, FVps<double, 6> map);

    FVps<double, 6> getExitFringeMap(double angle, double curve,
                                     const BMultipoleField &field,
                                     double scale, FVps<double, 6> map);

    FVps<double, 6> getTransformMap(const Euclid3D &euclid, double refLength,
                                    double beta, double scale, double p,
                                    double mass, FVps<double, 6> map);

    /// Construct the vector potential for a Rbend.
    Series getMultipoleMap(const BMultipoleField &);
    FVps<double, 6> getThinMultipoleMap(const BMultipoleField &field, double scale, FVps<double, 6> &m);

};

FVps<double, 6> RbendMap::getEntranceFringeMap(double angle, double curve,
        const BMultipoleField &field,
        double scale, FVps<double, 6> map) {
    // *** MISSING *** Higher order terms for entrance fringe.
    double hx = scale * field.normal(1);
    double ex = hx * tan(angle);
    double ey = hx * tan(angle + map[PX][0]);
    map[PX] += ex * map[X];
    map[PY] -= ey * map[Y];
    return map;
}

FVps<double, 6> RbendMap::getExitFringeMap(double angle, double curve,
        const BMultipoleField &field,
        double scale, FVps<double, 6> map) {
    // *** MISSING *** Higher order terms for exit fringe.
    double hx = scale * field.normal(1);
    double ex = hx * tan(angle);
    double ey = hx * tan(angle - map[PX][0]);
    map[PX] += ex * map[X];
    map[PY] -= ey * map[Y];
    return map;
}

FVps<double, 6> RbendMap::getBodyMap(const BMultipoleField &field, double length, double beta,
                                     double scale, double p, double mass,
                                     FVps<double, 6> map) {
    // Build Hamiltonian in local coordinates; substitution done later.
    // Step 1: Define variables.
    Series px = Series::makeVariable(PX);
    Series py = Series::makeVariable(PY);
    Series pt = Series::makeVariable(PT) + 1.0;

    // Step 2: Kinematic terms.
    Series x  = Series::makeVariable(X);
    Series pz = sqrt(pt * pt - px * px - py * py);
    double kin = mass / p;
    Series E  = sqrt(pt * pt + kin * kin) / beta;

    // Step 3: Vector potential  (1 + h*x) * As  in curved reference.
    Series As = getMultipoleMap(field) * scale;

    //cout << "As= " << As << endl;

    // Step 4: Finish Hamiltonian, substitute previous map,
    //         and apply result to the previous map.
    Series H = As + E - pz;

    //  cout << "H= " << H << endl;
    map = ExpMap(- H * length).substitute(map);
    return map;
}

Series  RbendMap::getHamiltonian(const BMultipoleField &field, double beta,
                                 double scale, double p, double mass) {
    // Build Hamiltonian in local coordinates; substitution done later.
    // Step 1: Define variables.
    Series px = Series::makeVariable(PX);
    Series py = Series::makeVariable(PY);
    Series pt = Series::makeVariable(PT) + 1.0;

    // Step 2: Kinematic terms.
    Series x  = Series::makeVariable(X);
    Series pz = sqrt(pt * pt - px * px - py * py);
    double kin = mass / p;
    Series E  = sqrt(pt * pt + kin * kin) / beta;

    // Step 3: Vector potential  (1 + h*x) * As  in curved reference.
    Series As = getMultipoleMap(field) * scale;

    // Step 4: Finish Hamiltonian
    Series H = As + E - pz;
    return H;
}

Series RbendMap::getMultipoleMap(const BMultipoleField &field) {
    int order = field.order();

    if(order > 0) {
        static const Series x = Series::makeVariable(X);
        static const Series y = Series::makeVariable(Y);
        Series kx = + field.normal(order) / double(order);
        Series ky = - field.skew(order)   / double(order);

        while(order > 1) {
            Series kxt = x * kx - y * ky;
            Series kyt = x * ky + y * kx;
            order--;
            kx = kxt + field.normal(order) / double(order);
            ky = kyt - field.skew(order)   / double(order);
        }

        return (x * kx - y * ky);
    } else {
        return Series(0.0);
    }



}

FVps<double, 6> RbendMap::getThinMultipoleMap(const BMultipoleField &field,  double scale,  FVps<double, 6> &m) {
    int order = field.order();

    if(order > 0) {
        static const Series x = Series::makeVariable(X);
        static const Series y = Series::makeVariable(Y);
        Series kx = + field.normal(order) / double(order);
        Series ky = - field.skew(order)   / double(order);

        while(order > 1) {
            Series kxt = x * kx - y * ky;
            Series kyt = x * ky + y * kx;
            order--;
            kx = kxt + field.normal(order) / double(order);
            ky = kyt - field.skew(order)   / double(order);
        }
        m[PX] -= kx * scale;
        m[PY] += ky * scale;
    }
    return m;
}


FVps<double, 6> RbendMap::getTransformMap(const Euclid3D &euclid, double refLength,
        double beta, double scale, double p,
        double mass, FVps<double, 6> map) {
    if(! euclid.isIdentity()) {
        double kin = mass / p;
        double refTime = refLength / beta;
        Series px = map[PX];
        Series py = map[PY];
        Series pt = map[PT] + 1.0;
        Series pz = sqrt(pt * pt - px * px - py * py);

        map[PX] = euclid.M(0, 0) * px + euclid.M(1, 0) * py + euclid.M(2, 0) * pz;
        map[PY] = euclid.M(0, 1) * px + euclid.M(1, 1) * py + euclid.M(2, 1) * pz;
        pz = euclid.M(0, 2) * px + euclid.M(1, 2) * py + euclid.M(2, 2) * pz;

        Series x = map[X] - euclid.getX();
        Series y = map[Y] - euclid.getY();
        Series x2 =
            euclid.M(0, 0) * x + euclid.M(1, 0) * y - euclid.M(2, 0) * euclid.getZ();
        Series y2 =
            euclid.M(0, 1) * x + euclid.M(1, 1) * y - euclid.M(2, 1) * euclid.getZ();
        Series s2 =
            euclid.M(0, 2) * x + euclid.M(1, 2) * y - euclid.M(2, 2) * euclid.getZ();
        Series sByPz = s2 / pz;

        Series E = sqrt(pt * pt + kin * kin);
        map[X] = x2 - sByPz * map[PX];
        map[Y] = y2 - sByPz * map[PY];
        map[T] += pt * (refTime / E  + sByPz);
        return map;
    } else
        return map;
}
