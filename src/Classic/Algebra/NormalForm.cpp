// ------------------------------------------------------------------------
// $RCSfile: NormalForm.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: NormalForm
//   Find normal form of a truncated Taylor series map.
//   Compute form representation of a map which can be a nil-potent,
//   static or dynamic symplectic map.
//
// Implementation of an algorithm described in
//   M. Berz, E. Forest and J. Irwin,
//   Particle Accelerators, 1989, Vol. 24, pp. 91-107.
//
// ------------------------------------------------------------------------
// Class category: Algebra
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:32 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Algebra/NormalForm.h"
#include "Algebra/DoubleEigen.h"
#include "Algebra/LieMap.h"
#include "Algebra/LUMatrix.h"
#include "Algebra/TpsMonomial.h"
#include "Physics/Physics.h"
#include "Utilities/LogicalError.h"
#include "Utilities/SizeError.h"
#include <algorithm>

using std::abs;
using std::arg;
using std::imag;
using std::pow;
using std::real;
using std::swap;


// Tolerance for accepting an eigenvalue.
namespace {
    const double tol = 1.0e-8;
}


// Class NormalForm
// ------------------------------------------------------------------------

NormalForm::NormalForm():
    dimension(0),
    freedom(0),
    A_Lie(),
    N_Lie(),
    lambda(),
    V()
{}


NormalForm::NormalForm(const NormalForm &form):
    dimension(form.dimension),
    freedom(form.freedom),
    A_Lie(form.A_Lie),
    N_Lie(form.N_Lie),
    lambda(form.lambda),
    V(form.V)
{}


NormalForm::NormalForm(const VpsInvMap<double> &M_scr):
    dimension(M_scr.getDimension()),
    freedom(dimension / 2),
    A_Lie(0.0),
    N_Lie(2, dimension),   // Reserve space for second-order.
    lambda(dimension),
    V(dimension, dimension) {
    int nVar = M_scr.getVariables();
    if(nVar != dimension) {
        throw SizeError("NormalForm::NormalForm()", "Map is not R^n --> R^n.");
    }
    if(dimension % 2 != 0) {
        throw LogicalError("NormalForm::NormalForm()", "Map dimension is odd.");
    }

    // Establish Jacobian matrix of map,
    // and find its eigenvalues and eigenvectors.
    // Order eigenvectors, attach modes to planes and find inverse of V.
    Matrix<double> V_inv;
    {
        Matrix<double> M = M_scr.linearTerms();
        DoubleEigen eigen(M, true);
        orderModes(eigen.eigenValues(), eigen.packedEigenVectors());
        LUMatrix<double> lu(V);
        V_inv = lu.inverse();
    }

    // Initialise.
    Matrix<double> Rot_inv(dimension, dimension, 1.0);
    Matrix<double> R_dir(dimension, dimension, 1.0);
    Matrix<double> I_dir(dimension, dimension, 1.0);
    Matrix<double> R_inv(dimension, dimension, 1.0);
    TpsMonomial pows(dimension);

    for(int i = 0; i < 2 * freedom; i += 2) {
        // Store linear part of normal form.
        double c1, c2;

        if(imag(lambda(i)) != 0.0) {
            c1 = c2 = arg(lambda(i)) / 2.0;
        } else {
            c1 = log(std::abs(lambda(i))) / 2.0;
            c2 = - c1;
        }

        pows[i] = 2;
        N_Lie.setCoefficient(pows, c1);
        pows[i] = 0;

        pows[i+1] = 2;
        N_Lie.setCoefficient(pows, c2);
        pows[i+1] = 0;

        // Set up the auxiliary matrices.
        R_dir(i, i) = R_dir(i, i + 1) = R_dir(i + 1, i) = 0.5;
        R_dir(i + 1, i + 1) = -0.5;
        R_inv(i, i) = R_inv(i, i + 1) = R_inv(i + 1, i) = 1.0;
        R_inv(i + 1, i + 1) = -1.0;

        if(imag(lambda(i)) != 0.0) {
            // Complex eigenvalue pair.
            Rot_inv(i, i)   = Rot_inv(i + 1, i + 1) = real(lambda(i));
            Rot_inv(i, i + 1) = -imag(lambda(i));
            Rot_inv(i + 1, i) =  imag(lambda(i));

            I_dir(i, i) = I_dir(i + 1, i + 1) = 0.0;
            I_dir(i, i + 1) = I_dir(i + 1, i) = 1.0;
        } else {
            // Real eigenvalue pair.
            Rot_inv(i, i) = Rot_inv(i + 1, i + 1) = real(lambda(i) + lambda(i + 1)) * 0.5;
            Rot_inv(i, i + 1) = Rot_inv(i + 1, i) = real(lambda(i) - lambda(i + 1)) * 0.5;
        }
    }

    // Fill in the values for the coasting modes.
    // Note that the matrices already have ones on the main diagonal.
    for(int i = 2 * freedom; i < dimension; i += 2) {
        // Store linear part of normal form.
        pows[i+1] = 2;
        N_Lie.setCoefficient(pows, 0.5);
        pows[i+1] = 0;
    }

    // Initialise the factored map.
    VpsInvMap<double> N_scr = V_inv * M_scr.substitute(V * Rot_inv);
    VpsInvMap<double> N_acc(dimension);

    // Remove non-resonant terms order by order.
    int maxOrder = M_scr.getTruncOrder();

    for(int omega = 2; omega < maxOrder; omega++) {
        // Compute the terms to be removed and store in f.
        VpsInvMap<double> temp = N_scr.substitute(N_acc).filter(omega, omega);
        Tps<double> f(omega + 1, dimension);

        for(int i = 0; i < dimension; i += 2) {
            f += temp[i+1].multiplyVariable(i) - temp[i].multiplyVariable(i + 1);
        }

        f /= double(omega + 1);

        // Set up the "phi" transformations.
        Tps<double> a(omega + 1, dimension);
        Tps<double> b(omega + 1, dimension);
        Tps<double> pi(omega + 1, dimension);
        Tps<double> t(omega + 1, dimension);

        for(int m = f.getSize(omega); m < f.getSize(omega + 1); m++) {
            const TpsMonomial &index = f.getExponents(m);
            complex<double> factor = 1.0;
            bool reject = true;
            int count = 0;

            for(int j = 0; j < 2 * freedom; j += 2) {
                if(index[j] != index[j+1]) reject = false;
                if(imag(lambda(j)) != 0.0) count += index[j+1];
                factor *= pow(lambda(j), int(index[j] - index[j+1]));
            }

            if(reject  ||  std::abs(1.0 - factor) < tol) {
                // Term cannot be removed.
                t.setCoefficient(m, 0.5);
            } else {
                // Term can be removed.
                factor = 1.0 / (1.0 - factor);
                a.setCoefficient(m, real(factor));
                b.setCoefficient(m, imag(factor));
            }

            pi.setCoefficient(m, pow(-1.0, int(count + 1) / 2));
        }

        // Compute cal_T^(-1) * f and T_omega.
        Tps<double> f1 = pi.scaleMonomials(f).substitute(R_dir);
        Tps<double> f2 = f1.substitute(I_dir);

        // Remove terms of order "omega" and accumulate map script(A).
        Tps<double> F_omega =
            (a.scaleMonomials(f1).substitute(R_inv) +
             b.scaleMonomials(f2).substitute(R_inv)).scaleMonomials(pi);

        if(F_omega != 0.0) {
            LieMap<double> A_omega =
                LieMap<double>::ExpMap(F_omega.substitute(Rot_inv));
            LieMap<double> A_invert = LieMap<double>::ExpMap(- F_omega);
            N_scr = A_invert.substitute(N_scr).substitute(A_omega);
            A_Lie -= F_omega;
        }

        // Accumulate map script(N).
        Tps<double> T_omega = t.scaleMonomials(f1 + f2).
                              substitute(R_inv).scaleMonomials(pi);

        if(T_omega != 0.0) {
            N_acc = LieMap<double>::ExpMap(- T_omega, N_acc);
            N_Lie -= T_omega;
        }
    }
}


NormalForm::~NormalForm()
{}


int NormalForm::degreesOfFreedom() const {
    return freedom;
}


Matrix<double> NormalForm::anharmonicity() const {
    Matrix<double> QQ(freedom, freedom);
    using Physics::pi;

    for(int mode1 = 0; mode1 < freedom; mode1++) {
        {
            TpsMonomial power1(dimension), power2(dimension), power3(dimension);
            power1[2*mode1] = power3[2*mode1+1] = 4;
            power2[2*mode1] = power2[2*mode1+1] = 2;
            QQ(mode1, mode1) =
                - (3.0 * (N_Lie.getCoefficient(power1) +
                          N_Lie.getCoefficient(power3)) +
                   N_Lie.getCoefficient(power2)) / (4.0 * pi);
        }

        for(int mode2 = mode1 + 1; mode2 < freedom; mode2++) {
            TpsMonomial power1(dimension), power2(dimension),
                        power3(dimension), power4(dimension);
            power1[2*mode1]   = power1[2*mode2]   = 2;
            power2[2*mode1+1] = power2[2*mode2]   = 2;
            power3[2*mode1]   = power3[2*mode2+1] = 2;
            power4[2*mode1+1] = power4[2*mode2+1] = 2;
            QQ(mode1, mode2) = QQ(mode2, mode1) =
                                   - (N_Lie.getCoefficient(power1) +
                                      N_Lie.getCoefficient(power2) +
                                      N_Lie.getCoefficient(power3) +
                                      N_Lie.getCoefficient(power4)) / (4.0 * pi);
        }
    }

    return QQ;
}


Tps<double> NormalForm::invariant(int mode) const {
    Tps<double> a(1, V.nrows());
    Tps<double> b(1, V.nrows());

    for(int j = 0; j < V.nrows(); j += 2) {
        a.setCoefficient(j + 1, - V(j + 1, 2 * mode));
        a.setCoefficient(j + 2,   V(j,  2 * mode));
        b.setCoefficient(j + 1, - V(j + 1, 2 * mode + 1));
        b.setCoefficient(j + 2,   V(j,  2 * mode + 1));
    }

    if(imag(lambda(2 * mode)) != 0.0) {
        return (a * a + b * b);
    } else {
        return (a * a - b * b);
    }
}


const Tps<double> &NormalForm::normalForm() const {
    return N_Lie;
}


const Tps<double> &NormalForm::normalisingMap() const {
    return A_Lie;
}


const Vector<complex<double> > &NormalForm::eigenValues() const {
    return lambda;
}


const Matrix<double> &NormalForm::eigenVectors() const {
    return V;
}


void NormalForm::orderModes(Vector<complex<double> > tlam, Matrix<double> tmat) {
    // Copy complex eigenvalue pairs.
    int nDim = dimension;
    int n_c = 0;
    int n_r = 0;

    for(int i = 0; i < dimension;) {
        if(imag(tlam(i)) != 0.0) {
            // Normalise this complex pair.
            double pb = 0.0;

            for(int j = 0; j < dimension; j += 2) {
                pb += tmat(j, i) * tmat(j + 1, i + 1) - tmat(j + 1, i) * tmat(j, i + 1);
            }

            lambda(n_c)   = tlam(i);
            lambda(n_c + 1) = tlam(i + 1);
            if(pb < 0.0) swap(lambda(i), lambda(i + 1));
            double fact1 = sqrt(std::abs(pb));
            double fact2 = pb / fact1;

            for(int j = 0; j < dimension; j++) {
                V(j, n_c)   = tmat(j, i)   / fact1;
                V(j, n_c + 1) = tmat(j, i + 1) / fact2;
            }

            i += 2;
            n_c += 2;
        } else {
            if(std::abs(real(tlam(i)) - 1.0) < tol) {
                // Move unit eigenvalues to upper end.
                for(int j = 0; j < dimension; j++) {
                    V(j, nDim - 1) = 0.0;
                }

                lambda(nDim - 1) = 1.0;
                V(nDim - 1, nDim - 1) = 1.0;
                nDim--;
            } else {
                // Compact storage of remaining real eigenvalues.
                if(n_r != i) {
                    tlam(n_r) = tlam(i);
                    std::copy(tmat.col_begin(i), tmat.col_end(i), tmat.col_begin(n_r));
                }

                n_r++;
            }

            i++;
        }
    }

    // Copy remaining real eigenvalue pairs.
    for(int i = 0; i < n_r;) {
        int m;
        for(m = i + 1; m < n_r; m++) {
            if(std::abs(tlam(i) * tlam(m) - 1.0) < tol) break;
        }

        if(m >= n_r) {
            throw LogicalError("FNormalForm::orderModes()",
                               "Cannot find pair of real eigenvalues.");
        }

        // Swap values to make pair.
        if(m != i + 1) {
            swap(tlam(m), tlam(i + 1));
            tmat.swapColumns(m, i + 1);
        }

        // Normalise this real pair and convert it to Sinh/Cosh form.
        double pb = 0.0;

        for(int j = 0; j < dimension; j += 2) {
            pb += tmat(j, i) * tmat(j + 1, i + 1) - tmat(j + 1, i) * tmat(j, i + 1);
        }


        pb = sqrt(2.0 * std::abs(pb));
        lambda(n_c)   = tlam(i);
        lambda(n_c + 1) = tlam(i + 1);

        for(int j = 0; j < dimension; j++) {
            V(j, n_c)   = (tmat(j, i) + tmat(j, i + 1)) / pb;
            V(j, n_c + 1) = (tmat(j, i) - tmat(j, i + 1)) / pb;
        }

        if(pb < 0.0) {
            swap(lambda(i), lambda(i + 1));
            V.swapColumns(i, i + 1);
        }

        n_c += 2;
        i += 2;
    }

    freedom = nDim / 2;
    if(nDim != 2 * freedom) {
        throw LogicalError("NormalForm::orderModes()", "Map dimension is odd.");
    }

    // Re-order eigenvector pairs by their main components.
    for(int i = 0; i < nDim; i += 2) {
        // Find eigenvector pair with larges component in direction i/2.
        double big = 0.0;
        int k = i;
        for(int j = i; j < dimension; j += 2) {
            double c = std::abs(V(j, i) * V(j + 1, i + 1) - V(j, i + 1) * V(j + 1, i));

            if(c > big) {
                big = c;
                k   = j;
            }
        }

        if(k != i) {
            // Move eigenvector pair to its place.
            swap(lambda(i),   lambda(k));
            swap(lambda(i + 1), lambda(k + 1));
            V.swapColumns(i, k);
            V.swapColumns(i + 1, k + 1);
        }

        if(imag(lambda(i)) != 0.0) {
            // Rotate complex eigenvectors to make their main components real.
            double re = V(i, i)   / sqrt(V(i, i) * V(i, i) + V(i, i + 1) * V(i, i + 1));
            double im = V(i, i + 1) / sqrt(V(i, i) * V(i, i) + V(i, i + 1) * V(i, i + 1));

            for(int j = 0; j < dimension; j++) {
                double real_part = re * V(j, i) + im * V(j, i + 1);
                double imag_part = re * V(j, i + 1) - im * V(j, i);
                V(j, i)   = real_part;
                V(j, i + 1) = imag_part;
            }
        }
    }
}
