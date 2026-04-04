#ifndef MSLANG_AFFINETRANSFORMATION_H
#define MSLANG_AFFINETRANSFORMATION_H

#include "Algorithms/Matrix.h"
#include "VectorMath.h"

#include <cmath>
#include <iostream>
#include <fstream>

namespace mslang {
     struct AffineTransformation {
        matrix3x3_t data_;
        
        AffineTransformation(const Vector_t<double, 3>& row0,
                             const Vector_t<double, 3>& row1):
          data_() {
            data_(0, 0) = row0[0];
            data_(0, 1) = row0[1];
            data_(0, 2) = row0[2];
            data_(1, 0) = row1[0];
            data_(1, 1) = row1[1];
            data_(1, 2) = row1[2];
            data_(2, 0) = 0.0;
            data_(2, 1) = 0.0;
            data_(2, 2) = 1.0;
          }

        AffineTransformation():
          data_() {
            // Default constructor already initializes to identity
          }
        
        // Access operators for compatibility
        double& operator()(int i, int j) {
            return data_(i, j);
        }
        
        double operator()(int i, int j) const {
            return data_(i, j);
        }

        AffineTransformation getInverse() const {
            AffineTransformation Ret;
            double det = (*this)(0, 0) * (*this)(1, 1) - (*this)(1, 0) * (*this)(0, 1);

            Ret(0, 0) = (*this)(1, 1) / det;
            Ret(1, 0) = -(*this)(1, 0) / det;
            Ret(0, 1) = -(*this)(0, 1) / det;
            Ret(1, 1) = (*this)(0, 0) / det;

            Ret(0, 2) = - Ret(0, 0) * (*this)(0, 2) - Ret(0, 1) * (*this)(1, 2);
            Ret(1, 2) = - Ret(1, 0) * (*this)(0, 2) - Ret(1, 1) * (*this)(1, 2);
            Ret(2, 2) = 1.0;

            return Ret;
        }

        Vector_t<double, 3> getOrigin() const {
            return Vector_t<double, 3>(-(*this)(0, 2), -(*this)(1, 2), 0.0);
        }

        double getAngle() const {
            return atan2((*this)(1, 0), (*this)(0, 0));
        }

        Vector_t<double, 3> transformTo(const Vector_t<double, 3> &v) const {
            vector_t<3> b(0.0); // Create a 3D vector
            b[0] = v(0);
            b[1] = v(1);
            b[2] = 1.0;
            vector_t<3> w = prod_matrix_vector(data_, b);
            return Vector_t<double, 3>(w[0], w[1], 0.0);
        }

        Vector_t<double, 3> transformFrom(const Vector_t<double, 3> &v) const {
            AffineTransformation inv = getInverse();
            return inv.transformTo(v);
        }

        AffineTransformation mult(const AffineTransformation &B) {
            AffineTransformation Ret;
            Ret.data_ = prod(data_, B.data_);
            return Ret;
        }
    };
}

#endif
