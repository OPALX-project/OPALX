#include "Algorithms/Vektor.h"
#include "AppTypes/Tenzor.h"

#include <string>
#include <iostream>
#include <fstream>
#include <list>

namespace mslang {
    typedef std::string::iterator iterator;

    struct BoundingBox {
        Vector_t center_m;
        double width_m;
        double height_m;

        BoundingBox():
            center_m(0.0),
            width_m(0.0),
            height_m(0.0)
        { }

        BoundingBox(const BoundingBox &right):
            center_m(right.center_m),
            width_m(right.width_m),
            height_m(right.height_m)
        { }

        BoundingBox(const Vector_t &llc,
                    const Vector_t &urc):
            center_m(0.5 * (llc + urc)),
            width_m(urc[0] - llc[0]),
            height_m(urc[1] - llc[1])
        { }

        bool isInside(const Vector_t &X) const {
            if (2 * std::abs(X[0] - center_m[0]) <= width_m &&
                2 * std::abs(X[1] - center_m[1]) <= height_m)
                return true;

            return false;
        }

        bool isInside(const BoundingBox &b) const {
            return (isInside(b.center_m + 0.5 * Vector_t( b.width_m,  b.height_m, 0.0)) &&
                    isInside(b.center_m + 0.5 * Vector_t(-b.width_m,  b.height_m, 0.0)) &&
                    isInside(b.center_m + 0.5 * Vector_t(-b.width_m, -b.height_m, 0.0)) &&
                    isInside(b.center_m + 0.5 * Vector_t( b.width_m, -b.height_m, 0.0)));
        }

        virtual void writeGnuplot(std::ofstream &out) const {
            std::vector<Vector_t> pts({Vector_t(center_m[0] + 0.5 * width_m, center_m[1] + 0.5 * height_m, 0),
                        Vector_t(center_m[0] - 0.5 * width_m, center_m[1] + 0.5 * height_m, 0),
                        Vector_t(center_m[0] - 0.5 * width_m, center_m[1] - 0.5 * height_m, 0),
                        Vector_t(center_m[0] + 0.5 * width_m, center_m[1] - 0.5 * height_m, 0)});
            unsigned int width = out.precision() + 8;
            for (unsigned int i = 0; i < 5; ++ i) {
                Vector_t & pt = pts[i % 4];

                out << std::setw(width) << pt[0]
                    << std::setw(width) << pt[1]
                    << std::endl;
            }
            out << std::endl;
        }

        void print(std::ostream &out) const {
            out << std::setw(18) << center_m[0] - 0.5 * width_m
                << std::setw(18) << center_m[1] - 0.5 * height_m
                << std::setw(18) << center_m[0] + 0.5 * width_m
                << std::setw(18) << center_m[1] + 0.5 * height_m
                << std::endl;
        }
    };

    std::ostream & operator<< (std::ostream &out, const BoundingBox &bb);

    struct AffineTransformation: public Tenzor<double, 3> {
        AffineTransformation(const Vector_t& row0,
                             const Vector_t& row1):
            Tenzor(row0[0], row0[1], row0[2], row1[0], row1[1], row1[2], 0.0, 0.0, 1.0) {
        }

        AffineTransformation():
            Tenzor(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0) { }

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

        Vector_t getOrigin() const {
            return Vector_t(-(*this)(0, 2), -(*this)(1, 2), 0.0);
        }

        double getAngle() const {
            return atan2((*this)(1, 0), (*this)(0, 0));
        }

        Vector_t transformTo(const Vector_t &v) const {
            const Tenzor<double, 3> &A = *static_cast<const Tenzor<double, 3>* >(this);
            Vector_t b(v[0], v[1], 1.0);
            Vector_t w = dot(A, b);

            return Vector_t(w[0], w[1], 0.0);
        }

        Vector_t transformFrom(const Vector_t &v) const {
            AffineTransformation inv = getInverse();
            return inv.transformTo(v);
        }

        AffineTransformation mult(const AffineTransformation &B) {
            AffineTransformation Ret;
            const Tenzor<double, 3> &A = *static_cast<const Tenzor<double, 3> *>(this);
            const Tenzor<double, 3> &BTenz = *static_cast<const Tenzor<double, 3> *>(&B);
            Tenzor<double, 3> &C = *static_cast<Tenzor<double, 3> *>(&Ret);

            C = dot(A, BTenz);

            return Ret;
        }
    };

    struct Base;

    struct Function {
        virtual ~Function() {};

        virtual void print(int indent) = 0;
        virtual void apply(std::vector<Base*> &bfuncs) = 0;
    };

    struct Base: public Function {
        AffineTransformation trafo_m;
        BoundingBox bb_m;

        Base():
            trafo_m()
        { }

        Base(const Base &right):
            trafo_m(right.trafo_m),
            bb_m(right.bb_m)
        { }

        virtual Base* clone() const = 0;
        virtual void writeGnuplot(std::ofstream &out) const = 0;
        virtual void computeBoundingBox() = 0;
        virtual bool isInside(const Vector_t &R) const = 0;
    };


    struct QuadTree {
        int level_m;
        std::list<Base*> objects_m;
        BoundingBox bb_m;
        QuadTree *nodes_m;

        QuadTree():
            level_m(0),
            bb_m(),
            nodes_m(0)
        { }

        QuadTree(int l, const BoundingBox &b):
            level_m(l),
            bb_m(b),
            nodes_m(0)
        { }

        QuadTree(const QuadTree &right);

        ~QuadTree();

        void operator=(const QuadTree &right);

        void transferIfInside(std::list<Base*> &objs);
        void buildUp();


        void writeGnuplot(std::ofstream &out) const {
            out << "# level: " << level_m << ", size: " << objects_m.size() << std::endl;
            bb_m.writeGnuplot(out);
            out << "# num holes: " << objects_m.size() << std::endl;
            for (const Base *obj: objects_m) {
                obj->writeGnuplot(out);
            }
            out << std::endl;

            if (nodes_m != 0) {
                for (unsigned int i = 0; i < 4u; ++ i) {
                    nodes_m[i].writeGnuplot(out);
                }
            }
        }

        bool isInside(const Vector_t &R) const;
    };

    bool parse(std::string str, Function* &fun);
}