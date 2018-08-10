#include "Utilities/MSLang.h"
#include "Utilities/PortableBitmapReader.h"
#include "Utilities/Mesher.h"
#include "Algorithms/Quaternion.h"
#include "Physics/Physics.h"

#include <boost/regex.hpp>
#include <boost/filesystem.hpp>

#include <iostream>
#include <string>
#include <fstream>
#include <streambuf>
#include <cstdlib>
#include <cmath>

namespace mslang {
    std::string UDouble("([0-9]+\\.?[0-9]*([Ee][+-]?[0-9]+)?)");
    std::string Double("(-?[0-9]+\\.?[0-9]*([Ee][+-]?[0-9]+)?)");
    std::string UInt("([0-9]+)");
    std::string FCall("([a-z_]*)\\((.*)");

    bool parse(iterator &it, const iterator &end, Function* &fun);

    std::ostream & operator<< (std::ostream &out, const BoundingBox &bb) {
        bb.print(out);
        return out;
    }

    struct Rectangle: public Base {
        double width_m;
        double height_m;

        Rectangle():
            Base(),
            width_m(0.0),
            height_m(0.0)
        { }

        Rectangle(const Rectangle &right):
            Base(right),
            width_m(right.width_m),
            height_m(right.height_m)
        { }

        virtual ~Rectangle() { }

        virtual void print(int indentwidth) {
            std::string indent(indentwidth, ' ');
            std::string indent2(indentwidth + 8, ' ');
            Vector_t origin = trafo_m.getOrigin();
            double angle = trafo_m.getAngle() * Physics::rad2deg;
            std::cout << indent << "rectangle, \n"
                      << indent2 << "w: " << width_m << ", \n"
                      << indent2 << "h: " << height_m << ", \n"
                      << indent2 << "origin: " << origin[0] << ", " << origin[1] << ",\n"
                      << indent2 << "angle: " << angle << "\n"
                      << indent2 << trafo_m(0, 0) << "\t" << trafo_m(0, 1) << "\t" << trafo_m(0, 2) << "\n"
                      << indent2 << trafo_m(1, 0) << "\t" << trafo_m(1, 1) << "\t" << trafo_m(1, 2) << "\n"
                      << indent2 << trafo_m(2, 0) << "\t" << trafo_m(2, 1) << "\t" << trafo_m(2, 2) << std::endl;
        }

        virtual void computeBoundingBox() {
            std::vector<Vector_t> corners({Vector_t(0.5 * width_m, 0.5 * height_m, 0),
                        Vector_t(-0.5 * width_m, 0.5 * height_m, 0),
                        Vector_t(-0.5 * width_m, -0.5 * height_m, 0),
                        Vector_t(0.5 * width_m, -0.5 * height_m, 0)});

            for (Vector_t &v: corners) {
                v = trafo_m.transformFrom(v);
            }

            Vector_t llc = corners[0], urc = corners[0];
            for (unsigned int i = 1; i < 4; ++ i) {
                if (corners[i][0] < llc[0]) llc[0] = corners[i][0];
                else if (corners[i][0] > urc[0]) urc[0] = corners[i][0];

                if (corners[i][1] < llc[1]) llc[1] = corners[i][1];
                else if (corners[i][1] > urc[1]) urc[1] = corners[i][1];
            }

            bb_m = BoundingBox(llc, urc);

            for (auto item: divisor_m) {
                item->computeBoundingBox();
            }
        }

        virtual bool isInside(const Vector_t &R) const {
            if (!bb_m.isInside(R)) return false;

            Vector_t X = trafo_m.transformTo(R);
            if (2 * std::abs(X[0]) <= width_m &&
                2 * std::abs(X[1]) <= height_m) {
                for (auto item: divisor_m) {
                    if (item->isInside(R))
                        return false;
                }
                return true;
            }

            return false;
        }

        virtual void writeGnuplot(std::ofstream &out) const {
            std::vector<Vector_t> pts({Vector_t(0.5 * width_m, 0.5 * height_m, 0),
                        Vector_t(-0.5 * width_m, 0.5 * height_m, 0),
                        Vector_t(-0.5 * width_m, -0.5 * height_m, 0),
                        Vector_t(0.5 * width_m, -0.5 * height_m, 0)});
            unsigned int width = out.precision() + 8;
            for (unsigned int i = 0; i < 5; ++ i) {
                Vector_t pt = pts[i % 4];
                pt = trafo_m.transformFrom(pt);

                out << std::setw(width) << pt[0]
                    << std::setw(width) << pt[1]
                    << std::endl;
            }
            out << std::endl;

            for (auto item: divisor_m) {
                item->writeGnuplot(out);
            }

            // bb_m.writeGnuplot(out);
        }

        virtual void apply(std::vector<Base*> &bfuncs) {
            bfuncs.push_back(this->clone());
        }

        virtual Base* clone() const {
            Rectangle *rect = new Rectangle;
            rect->width_m = width_m;
            rect->height_m = height_m;
            rect->trafo_m = trafo_m;

            for (auto item: divisor_m) {
                rect->divisor_m.push_back(item->clone());
            }
            return rect;
        }

        virtual void divideBy(std::vector<Base*> &divisors) {
            for (auto item: divisors) {
                if (bb_m.doesIntersect(item->bb_m)) {
                    divisor_m.push_back(item->clone());
                }
            }
        }

        static
        bool parse_detail(iterator &it, const iterator &end, Function* fun) {
            std::string str(it, end);
            boost::regex argumentList(UDouble + "," + UDouble + "(\\).*)");
            boost::smatch what;

            if (!boost::regex_match(str, what, argumentList)) return false;

            Rectangle *rect = static_cast<Rectangle*>(fun);
            rect->width_m = atof(std::string(what[1]).c_str());
            rect->height_m = atof(std::string(what[3]).c_str());

            std::string fullMatch = what[0];
            std::string rest = what[5];
            it += (fullMatch.size() - rest.size() + 1);

            return true;

        }
    };

    struct Ellipse: public Base {
        double width_m;
        double height_m;

        Ellipse():
            Base(),
            width_m(0.0),
            height_m(0.0)
        { }

        Ellipse(const Ellipse &right):
            Base(right),
            width_m(right.width_m),
            height_m(right.height_m)
        { }

        virtual ~Ellipse() { }

        virtual void print(int indentwidth) {
            std::string indent(indentwidth, ' ');
            std::string indent2(indentwidth + 8, ' ');
            Vector_t origin = trafo_m.getOrigin();
            double angle = trafo_m.getAngle() * Physics::rad2deg;
            std::cout << indent << "ellipse, \n"
                      << indent2 << "w: " << width_m << ", \n"
                      << indent2 << "h: " << height_m << ", \n"
                      << indent2 << "origin: " << origin[0] << ", " << origin[1] << ",\n"
                      << indent2 << "angle: " << angle << "\n"
                      << indent2 << std::setw(14) << trafo_m(0, 0) << std::setw(14) << trafo_m(0, 1) << std::setw(14) << trafo_m(0, 2) << "\n"
                      << indent2 << std::setw(14) << trafo_m(1, 0) << std::setw(14) << trafo_m(1, 1) << std::setw(14) << trafo_m(1, 2) << "\n"
                      << indent2 << std::setw(14) << trafo_m(2, 0) << std::setw(14) << trafo_m(2, 1) << std::setw(14) << trafo_m(2, 2)
                      << std::endl;
        }

        virtual void writeGnuplot(std::ofstream &out) const {
            const unsigned int N = 101;
            const double dp = Physics::two_pi / (N - 1);
            const unsigned int colwidth = out.precision() + 8;

            double phi = 0;
            for (unsigned int i = 0; i < N; ++ i, phi += dp) {
                Vector_t pt(0.0);
                pt[0] = std::copysign(sqrt(std::pow(height_m * width_m * 0.25, 2) /
                                           (std::pow(height_m * 0.5, 2) +
                                            std::pow(width_m * 0.5 * tan(phi), 2))),
                                      cos(phi));
                pt[1] = pt[0] * tan(phi);
                pt = trafo_m.transformFrom(pt);

                out << std::setw(colwidth) << pt[0]
                    << std::setw(colwidth) << pt[1]
                    << std::endl;
            }
            out << std::endl;

            for (auto item: divisor_m) {
                item->writeGnuplot(out);
            }

            // bb_m.writeGnuplot(out);
        }

        virtual void apply(std::vector<Base*> &bfuncs) {
            bfuncs.push_back(this->clone());
        }

        virtual Base* clone() const{
            Ellipse *elps = new Ellipse;
            elps->width_m = width_m;
            elps->height_m = height_m;
            elps->trafo_m = trafo_m;

            for (auto item: divisor_m) {
                elps->divisor_m.push_back(item->clone());
            }

            return elps;
        }

        virtual void computeBoundingBox() {
            Vector_t llc(0.0), urc(0.0);
            const Vector_t e_x(1.0, 0.0, 0.0), e_y(0.0, 1.0, 0.0);
            const Vector_t center = trafo_m.transformFrom(Vector_t(0.0));
            const Vector_t e_xp = trafo_m.transformFrom(e_x) - center;
            const Vector_t e_yp = trafo_m.transformFrom(e_y) - center;
            const double &M11 = e_xp[0];
            const double &M12 = e_yp[0];
            const double &M21 = e_xp[1];
            const double &M22 = e_yp[1];

            double t = atan2(height_m * M12, width_m * M11);
            double halfwidth = 0.5 * (M11 * width_m * cos(t) +
                                      M12 * height_m * sin(t));
            llc[0] = center[0] - std::abs(halfwidth);
            urc[0] = center[0] + std::abs(halfwidth);

            t = atan2(height_m * M22, width_m * M21);

            double halfheight = 0.5 * (M21 * width_m * cos(t) +
                                       M22 * height_m * sin(t));

            llc[1] = center[1] - std::abs(halfheight);
            urc[1] = center[1] + std::abs(halfheight);

            bb_m = BoundingBox(llc, urc);

            for (auto item: divisor_m) {
                item->computeBoundingBox();
            }
        }

        virtual bool isInside(const Vector_t &R) const {
            if (!bb_m.isInside(R)) return false;

            Vector_t X = trafo_m.transformTo(R);
            if (4 * (std::pow(X[0] / width_m, 2) + std::pow(X[1] / height_m, 2)) <= 1) {

                for (auto item: divisor_m) {
                    if (item->isInside(R)) return false;
                }

                return true;
            }

            return false;
        }

        virtual void divideBy(std::vector<Base*> &divisors) {
            for (auto item: divisors) {
                if (bb_m.doesIntersect(item->bb_m)) {
                    divisor_m.push_back(item->clone());
                }
            }
        }

        static
        bool parse_detail(iterator &it, const iterator &end, Function* fun) {
            std::string str(it, end);
            boost::regex argumentList(UDouble + "," + UDouble + "(\\).*)");
            boost::smatch what;

            if (!boost::regex_match(str, what, argumentList)) return false;

            Ellipse *elps = static_cast<Ellipse*>(fun);
            elps->width_m = atof(std::string(what[1]).c_str());
            elps->height_m = atof(std::string(what[3]).c_str());

            std::string fullMatch = what[0];
            std::string rest = what[5];
            it += (fullMatch.size() - rest.size() + 1);

            return true;
        }
    };

    void Triangle::print(int indentwidth) {
        std::string indent(indentwidth, ' ');
        std::string indentbase(4, ' ');
        Vector_t origin = trafo_m.getOrigin();
        double angle = trafo_m.getAngle() * Physics::rad2deg;

        std::cout << indent << "triangle, \n";

        for (unsigned int i = 0; i < 3u; ++ i) {
            std::cout << indent << indentbase << "node " << i << ": " << nodes_m[i] << "\n";
        }
        std::cout << indent << indentbase << "origin: " << origin[0] << ", " << origin[1] << ",\n"
                  << indent << indentbase << "angle: " << angle << "\n"
                  << indent << indentbase << trafo_m(0, 0) << "\t" << trafo_m(0, 1) << "\t" << trafo_m(0, 2) << "\n"
                  << indent << indentbase << trafo_m(1, 0) << "\t" << trafo_m(1, 1) << "\t" << trafo_m(1, 2) << "\n"
                  << indent << indentbase << trafo_m(2, 0) << "\t" << trafo_m(2, 1) << "\t" << trafo_m(2, 2) << std::endl;
    }
    void Triangle::apply(std::vector<Base*> &bfuncs) {
        bfuncs.push_back(this->clone());
    }

    Base* Triangle::clone() const {
        Triangle *tri = new Triangle(*this);
        tri->trafo_m = trafo_m;

        for (auto item: divisor_m) {
            tri->divisor_m.push_back(item->clone());
        }

        return tri;
    }

    void Triangle::writeGnuplot(std::ofstream &out) const {
        unsigned int width = out.precision() + 8;

        for (unsigned int i = 0; i < 4u; ++ i) {
            Vector_t corner = trafo_m.transformFrom(nodes_m[i % 3u]);

            out << std::setw(width) << corner[0]
                << std::setw(width) << corner[1]
                << std::endl;
        }
        out << std::endl;

        bb_m.writeGnuplot(out);

    }

    void Triangle::computeBoundingBox() {
        std::vector<Vector_t> corners;
        for (unsigned int i = 0; i < 3u; ++ i) {
            corners.push_back(trafo_m.transformFrom(nodes_m[i]));
        }

        Vector_t llc = corners[0], urc = corners[0];
        for (unsigned int i = 1u; i < 3u; ++ i) {
            if (corners[i][0] < llc[0]) llc[0] = corners[i][0];
            else if (corners[i][0] > urc[0]) urc[0] = corners[i][0];

            if (corners[i][1] < llc[1]) llc[1] = corners[i][1];
            else if (corners[i][1] > urc[1]) urc[1] = corners[i][1];
        }

        bb_m = BoundingBox(llc, urc);

        for (auto item: divisor_m) {
            item->computeBoundingBox();
        }
    }

    double Triangle::crossProduct(const Vector_t &pt, unsigned int nodeNum) const {
        nodeNum = nodeNum % 3u;
        unsigned int nextNode = (nodeNum + 1) % 3u;
        Vector_t nodeToPt = pt - nodes_m[nodeNum];
        Vector_t nodeToNext = nodes_m[nextNode] - nodes_m[nodeNum];

        return nodeToPt[0] * nodeToNext[1] - nodeToPt[1] * nodeToNext[0];

    }

    bool Triangle::isInside(const Vector_t &R) const {
        Vector_t X = trafo_m.transformTo(R);

        bool test0 = (crossProduct(X, 0) <= 0.0);
        bool test1 = (crossProduct(X, 1) <= 0.0);
        bool test2 = (crossProduct(X, 2) <= 0.0);

        if (!(test0 && test1 && test2)) return false;

        for (auto item: divisor_m)
            if (item->isInside(R)) return false;

        return true;
    }

    void Triangle::orientNodesCCW() {
        if (crossProduct(nodes_m[0], 1) > 0.0) {
            std::swap(nodes_m[1], nodes_m[2]);
        }
    }

    void Triangle::divideBy(std::vector<Base*> &divisors) {
        for (auto item: divisors) {
            if (bb_m.doesIntersect(item->bb_m)) {
                divisor_m.push_back(item->clone());
            }
        }
    }

    struct Repeat: public Function {
        Function* func_m;
        unsigned int N_m;
        double shiftx_m;
        double shifty_m;
        double rot_m;

        virtual ~Repeat() {
            delete func_m;
        }

        virtual void print(int indentwidth) {
            std::string indent(indentwidth, ' ');
            std::string indent2(indentwidth + 8, ' ');
            std::cout << indent << "repeat, " << std::endl;
            func_m->print(indentwidth + 8);
            std::cout << ",\n"
                      << indent2 << "N: " << N_m << ", \n"
                      << indent2 << "dx: " << shiftx_m << ", \n"
                      << indent2 << "dy: " << shifty_m;
        }

        virtual void apply(std::vector<Base*> &bfuncs) {
            AffineTransformation trafo(Vector_t(cos(rot_m), sin(rot_m), -shiftx_m),
                                       Vector_t(-sin(rot_m), cos(rot_m), -shifty_m));

            func_m->apply(bfuncs);
            const unsigned int size = bfuncs.size();

            AffineTransformation current_trafo = trafo;
            for (unsigned int i = 0; i < N_m; ++ i) {
                for (unsigned int j = 0; j < size; ++ j) {
                    Base *obj = bfuncs[j]->clone();
                    obj->trafo_m = obj->trafo_m.mult(current_trafo);
                    bfuncs.push_back(obj);
                }

                current_trafo = current_trafo.mult(trafo);
            }
        }

        static
        bool parse_detail(iterator &it, const iterator &end, Function* &fun) {
            Repeat *rep = static_cast<Repeat*>(fun);
            if (!parse(it, end, rep->func_m)) return false;

            boost::regex argumentListTrans("," + UInt + "," + Double + "," + Double + "\\)(.*)");
            boost::regex argumentListRot("," + UInt + "," + Double + "\\)(.*)");
            boost::smatch what;

            std::string str(it, end);
            if (boost::regex_match(str, what, argumentListTrans)) {
                rep->N_m = atof(std::string(what[1]).c_str());
                rep->shiftx_m = atof(std::string(what[2]).c_str());
                rep->shifty_m = atof(std::string(what[4]).c_str());
                rep->rot_m = 0.0;

                std::string fullMatch = what[0];
                std::string rest = what[6];

                it += (fullMatch.size() - rest.size());

                return true;
            }

            if (boost::regex_match(str, what, argumentListRot)) {
                rep->N_m = atof(std::string(what[1]).c_str());
                rep->shiftx_m = 0.0;
                rep->shifty_m = 0.0;
                rep->rot_m = atof(std::string(what[2]).c_str());

                std::string fullMatch = what[0];
                std::string rest = what[4];

                it += (fullMatch.size() - rest.size());

                return true;
            }

            return false;
        }
    };

    struct Translation: public Function {
        Function* func_m;
        double shiftx_m;
        double shifty_m;

        virtual ~Translation() {
            delete func_m;
        }

        virtual void print(int indentwidth) {
            std::string indent(indentwidth, ' ');
            std::string indent2(indentwidth + 8, ' ');
            std::cout << indent << "translate, " << std::endl;
            func_m->print(indentwidth + 8);
            std::cout << ",\n"
                      << indent2 << "dx: " << shiftx_m << ", \n"
                      << indent2 << "dy: " << shifty_m;
        }

        void applyTranslation(std::vector<Base*> &bfuncs) {
            AffineTransformation shift(Vector_t(1.0, 0.0, -shiftx_m),
                                       Vector_t(0.0, 1.0, -shifty_m));

            const unsigned int size = bfuncs.size();
            for (unsigned int j = 0; j < size; ++ j) {
                Base *obj = bfuncs[j];
                obj->trafo_m = obj->trafo_m.mult(shift);

                if (obj->divisor_m.size() > 0)
                    applyTranslation(obj->divisor_m);
            }
        }

        virtual void apply(std::vector<Base*> &bfuncs) {
            func_m->apply(bfuncs);
            applyTranslation(bfuncs);
        }

        static
        bool parse_detail(iterator &it, const iterator &end, Function* &fun) {
            Translation *trans = static_cast<Translation*>(fun);
            if (!parse(it, end, trans->func_m)) return false;

            boost::regex argumentList("," + Double + "," + Double + "\\)(.*)");
            boost::smatch what;

            std::string str(it, end);
            if (!boost::regex_match(str, what, argumentList)) return false;

            trans->shiftx_m = atof(std::string(what[1]).c_str());
            trans->shifty_m = atof(std::string(what[3]).c_str());

            std::string fullMatch = what[0];
            std::string rest = what[5];

            it += (fullMatch.size() - rest.size());

            return true;
        }
    };


    struct Rotation: public Function {
        Function* func_m;
        double angle_m;

        virtual ~Rotation() {
            delete func_m;
        }

        virtual void print(int indentwidth) {
            std::string indent(indentwidth, ' ');
            std::string indent2(indentwidth + 8, ' ');
            std::cout << indent << "rotate, " << std::endl;
            func_m->print(indentwidth + 8);
            std::cout << ",\n"
                      << indent2 << "angle: " << angle_m;
        }

        void applyRotation(std::vector<Base*> &bfuncs) {

            AffineTransformation rotation(Vector_t(cos(angle_m), sin(angle_m), 0.0),
                                          Vector_t(-sin(angle_m), cos(angle_m), 0.0));

            const unsigned int size = bfuncs.size();

            for (unsigned int j = 0; j < size; ++ j) {
                Base *obj = bfuncs[j];
                obj->trafo_m = obj->trafo_m.mult(rotation);

                if (obj->divisor_m.size() > 0)
                    applyRotation(obj->divisor_m);
            }
        }

        virtual void apply(std::vector<Base*> &bfuncs) {
            func_m->apply(bfuncs);
            applyRotation(bfuncs);
        }

        static
        bool parse_detail(iterator &it, const iterator &end, Function* &fun) {
            Rotation *rot = static_cast<Rotation*>(fun);
            if (!parse(it, end, rot->func_m)) return false;

            boost::regex argumentList("," + Double + "\\)(.*)");
            boost::smatch what;

            std::string str(it, end);
            if (!boost::regex_match(str, what, argumentList)) return false;

            rot->angle_m = atof(std::string(what[1]).c_str());

            std::string fullMatch = what[0];
            std::string rest = what[3];

            it += (fullMatch.size() - rest.size());

            return true;
        }
    };

    struct Shear: public Function {
        Function* func_m;
        double angleX_m;
        double angleY_m;

        virtual ~Shear() {
            delete func_m;
        }

        virtual void print(int indentwidth) {
            std::string indent(indentwidth, ' ');
            std::string indent2(indentwidth + 8, ' ');
            std::cout << indent << "shear, " << std::endl;
            func_m->print(indentwidth + 8);
            if (std::abs(angleX_m) > 0.0) {
                std::cout << ",\n"
                          << indent2 << "angle X: " << angleX_m;
            } else {
                std::cout << ",\n"
                          << indent2 << "angle Y: " << angleY_m;
            }
        }

        void applyShear(std::vector<Base*> &bfuncs) {
            AffineTransformation shear(Vector_t(1.0, tan(angleX_m), 0.0),
                                       Vector_t(-tan(angleY_m), 1.0, 0.0));

            const unsigned int size = bfuncs.size();

            for (unsigned int j = 0; j < size; ++ j) {
                Base *obj = bfuncs[j];
                obj->trafo_m = obj->trafo_m.mult(shear);

                if (obj->divisor_m.size() > 0)
                    applyShear(obj->divisor_m);
            }
        }

        virtual void apply(std::vector<Base*> &bfuncs) {
            func_m->apply(bfuncs);
            applyShear(bfuncs);
        }

        static
        bool parse_detail(iterator &it, const iterator &end, Function* &fun) {
            Shear *shr = static_cast<Shear*>(fun);
            if (!parse(it, end, shr->func_m)) return false;

            boost::regex argumentList("," + Double + "," + Double + "\\)(.*)");
            boost::smatch what;

            std::string str(it, end);
            if (!boost::regex_match(str, what, argumentList)) return false;

            shr->angleX_m = atof(std::string(what[1]).c_str());
            shr->angleY_m = atof(std::string(what[3]).c_str());

            if (std::abs(shr->angleX_m) > 0.0 && std::abs(shr->angleY_m) > 0.0)
                return false;

            std::string fullMatch = what[0];
            std::string rest = what[5];

            it += (fullMatch.size() - rest.size());

            return true;
        }
    };

    struct Union: public Function {
        std::vector<Function*> funcs_m;

        virtual ~Union () {
            for (Function* func: funcs_m) {
                delete func;
            }
        }

        virtual void print(int indentwidth) {
            std::string indent(indentwidth, ' ');
            std::string indent2(indentwidth + 8, ' ');
            std::string indent3(indentwidth + 16, ' ');
            std::cout << indent << "union, " << std::endl;
            std::cout << indent2 << "funcs: {\n";
            funcs_m.front()->print(indentwidth + 16);
            for (unsigned int i = 1; i < funcs_m.size(); ++ i) {
                std::cout << "\n"
                          << indent3 << "," << std::endl;
                funcs_m[i]->print(indentwidth + 16);
            }
            std::cout << "\n"
                      << indent2 << "} ";
        }

        virtual void apply(std::vector<Base*> &bfuncs) {
            for (unsigned int i = 0; i < funcs_m.size(); ++ i) {
                std::vector<Base*> children;
                Function *func = funcs_m[i];
                func->apply(children);
                bfuncs.insert(bfuncs.end(), children.begin(), children.end());
            }
        }

        static
        bool parse_detail(iterator &it, const iterator &end, Function* &fun) {
            Union *unin = static_cast<Union*>(fun);
            unin->funcs_m.push_back(NULL);
            if (!parse(it, end, unin->funcs_m.back())) return false;

            boost::regex argumentList("(,[a-z]+\\(.*)");
            boost::regex endParenthesis("\\)(.*)");
            boost::smatch what;

            std::string str(it, end);
            while (boost::regex_match(str, what, argumentList)) {
                iterator it2 = it + 1;
                unin->funcs_m.push_back(NULL);

                if (!parse(it2, end, unin->funcs_m.back())) return false;

                it = it2;
                str = std::string(it, end);
            }

            str = std::string(it, end);
            if (!boost::regex_match(str, what, endParenthesis)) return false;

            std::string fullMatch = what[0];
            std::string rest = what[1];

            it += (fullMatch.size() - rest.size());

            return true;
        }
    };

    struct Mask: public Function {
        static
        bool parse_detail(iterator &it, const iterator &end, Function* &fun) {
            Mask *pixmap = static_cast<Mask*>(fun);

            std::string str(it, end);
            boost::regex argument("'([^,]+)'," + UDouble + "," + UDouble + "(\\).*)");
            boost::smatch what;
            if (!boost::regex_match(str, what, argument)) return false;

            std::string filename = what[1];
            if (!boost::filesystem::exists(filename)) {
                ERRORMSG("file '" << filename << "' doesn't exists" << endl);
                return false;
            }

            PortableBitmapReader reader(filename);
            unsigned int width = reader.getWidth();
            unsigned int height = reader.getHeight();
            double pixel_width = atof(std::string(what[2]).c_str()) / width;
            double pixel_height = atof(std::string(what[4]).c_str()) / height;

            for (unsigned int i = 0; i < height; ++ i) {
                for (unsigned int j = 0; j < width; ++ j) {
                    if (reader.isBlack(i, j)) {
                        Rectangle rect;
                        rect.width_m = pixel_width;
                        rect.height_m = pixel_height;
                        rect.trafo_m = AffineTransformation(Vector_t(1, 0, (0.5 * width - j) * pixel_width),
                                                            Vector_t(0, 1, (i - 0.5 * height) * pixel_height));

                        pixmap->pixels_m.push_back(rect);
                        pixmap->pixels_m.back().computeBoundingBox();
                    }
                }
            }

            std::string fullMatch = what[0];
            std::string rest = what[6];
            it += (fullMatch.size() - rest.size() + 1);

            return true;
        }

        virtual void print(int ident) {
            for (auto pix: pixels_m) pix.print(ident);
        }

        virtual void apply(std::vector<Base*> &bfuncs) {
            for (auto pix: pixels_m) pix.apply(bfuncs);
        }

        std::vector<Rectangle> pixels_m;
    };

    struct Difference: public Function {
        Function* dividend_m;
        Function* divisor_m;

        Difference()
        { }

        Difference(const Difference &right):
            dividend_m(right.dividend_m),
            divisor_m(right.divisor_m)
        { }

        virtual ~Difference() {
            delete dividend_m;
            delete divisor_m;
        }

        virtual void print(int indentwidth) {
            std::string indent(' ', indentwidth);
            std::cout << indent << "Difference\n"
                      << indent << "    nominators\n";
            dividend_m->print(indentwidth + 8);

            std::cout << indent << "    denominators\n";
            divisor_m->print(indentwidth + 8);
        }

        virtual void apply(std::vector<Base*> &bfuncs) {
            std::vector<Base*> nom, denom;

            dividend_m->apply(nom);
            divisor_m->apply(denom);
            for (auto item: nom) {
                item->divideBy(denom);
                bfuncs.push_back(item->clone());
            }
        }

        static
        bool parse_detail(iterator &it, const iterator &end, Function* &fun) {
            Difference *dif = static_cast<Difference*>(fun);
            if (!parse(it, end, dif->dividend_m)) return false;

            boost::regex argumentList("(,[a-z]+\\(.*)");
            boost::regex endParenthesis("\\)(.*)");
            boost::smatch what;

            std::string str(it, end);
            if (!boost::regex_match(str, what, argumentList)) return false;

            iterator it2 = it + 1;
            if (!parse(it2, end, dif->divisor_m)) return false;

            it = it2;
            str = std::string(it, end);
            if (!boost::regex_match(str, what, endParenthesis)) return false;

            std::string fullMatch = what[0];
            std::string rest = what[1];

            it += (fullMatch.size() - rest.size());

            return true;
        }
    };

    struct SymmetricDifference: public Function {
        Function *firstOperand_m;
        Function *secondOperand_m;

        SymmetricDifference()
        { }

        SymmetricDifference(const SymmetricDifference &right):
            firstOperand_m(right.firstOperand_m),
            secondOperand_m(right.secondOperand_m)
        { }

        virtual ~SymmetricDifference() {
            delete firstOperand_m;
            delete secondOperand_m;
        }

        virtual void print(int indentwidth) {
            std::string indent(' ', indentwidth);
            std::cout << indent << "Symmetric division\n"
                      << indent << "    first operand\n";
            firstOperand_m->print(indentwidth + 8);

            std::cout << indent << "    second operand\n";
            secondOperand_m->print(indentwidth + 8);
        }

        virtual void apply(std::vector<Base*> &bfuncs) {
            std::vector<Base*> first, second;

            firstOperand_m->apply(first);
            secondOperand_m->apply(second);
            for (auto item: first) {
                item->divideBy(second);
                bfuncs.push_back(item->clone());
            }

            for (auto item: second)
                delete item;
            for (auto item: first)
                delete item;

            first.clear();
            second.clear();

            secondOperand_m->apply(first);
            firstOperand_m->apply(second);
            for (auto item: first) {
                item->divideBy(second);
                bfuncs.push_back(item->clone());
            }

            first.clear();
            for (auto item: second)
                delete item;
            second.clear();
        }

        static
        bool parse_detail(iterator &it, const iterator &end, Function* &fun) {
            SymmetricDifference *dif = static_cast<SymmetricDifference*>(fun);
            if (!parse(it, end, dif->firstOperand_m)) return false;

            boost::regex argumentList("(,[a-z]+\\(.*)");
            boost::regex endParenthesis("\\)(.*)");
            boost::smatch what;

            std::string str(it, end);
            if (!boost::regex_match(str, what, argumentList)) return false;

            iterator it2 = it + 1;
            if (!parse(it2, end, dif->secondOperand_m)) return false;

            it = it2;
            str = std::string(it, end);
            if (!boost::regex_match(str, what, endParenthesis)) return false;

            std::string fullMatch = what[0];
            std::string rest = what[1];

            it += (fullMatch.size() - rest.size());

            return true;
        }
    };

    struct Intersection: public Function {
        Function *firstOperand_m;
        Function *secondOperand_m;

        Intersection()
        { }

        Intersection(const Intersection &right):
            firstOperand_m(right.firstOperand_m),
            secondOperand_m(right.secondOperand_m)
        { }

        virtual ~Intersection() {
            delete firstOperand_m;
            delete secondOperand_m;
        }

        virtual void print(int indentwidth) {
            std::string indent(' ', indentwidth);
            std::cout << indent << "Intersection\n"
                      << indent << "    first operand\n";
            firstOperand_m->print(indentwidth + 8);

            std::cout << indent << "    second operand\n";
            secondOperand_m->print(indentwidth + 8);
        }

        virtual void apply(std::vector<Base*> &bfuncs) {
            std::vector<Base*> first, firstrep, second;

            firstOperand_m->apply(first);
            firstOperand_m->apply(firstrep);
            secondOperand_m->apply(second);

            for (auto item: firstrep) {
                item->divideBy(second);
            }

            for (auto item: first) {
                item->divideBy(firstrep);
                bfuncs.push_back(item->clone());
            }

            for (auto item: first)
                delete item;
            for (auto item: firstrep)
                delete item;
            for (auto item: second)
                delete item;
        }

        static
        bool parse_detail(iterator &it, const iterator &end, Function* &fun) {
            Intersection *inter = static_cast<Intersection*>(fun);
            if (!parse(it, end, inter->firstOperand_m)) return false;

            boost::regex argumentList("(,[a-z]+\\(.*)");
            boost::regex endParenthesis("\\)(.*)");
            boost::smatch what;

            std::string str(it, end);
            if (!boost::regex_match(str, what, argumentList)) return false;

            iterator it2 = it + 1;
            if (!parse(it2, end, inter->secondOperand_m)) return false;

            it = it2;
            str = std::string(it, end);
            if (!boost::regex_match(str, what, endParenthesis)) return false;

            std::string fullMatch = what[0];
            std::string rest = what[1];

            it += (fullMatch.size() - rest.size());

            return true;
        }
    };

    struct Polygon: public Function {
        std::vector<Triangle> triangles_m;

        void triangulize(std::vector<Vector_t> &nodes) {
            Mesher mesher(nodes);
            triangles_m = mesher.getTriangles();
        }

        static
        bool parse_detail(iterator &it, const iterator &end, Function* &fun) {
            Polygon *poly = static_cast<Polygon*>(fun);

            std::vector<Vector_t> nodes;
            std::string str(it, end);
            boost::regex argument(Double + "," + Double + ";(.*)");
            boost::regex argumentEnd(Double + "," + Double + "(\\).*)");

            boost::smatch what;
            while (boost::regex_match(str, what, argument)) {
                Vector_t p(atof(std::string(what[1]).c_str()),
                           atof(std::string(what[3]).c_str()),
                           1.0);
                nodes.push_back(p);

                std::string fullMatch = what[0];
                std::string rest = what[5];
                it += (fullMatch.size() - rest.size());

                str = std::string(it, end);
            }

            if (!boost::regex_match(str, what, argumentEnd) ||
                nodes.size() < 2) return false;

            Vector_t p(atof(std::string(what[1]).c_str()),
                       atof(std::string(what[3]).c_str()),
                       1.0);
            nodes.push_back(p);


            std::string fullMatch = what[0];
            std::string rest = what[5];
            it += (fullMatch.size() - rest.size() + 1);

            str = std::string(it, end);

            poly->triangulize(nodes);

            return true;
        }

        virtual void print(int ident) {
            // for (auto pix: pixels_m) pix.print(ident);
        }

        virtual void apply(std::vector<Base*> &bfuncs) {
            for (Triangle &tri: triangles_m)
                bfuncs.push_back(tri.clone());
        }
    };

    QuadTree::QuadTree(const QuadTree &right):
        level_m(right.level_m),
        objects_m(right.objects_m.begin(),
                  right.objects_m.end()),
        bb_m(right.bb_m),
        nodes_m(0)
    {
        if (right.nodes_m != 0) {
            nodes_m = new QuadTree[4];
            for (unsigned int i = 0; i < 4u; ++ i) {
                nodes_m[i] = right.nodes_m[i];
            }
        }
    }

    QuadTree::~QuadTree() {
        for (Base *&obj: objects_m)
            obj = 0; // memory isn't handled by QuadTree class

        if (nodes_m != 0) {
            delete[] nodes_m;
        }
        nodes_m = 0;
    }

    void QuadTree::operator=(const QuadTree &right) {
        level_m = right.level_m;
        objects_m.insert(objects_m.end(),
                         right.objects_m.begin(),
                         right.objects_m.end());
        bb_m = right.bb_m;

        if (nodes_m != 0) delete[] nodes_m;
        nodes_m = 0;

        if (right.nodes_m != 0) {
            nodes_m = new QuadTree[4];
            for (unsigned int i = 0; i < 4u; ++ i) {
                nodes_m[i] = right.nodes_m[i];
            }
        }
    }

    void QuadTree::transferIfInside(std::list<Base*> &objs) {
        for (Base* &obj: objs) {
            if (bb_m.isInside(obj->bb_m)) {
                objects_m.push_back(obj);
                obj = 0;
            }
        }

        objs.remove_if([](const Base *obj) { return obj == 0; });
    }

    void QuadTree::buildUp() {
        QuadTree *next = new QuadTree[4];
        next[0] = QuadTree(level_m + 1,
                           BoundingBox(bb_m.center_m,
                                       Vector_t(bb_m.center_m[0] + 0.5 * bb_m.width_m,
                                                bb_m.center_m[1] + 0.5 * bb_m.height_m,
                                                0.0)));
        next[1] = QuadTree(level_m + 1,
                           BoundingBox(Vector_t(bb_m.center_m[0],
                                                bb_m.center_m[1] - 0.5 * bb_m.height_m,
                                                0.0),
                                       Vector_t(bb_m.center_m[0] + 0.5 * bb_m.width_m,
                                                bb_m.center_m[1],
                                                0.0)));
        next[2] = QuadTree(level_m + 1,
                           BoundingBox(Vector_t(bb_m.center_m[0] - 0.5 * bb_m.width_m,
                                                bb_m.center_m[1],
                                                0.0),
                                       Vector_t(bb_m.center_m[0],
                                                bb_m.center_m[1] + 0.5 * bb_m.height_m,
                                                0.0)));
        next[3] = QuadTree(level_m + 1,
                           BoundingBox(Vector_t(bb_m.center_m[0] - 0.5 * bb_m.width_m,
                                                bb_m.center_m[1] - 0.5 * bb_m.height_m,
                                                0.0),
                                       bb_m.center_m));

        for (unsigned int i = 0; i < 4u; ++ i) {
            next[i].transferIfInside(objects_m);
        }

        bool allEmpty = true;
        for (unsigned int i = 0; i < 4u; ++ i) {
            if (next[i].objects_m.size() != 0) {
                allEmpty = false;
                break;
            }
        }

        if (allEmpty) {
            for (unsigned int i = 0; i < 4u; ++ i) {
                objects_m.merge(next[i].objects_m);
            }

            delete[] next;
            return;
        }

        for (unsigned int i = 0; i < 4u; ++ i) {
            next[i].buildUp();
        }

        nodes_m = next;
    }

    bool QuadTree::isInside(const Vector_t &R) const {
        if (nodes_m != 0) {
            Vector_t X = R - bb_m.center_m;
            unsigned int idx = (X[1] >= 0.0 ? 0: 1);
            idx += (X[0] >= 0.0 ? 0: 2);

            if (nodes_m[idx].isInside(R)) {
                return true;
            }
        }

        for (Base* obj: objects_m) {
            if (obj->isInside(R)) {
                return true;
            }
        }

        return false;
    }

    bool parse(std::string str, Function* &fun) {
        iterator it = str.begin();
        iterator end = str.end();
        if (!parse(it, end, fun)) {
            std::cout << "parsing failed here:" << std::string(it, end) << std::endl;
            return false;
        }

        return true;
    }

    bool parse(iterator &it, const iterator &end, Function* &fun) {
        boost::regex functionCall(FCall);
        boost::smatch what;

        std::string str(it, end);
        if( !boost::regex_match(str , what, functionCall ) ) return false;

        std::string identifier = what[1];
        std::string arguments = what[2];
        unsigned int shift = identifier.size() + 1;

        if (identifier == "rectangle") {
            fun = new Rectangle;
            it += shift;
            if (!Rectangle::parse_detail(it, end, fun)) return false;

            return true;
        } else if (identifier == "ellipse") {
            fun = new Ellipse;
            it += shift;
            if (!Ellipse::parse_detail(it, end, fun)) return false;

            return true;
        } else if (identifier == "polygon") {
            fun = new Polygon;
            it += shift;
            if (!Polygon::parse_detail(it, end, fun)) return false;

            return true;
        } else if (identifier == "mask") {
            fun = new Mask;
            it += shift;

            return Mask::parse_detail(it, end, fun);
        } else if (identifier == "repeat") {
            fun = new Repeat;
            it += shift;
            if (!Repeat::parse_detail(it, end, fun)) return false;

            return true;
        } else if (identifier == "rotate") {
            fun = new Rotation;
            it += shift;
            if (!Rotation::parse_detail(it, end, fun)) return false;

            return true;
        } else if (identifier == "translate") {
            fun = new Translation;
            it += shift;
            if (!Translation::parse_detail(it, end, fun)) return false;

            return true;
        } else if (identifier == "shear") {
            fun = new Shear;
            it += shift;
            if (!Shear::parse_detail(it, end, fun)) return false;

            return true;
        } else if (identifier == "union") {
            fun = new Union;
            it += shift;
            if (!Union::parse_detail(it, end, fun)) return false;

            return true;
        } else if (identifier == "difference") {
            fun = new Difference;
            it += shift;
            if (!Difference::parse_detail(it, end, fun)) return false;

            return true;
        } else if (identifier == "symmetric_difference") {

            fun = new SymmetricDifference;
            it += shift;
            if (!SymmetricDifference::parse_detail(it, end, fun)) return false;

            return true;
        } else if (identifier == "intersection") {
            fun = new Intersection;
            it += shift;
            if (!Intersection::parse_detail(it, end, fun)) return false;

            return true;
        }

        return (it == end);

    }
}