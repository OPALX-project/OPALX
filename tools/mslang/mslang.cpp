#include "Algorithms/CoordinateSystemTrafo.h"
#include "Physics/Physics.h"

#include <boost/regex.hpp>

#include <iostream>
#include <string>
#include <fstream>
#include <streambuf>
#include <cstdlib>
#include <cmath>
#include <tuple>

std::string UDouble("([0-9]+\\.?[0-9]*([Ee][+-]?[0-9]+)?)");
std::string Double("(-?[0-9]+\\.?[0-9]*([Ee][+-]?[0-9]+)?)");
std::string UInt("([0-9]+)");
std::string FCall("([a-z]*)\\((.*)");

typedef std::string::iterator iterator;

struct Base;

struct Function {
    virtual ~Function() {};

    virtual void print(int indent) = 0;
    virtual void apply(std::vector<Base*> &bfuncs) = 0;
};

bool parse(iterator &it, const iterator &end, Function* &fun);

struct Base: public Function {
    CoordinateSystemTrafo trafo_m;
    // std::tuple<unsigned int,
    //            double,
    //            double> repeat_m;
    Base():
        trafo_m(Vector_t(0.0),
                Quaternion(1.0, 0.0, 0.0, 0.0))
        // , repeat_m(std::make_tuple(0u, 0.0, 0.0))
    { }

    virtual Base* clone() = 0;
    virtual void writeGnuplot(std::ofstream &out) const = 0;
};

struct Rectangle: public Base {
    double width_m;
    double height_m;

    Rectangle():
        Base(),
        width_m(0.0),
        height_m(0.0)
    { }

    virtual ~Rectangle() { }

    virtual void print(int indentwidth) {
        std::string indent(indentwidth, ' ');
        std::string indent2(indentwidth + 8, ' ');
        Vector_t origin = trafo_m.getOrigin();
        double angle = acos(trafo_m.getRotation()[0]) * 2 * Physics::rad2deg;
        std::cout << indent << "rectangle, \n"
                  << indent2 << "w: " << width_m << ", \n"
                  << indent2 << "h: " << height_m << ", \n"
                  << indent2 << "origin: " << origin[0] << ", " << origin[1] << ",\n"
                  << indent2 << "angle: " << angle;
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
    }

    virtual void apply(std::vector<Base*> &bfuncs) {
        bfuncs.push_back(this->clone());
    }

    virtual Base* clone() {
        Rectangle *rect = new Rectangle;
        rect->width_m = width_m;
        rect->height_m = height_m;
        rect->trafo_m = trafo_m;

        return rect;
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
        width_m(0.0),
        height_m(0.0)
    { }

    virtual ~Ellipse() { }

    virtual void print(int indentwidth) {
        std::string indent(indentwidth, ' ');
        std::string indent2(indentwidth + 8, ' ');
        Vector_t origin = trafo_m.getOrigin();
        double angle = acos(trafo_m.getRotation()[0]) * 2 * Physics::rad2deg;
        std::cout << indent << "ellipse, \n"
                  << indent2 << "w: " << width_m << ", \n"
                  << indent2 << "h: " << height_m << ", \n"
                  << indent2 << "origin: " << origin[0] << ", " << origin[1] << ",\n"
                  << indent2 << "angle: " << angle;
    }

    virtual void writeGnuplot(std::ofstream &out) const {
        std::vector<Vector_t> pts({Vector_t(0.5 * width_m, 0.5 * height_m, 0),
                                   Vector_t(-0.5 * width_m, 0.5 * height_m, 0),
                                   Vector_t(-0.5 * width_m, -0.5 * height_m, 0),
                                   Vector_t(0.5 * width_m, -0.5 * height_m, 0)});
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

        for (unsigned int i = 0; i < 5; ++ i) {
            Vector_t pt = pts[i % 4];
            pt = trafo_m.transformFrom(pt);

            out << std::setw(colwidth) << pt[0]
                << std::setw(colwidth) << pt[1]
                << std::endl;
        }
        out << std::endl;
    }

    virtual void apply(std::vector<Base*> &bfuncs) {
        bfuncs.push_back(this->clone());
    }

    virtual Base* clone() {
        Ellipse *elps = new Ellipse;
        elps->width_m = width_m;
        elps->height_m = height_m;
        elps->trafo_m = trafo_m;

        return elps;
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

struct Repeat: public Function {
    Function* func_m;
    unsigned int N_m;
    double shiftx_m;
    double shifty_m;

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
        CoordinateSystemTrafo shift(Vector_t(shiftx_m, shifty_m, 0.0),
                                    Quaternion(1.0, 0.0, 0.0, 0.0));
        func_m->apply(bfuncs);
        const unsigned int size = bfuncs.size();

        CoordinateSystemTrafo current_shift = shift;
        for (unsigned int i = 0; i < N_m; ++ i) {
            for (unsigned int j = 0; j < size; ++ j) {
                Base *obj = bfuncs[j]->clone();
                obj->trafo_m *= current_shift;
                bfuncs.push_back(obj);
            }

            current_shift *= shift;
        }
    }

    static
    bool parse_detail(iterator &it, const iterator &end, Function* &fun) {
        Repeat *rep = static_cast<Repeat*>(fun);
        if (!parse(it, end, rep->func_m)) return false;

        boost::regex argumentList("," + UInt + "," + Double + "," + Double + "\\)(.*)");
        boost::smatch what;

        std::string str(it, end);
        if (!boost::regex_match(str, what, argumentList)) return false;

        rep->N_m = atof(std::string(what[1]).c_str());
        rep->shiftx_m = atof(std::string(what[2]).c_str());
        rep->shifty_m = atof(std::string(what[4]).c_str());

        std::string fullMatch = what[0];
        std::string rest = what[6];

        it += (fullMatch.size() - rest.size());

        return true;
    }
};

struct Translate: public Function {
    Function* func_m;
    double shiftx_m;
    double shifty_m;

    virtual ~Translate() {
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

    virtual void apply(std::vector<Base*> &bfuncs) {
        CoordinateSystemTrafo shift(Vector_t(shiftx_m, shifty_m, 0.0),
                                    Quaternion(1.0, 0.0, 0.0, 0.0));
        func_m->apply(bfuncs);
        const unsigned int size = bfuncs.size();

        for (unsigned int j = 0; j < size; ++ j) {
            Base *obj = bfuncs[j];
            obj->trafo_m *= shift;// * obj->trafo_m;
        }
    }

    static
    bool parse_detail(iterator &it, const iterator &end, Function* &fun) {
        Translate *trans = static_cast<Translate*>(fun);
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


struct Rotate: public Function {
    Function* func_m;
    double angle_m;

    virtual ~Rotate() {
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

    virtual void apply(std::vector<Base*> &bfuncs) {
        CoordinateSystemTrafo rotation(Vector_t(0.0, 0.0, 0.0),
                                       Quaternion(cos(-0.5 * angle_m), 0.0, 0.0, sin(-0.5 * angle_m)));
        func_m->apply(bfuncs);
        const unsigned int size = bfuncs.size();

        for (unsigned int j = 0; j < size; ++ j) {
            Base *obj = bfuncs[j];
            obj->trafo_m *= rotation;
        }
    }

    static
    bool parse_detail(iterator &it, const iterator &end, Function* &fun) {
        Rotate *rot = static_cast<Rotate*>(fun);
        if (!parse(it, end, rot->func_m)) return false;
        std::cout << std::string(it, end) << std::endl;
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


bool parse(std::string str, Function* &fun) {
    str = boost::regex_replace(str, boost::regex("\\s"), std::string(""), boost::match_default | boost::format_all);
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
        iterator it2 = it + shift;
        if (!Rectangle::parse_detail(it2, end, fun)) return false;

        it = it2;
        return true;
    } else if (identifier == "ellipse") {
        fun = new Ellipse;
        iterator it2 = it + shift;
        if (!Ellipse::parse_detail(it2, end, fun)) return false;

        it = it2;
        return true;
    } else if (identifier == "repeat") {
        fun = new Repeat;
        iterator it2 = it + shift;
        if (!Repeat::parse_detail(it2, end, fun)) return false;

        it = it2;

        return true;
    } else if (identifier == "rotate") {
        fun = new Rotate;
        // iterator it2 =
            it += shift;
        if (!Rotate::parse_detail(it, end, fun)) return false;

        // it = it2;

        return true;
    } else if (identifier == "translate") {
        fun = new Translate;
        iterator it2 = it + shift;
        if (!Translate::parse_detail(it2, end, fun)) return false;

        it = it2;

        return true;
    } else if (identifier == "union") {
        fun = new Union;
        iterator it2 = it + shift;
        if (!Union::parse_detail(it2, end, fun)) return false;

        it = it2;

        return true;
    }


    return (it == end);

}

int main(int argc, char *argv[])
{
    if (argc < 2) {
        std::cout << "please provide the name of the file that contains your code" << std::endl;
        return 1;
    }
    Function *fun;

    std::ifstream t(argv[1]);
    std::string str((std::istreambuf_iterator<char>(t)),
                    std::istreambuf_iterator<char>());

    // std::string str("repeat( translate(union(rectangle(0.1, 0.1), ellipse(0.1, 0.1)), -0.01, -0.02), 2, 0.1, 0.2)");

    if (parse(str, fun)) {
        fun->print(0);
        std::cout << "\n" << std::endl;

        std::vector<Base*> baseBlocks;
        fun->apply(baseBlocks);

        std::ofstream out("test.gpl");
        for (Base* bfun: baseBlocks) {
            bfun->print(0);
            std::cout << std::endl;
            bfun->writeGnuplot(out);
        }

        for (Base* func: baseBlocks) {
            delete func;
        }
    }

    delete fun;

    return 0;
}