#include "Algorithms/CoordinateSystemTrafo.h"
#include "Physics/Physics.h"

#include <boost/regex.hpp>

#include <iostream>
#include <string>
#include <fstream>
#include <streambuf>
#include <cstdlib>

std::string UDouble("([0-9]+\\.?[0-9]*([Ee][+-]?[0-9]+)?)");
std::string Double("(-?[0-9]+\\.?[0-9]*([Ee][+-]?[0-9]+)?)");
std::string UInt("([0-9]+)");
std::string FCall("([a-z]*)\\((.*)");

struct base;

struct function {
    virtual ~function() {};

    virtual void print(int indent) = 0;
    virtual void apply(std::vector<base*> &bfuncs) = 0;
};

typedef std::string::iterator iterator;

bool parse(iterator &it, const iterator &end, function* &fun);

struct base: public function {
    CoordinateSystemTrafo trafo;
    base():
        trafo(Vector_t(0.0),
              Quaternion(1.0, 0.0, 0.0, 0.0))
    { }

    virtual base* clone() = 0;
};

struct rectangle: public base {
    double width;
    double height;

    rectangle():
        base(),
        width(0.0),
        height(0.0)
    { }

    virtual ~rectangle() { }

    virtual void print(int indentwidth) {
        std::string indent(indentwidth, ' ');
        std::string indent2(indentwidth + 8, ' ');
        Vector_t origin = trafo.getOrigin();
        double angle = acos(trafo.getRotation()[0]) * 2 * Physics::rad2deg;
        std::cout << indent << "rectangle, \n"
                  << indent2 << "w: " << width << ", \n"
                  << indent2 << "h: " << height << ", \n"
                  << indent2 << "origin: " << origin[0] << ", " << origin[1] << ",\n"
                  << indent2 << "angle: " << angle;
    }

    virtual void apply(std::vector<base*> &bfuncs) {
        bfuncs.push_back(this->clone());
    }

    virtual base* clone() {
        rectangle *rect = new rectangle;
        rect->width = width;
        rect->height = height;
        rect->trafo = trafo;

        return rect;
    }

    static
    bool parse_detail(iterator &it, const iterator &end, function* fun) {
        std::string str(it, end);
        boost::regex argumentList(UDouble + "," + UDouble + "(\\).*)");
        boost::smatch what;

        if (!boost::regex_match(str, what, argumentList)) return false;

        rectangle *rect = static_cast<rectangle*>(fun);
        rect->width = atof(std::string(what[1]).c_str());
        rect->height = atof(std::string(what[3]).c_str());

        std::string fullMatch = what[0];
        std::string rest = what[5];
        it += (fullMatch.size() - rest.size() + 1);

        return true;

    }
};

struct ellipse: public base {
    double width;
    double height;

    ellipse():
        width(0.0),
        height(0.0)
    { }

    virtual ~ellipse() { }

    virtual void print(int indentwidth) {
        std::string indent(indentwidth, ' ');
        std::string indent2(indentwidth + 8, ' ');
        Vector_t origin = trafo.getOrigin();
        double angle = acos(trafo.getRotation()[0]) * 2 * Physics::rad2deg;
        std::cout << indent << "ellipse, \n"
                  << indent2 << "w: " << width << ", \n"
                  << indent2 << "h: " << height << ", \n"
                  << indent2 << "origin: " << origin[0] << ", " << origin[1] << ",\n"
                  << indent2 << "angle: " << angle;
    }

    virtual void apply(std::vector<base*> &bfuncs) {
        bfuncs.push_back(this->clone());
    }

    virtual base* clone() {
        ellipse *elps = new ellipse;
        elps->width = width;
        elps->height = height;
        elps->trafo = trafo;

        return elps;
    }

    static
    bool parse_detail(iterator &it, const iterator &end, function* fun) {
        std::string str(it, end);
        boost::regex argumentList(UDouble + "," + UDouble + "(\\).*)");
        boost::smatch what;

        if (!boost::regex_match(str, what, argumentList)) return false;

        ellipse *elps = static_cast<ellipse*>(fun);
        elps->width = atof(std::string(what[1]).c_str());
        elps->height = atof(std::string(what[3]).c_str());

        std::string fullMatch = what[0];
        std::string rest = what[5];
        it += (fullMatch.size() - rest.size() + 1);

        return true;

    }
};

struct repeat: public function {
    function* func;
    unsigned int N;
    double shiftx;
    double shifty;

    virtual ~repeat() {
        delete func;
    }

    virtual void print(int indentwidth) {
        std::string indent(indentwidth, ' ');
        std::string indent2(indentwidth + 8, ' ');
        std::cout << indent << "repeat, " << std::endl;
        func->print(indentwidth + 8);
        std::cout << ",\n"
                  << indent2 << "N: " << N << ", \n"
                  << indent2 << "dx: " << shiftx << ", \n"
                  << indent2 << "dy: " << shifty;
    }

    virtual void apply(std::vector<base*> &bfuncs) {
        CoordinateSystemTrafo shift(Vector_t(shiftx, shifty, 0.0),
                                    Quaternion(1.0, 0.0, 0.0, 0.0));
        func->apply(bfuncs);
        const unsigned int size = bfuncs.size();

        CoordinateSystemTrafo current_shift = shift;
        for (unsigned int i = 0; i < N; ++ i) {
            for (unsigned int j = 0; j < size; ++ j) {
                base *obj = bfuncs[j]->clone();
                obj->trafo = current_shift * obj->trafo;
                bfuncs.push_back(obj);
            }

            current_shift *= shift;
        }
    }

    static
    bool parse_detail(iterator &it, const iterator &end, function* &fun) {
        repeat *rep = static_cast<repeat*>(fun);
        if (!parse(it, end, rep->func)) return false;

        boost::regex argumentList("," + UInt + "," + Double + "," + Double + "\\)(.*)");
        boost::smatch what;

        std::string str(it, end);
        if (!boost::regex_match(str, what, argumentList)) return false;

        rep->N = atof(std::string(what[1]).c_str());
        rep->shiftx = atof(std::string(what[2]).c_str());
        rep->shifty = atof(std::string(what[4]).c_str());

        std::string fullMatch = what[0];
        std::string rest = what[6];

        it += (fullMatch.size() - rest.size());

        return true;
    }
};

struct translate: public function {
    function* func;
    double shiftx;
    double shifty;

    virtual ~translate() {
        delete func;
    }

    virtual void print(int indentwidth) {
        std::string indent(indentwidth, ' ');
        std::string indent2(indentwidth + 8, ' ');
        std::cout << indent << "translate, " << std::endl;
        func->print(indentwidth + 8);
        std::cout << ",\n"
                  << indent2 << "dx: " << shiftx << ", \n"
                  << indent2 << "dy: " << shifty;
    }

    virtual void apply(std::vector<base*> &bfuncs) {
        CoordinateSystemTrafo shift(Vector_t(shiftx, shifty, 0.0),
                                    Quaternion(1.0, 0.0, 0.0, 0.0));
        func->apply(bfuncs);
        const unsigned int size = bfuncs.size();

        for (unsigned int j = 0; j < size; ++ j) {
            base *obj = bfuncs[j];
            CoordinateSystemTrafo trafo = obj->trafo;
            obj->trafo = shift * obj->trafo;
        }
    }

    static
    bool parse_detail(iterator &it, const iterator &end, function* &fun) {
        translate *trans = static_cast<translate*>(fun);
        if (!parse(it, end, trans->func)) return false;

        boost::regex argumentList("," + Double + "," + Double + "\\)(.*)");
        boost::smatch what;

        std::string str(it, end);
        if (!boost::regex_match(str, what, argumentList)) return false;

        trans->shiftx = atof(std::string(what[1]).c_str());
        trans->shifty = atof(std::string(what[3]).c_str());

        std::string fullMatch = what[0];
        std::string rest = what[5];

        it += (fullMatch.size() - rest.size());

        return true;
    }
};


struct rotate: public function {
    function* func;
    double angle;

    virtual ~rotate() {
        delete func;
    }

    virtual void print(int indentwidth) {
        std::string indent(indentwidth, ' ');
        std::string indent2(indentwidth + 8, ' ');
        std::cout << indent << "rotate, " << std::endl;
        func->print(indentwidth + 8);
        std::cout << ",\n"
                  << indent2 << "angle: " << angle;
    }

    virtual void apply(std::vector<base*> &bfuncs) {
        CoordinateSystemTrafo rotation(Vector_t(0.0, 0.0, 0.0),
                                       Quaternion(cos(0.5 * angle), 0.0, 0.0, sin(0.5 * angle)));
        func->apply(bfuncs);
        const unsigned int size = bfuncs.size();

        for (unsigned int j = 0; j < size; ++ j) {
            base *obj = bfuncs[j];
            obj->trafo = rotation * obj->trafo;
        }
    }

    static
    bool parse_detail(iterator &it, const iterator &end, function* &fun) {
        rotate *rot = static_cast<rotate*>(fun);
        if (!parse(it, end, rot->func)) return false;
        std::cout << std::string(it, end) << std::endl;
        boost::regex argumentList("," + Double + "\\)(.*)");
        boost::smatch what;

        std::string str(it, end);
        if (!boost::regex_match(str, what, argumentList)) return false;

        rot->angle = atof(std::string(what[1]).c_str());

        std::string fullMatch = what[0];
        std::string rest = what[3];

        it += (fullMatch.size() - rest.size());

        return true;
    }
};

struct unionf: public function {
    std::vector<function*> funcs;

    virtual ~unionf () {
        for (function* func: funcs) {
            delete func;
        }
    }

    virtual void print(int indentwidth) {
        std::string indent(indentwidth, ' ');
        std::string indent2(indentwidth + 8, ' ');
        std::string indent3(indentwidth + 16, ' ');
        std::cout << indent << "union, " << std::endl;
        std::cout << indent2 << "funcs: {\n";
        funcs.front()->print(indentwidth + 16);
        for (unsigned int i = 1; i < funcs.size(); ++ i) {
            std::cout << "\n"
                      << indent3 << "," << std::endl;
            funcs[i]->print(indentwidth + 16);
        }
        std::cout << "\n"
                  << indent2 << "} ";
    }

    virtual void apply(std::vector<base*> &bfuncs) {
        for (unsigned int i = 0; i < funcs.size(); ++ i) {
            std::vector<base*> children;
            function *func = funcs[i];
            func->apply(children);
            bfuncs.insert(bfuncs.end(), children.begin(), children.end());
        }
    }

    static
    bool parse_detail(iterator &it, const iterator &end, function* &fun) {
        unionf *unin = static_cast<unionf*>(fun);
        unin->funcs.push_back(NULL);
        if (!parse(it, end, unin->funcs.back())) return false;

        boost::regex argumentList("(,[a-z]+\\(.*)");
        boost::regex endParenthesis("\\)(.*)");
        boost::smatch what;

        std::string str(it, end);
        while (boost::regex_match(str, what, argumentList)) {
            iterator it2 = it + 1;
            unin->funcs.push_back(NULL);

            if (!parse(it2, end, unin->funcs.back())) return false;

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

bool parse(std::string str, function* &fun) {
    str = boost::regex_replace(str, boost::regex("\\s"), std::string(""), boost::match_default | boost::format_all);
    iterator it = str.begin();
    iterator end = str.end();
    if (!parse(it, end, fun)) {
        std::cout << "parsing failed here:" << std::string(it, end) << std::endl;
        return false;
    }

    return true;
}

bool parse(iterator &it, const iterator &end, function* &fun) {
    boost::regex functionCall(FCall);
    boost::smatch what;

    std::string str(it, end);
    if( !boost::regex_match(str , what, functionCall ) ) return false;

    std::string identifier = what[1];
    std::string arguments = what[2];
    unsigned int shift = identifier.size() + 1;

    if (identifier == "rectangle") {
        fun = new rectangle;
        iterator it2 = it + shift;
        if (!rectangle::parse_detail(it2, end, fun)) return false;

        it = it2;
        return true;
    } else if (identifier == "ellipse") {
        fun = new ellipse;
        iterator it2 = it + shift;
        if (!ellipse::parse_detail(it2, end, fun)) return false;

        it = it2;
        return true;
    } else if (identifier == "repeat") {
        fun = new repeat;
        iterator it2 = it + shift;
        if (!repeat::parse_detail(it2, end, fun)) return false;

        it = it2;

        return true;
    } else if (identifier == "rotate") {
        fun = new rotate;
        // iterator it2 =
            it += shift;
        if (!rotate::parse_detail(it, end, fun)) return false;

        // it = it2;

        return true;
    } else if (identifier == "translate") {
        fun = new translate;
        iterator it2 = it + shift;
        if (!translate::parse_detail(it2, end, fun)) return false;

        it = it2;

        return true;
    } else if (identifier == "union") {
        fun = new unionf;
        iterator it2 = it + shift;
        if (!unionf::parse_detail(it2, end, fun)) return false;

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
    function *fun;

    std::ifstream t(argv[1]);
    std::string str((std::istreambuf_iterator<char>(t)),
                    std::istreambuf_iterator<char>());

    // std::string str("repeat( translate(union(rectangle(0.1, 0.1), ellipse(0.1, 0.1)), -0.01, -0.02), 2, 0.1, 0.2)");

    if (parse(str, fun)) {
        fun->print(0);
        std::cout << "\n" << std::endl;

        std::vector<base*> baseBlocks;
        fun->apply(baseBlocks);

        for (base* bfun: baseBlocks) {
            bfun->print(0);
            std::cout << std::endl;
        }

        for (base* func: baseBlocks) {
            delete func;
        }
    }

    delete fun;

    return 0;
}