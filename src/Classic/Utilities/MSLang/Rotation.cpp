#include "Utilities/MSLang/Rotation.h"

#include <boost/regex.hpp>

namespace mslang {
    void Rotation::print(int indentwidth) {
        std::string indent(indentwidth, ' ');
        std::string indent2(indentwidth + 8, ' ');
        std::cout << indent << "rotate, " << std::endl;
        func_m->print(indentwidth + 8);
        std::cout << ",\n"
                  << indent2 << "angle: " << angle_m;
    }

    void Rotation::applyRotation(std::vector<Base*> &bfuncs) {

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

    void Rotation::apply(std::vector<Base*> &bfuncs) {
        func_m->apply(bfuncs);
        applyRotation(bfuncs);
    }

    bool Rotation::parse_detail(iterator &it, const iterator &end, Function* &fun) {
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
}