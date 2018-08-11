#include "Utilities/MSLang/Shear.h"

#include <boost/regex.hpp>

namespace mslang {
    void Shear::print(int indentwidth) {
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

    void Shear::applyShear(std::vector<Base*> &bfuncs) {
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

    void Shear::apply(std::vector<Base*> &bfuncs) {
        func_m->apply(bfuncs);
        applyShear(bfuncs);
    }

    bool Shear::parse_detail(iterator &it, const iterator &end, Function* &fun) {
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
}