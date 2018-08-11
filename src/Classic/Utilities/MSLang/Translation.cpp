#include "Utilities/MSLang/Translation.h"
#include "Physics/Physics.h"

#include <boost/regex.hpp>
#include <string>

namespace mslang {
    void Translation::print(int indentwidth) {
        std::string indent(indentwidth, ' ');
        std::string indent2(indentwidth + 8, ' ');
        std::cout << indent << "translate, " << std::endl;
        func_m->print(indentwidth + 8);
        std::cout << ",\n"
                  << indent2 << "dx: " << shiftx_m << ", \n"
                  << indent2 << "dy: " << shifty_m;
    }

    void Translation::applyTranslation(std::vector<Base*> &bfuncs) {
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

    void Translation::apply(std::vector<Base*> &bfuncs) {
        func_m->apply(bfuncs);
        applyTranslation(bfuncs);
    }

    bool Translation::parse_detail(iterator &it, const iterator &end, Function* &fun) {
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
}