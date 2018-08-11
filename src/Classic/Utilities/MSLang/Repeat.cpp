#include "Utilities/MSLang/Repeat.h"

#include <boost/regex.hpp>

namespace mslang {
    void Repeat::print(int indentwidth) {
        std::string indent(indentwidth, ' ');
        std::string indent2(indentwidth + 8, ' ');
        std::cout << indent << "repeat, " << std::endl;
        func_m->print(indentwidth + 8);
        std::cout << ",\n"
                  << indent2 << "N: " << N_m << ", \n"
                  << indent2 << "dx: " << shiftx_m << ", \n"
                  << indent2 << "dy: " << shifty_m;
    }

    void Repeat::apply(std::vector<Base*> &bfuncs) {
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

    bool Repeat::parse_detail(iterator &it, const iterator &end, Function* &fun) {
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
}