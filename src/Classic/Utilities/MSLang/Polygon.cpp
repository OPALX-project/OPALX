#include "Utilities/MSLang/Polygon.h"
#include "Utilities/Mesher.h"
#include "Physics/Physics.h"

#include <boost/regex.hpp>

namespace mslang {
    void Polygon::triangulize(std::vector<Vector_t> &nodes) {
        Mesher mesher(nodes);
        triangles_m = mesher.getTriangles();
    }

    bool Polygon::parse_detail(iterator &it, const iterator &end, Function* &fun) {
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

    void Polygon::print(int ident) {
        // for (auto pix: pixels_m) pix.print(ident);
    }

    void Polygon::apply(std::vector<Base*> &bfuncs) {
        for (Triangle &tri: triangles_m)
            bfuncs.push_back(tri.clone());
    }
}