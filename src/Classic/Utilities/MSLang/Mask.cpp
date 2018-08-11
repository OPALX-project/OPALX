#include "Utilities/MSLang/Mask.h"
#include "Utilities/PortableBitmapReader.h"

#include <boost/regex.hpp>
#include <boost/filesystem.hpp>

namespace mslang {
    bool Mask::parse_detail(iterator &it, const iterator &end, Function* &fun) {
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

    void Mask::print(int ident) {
        for (auto pix: pixels_m) pix.print(ident);
    }

    void Mask::apply(std::vector<Base*> &bfuncs) {
        for (auto pix: pixels_m) pix.apply(bfuncs);
    }
}