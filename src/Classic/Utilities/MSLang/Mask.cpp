#include "Utilities/MSLang/Mask.h"
#include "Utilities/MSLang/ArgumentExtractor.h"
#include "Utilities/MSLang/matheval.h"
#include "Utilities/PortableBitmapReader.h"

#include <boost/regex.hpp>
#include <boost/filesystem.hpp>

namespace mslang {
    bool Mask::parse_detail(iterator &it, const iterator &end, Function* &fun) {
        Mask *pixmap = static_cast<Mask*>(fun);

        ArgumentExtractor arguments(std::string(it, end));
        std::string filename = arguments.get(0);
        if (filename[0] == '\'' && filename.back() == '\'') {
            filename = filename.substr(1, filename.length() - 2);
        }

        if (!boost::filesystem::exists(filename)) {
            ERRORMSG("file '" << filename << "' doesn't exists" << endl);
            return false;
        }

        PortableBitmapReader reader(filename);
        unsigned int width = reader.getWidth();
        unsigned int height = reader.getHeight();

        double pixel_width;
        double pixel_height;
        try {
            pixel_width = parseMathExpression(arguments.get(1)) / width;
            pixel_height = parseMathExpression(arguments.get(2)) / height;
        } catch (std::runtime_error &e) {
            std::cout << e.what() << std::endl;
            return false;
        }

        if (pixel_width < 0.0) {
            std::cout << "Mask: a negative width provided '"
                      << arguments.get(0) << " = " << pixel_width * width << "'"
                      << std::endl;
            return false;
        }

        if (pixel_height < 0.0) {
            std::cout << "Mask: a negative height provided '"
                      << arguments.get(1) << " = " << pixel_height * height << "'"
                      << std::endl;
            return false;
        }

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

        it += (arguments.getLengthConsumed() + 1);

        return true;
    }

    void Mask::print(int ident) {
        for (auto pix: pixels_m) pix.print(ident);
    }

    void Mask::apply(std::vector<Base*> &bfuncs) {
        for (auto pix: pixels_m) pix.apply(bfuncs);
    }
}