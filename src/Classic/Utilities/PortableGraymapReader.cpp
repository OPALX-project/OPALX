#include "Utilities/PortableGraymapReader.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <cstdint>
// #include <set>

PortableGraymapReader::PortableGraymapReader(const std::string &input) {

    std::ifstream in(input);
    readHeader(in);
    pixels_m.resize(width_m * height_m);

    if (type_m == ASCII) {
        readImageAscii(in);
    } else {
        readImageBinary(in);
    }
}

std::string PortableGraymapReader::getNextPart(std::istream &in) {
    char c;
    do {
        c = in.get();
        if (c == '#') {
            do {
                c = in.get();
            } while (c != '\n');
        } else if (!(c == ' ' ||
                     c == '\t' ||
                     c == '\n' ||
                     c == '\r')) {
            in.putback(c);
            break;
        }
    } while (true);

    std::string nextPart;
    in >> nextPart;

    return nextPart;
}

void PortableGraymapReader::readHeader(std::istream &in) {
    std::string magicValue = getNextPart(in);

    if (magicValue == "P2") {
        type_m = ASCII;
    } else if (magicValue == "P5") {
        type_m = BINARY;
    } else {
        throw OpalException("PortableGraymapReader::readHeader",
                            "Unknown magic value: '" + magicValue + "'");
    }

    {
        std::string tmp = getNextPart(in);
        std::istringstream conv;
        conv.str(tmp);
        conv >> width_m;
    }

    {
        std::string tmp = getNextPart(in);
        std::istringstream conv;
        conv.str(tmp);
        conv >> height_m;
    }

    {
        std::string tmp = getNextPart(in);
        std::istringstream conv;
        conv.str(tmp);
        conv >> depth_m;
    }
    // std::cout << width_m << "\t" << height_m << "\t" << depth_m << std::endl;
    char tmp;
    in.read(&tmp, 1);
}

void PortableGraymapReader::readImageAscii(std::istream &in) {
    // std::set<unsigned int> values;
    unsigned int size = height_m * width_m;
    unsigned int i = 0;
    while (i < size) {
        unsigned int c;
        in >> c;

        pixels_m[i] = c / (float)depth_m;
        // values.insert(c);
        ++ i;
    }

    // for (auto val: values) {
    //     std::cout << val << "\t" << val / (float)depth_m << std::endl;
    // }
}

void PortableGraymapReader::readImageBinary(std::istream &in) {
    // static const unsigned int sizeChar = sizeof(char) * 8;

    unsigned int numPixels = 0;
    char c[] = {0, 0};
    unsigned char uc;
    for (unsigned int row = 0; row < height_m; ++ row) {
        for (unsigned int col = 0; col < width_m; ++ col) {
            if ( depth_m > 255) {
                in.read(&c[1], 1);
                in.read(&c[0], 1);
                const uint16_t *val = reinterpret_cast<const uint16_t*>(&c[0]);
                pixels_m[numPixels] = val[0] / (float)depth_m;
            } else {
                in.read(&c[0], 1);
                uc = c[0];
                pixels_m[numPixels] = uc / (float)depth_m;
            }
            ++ numPixels;
        }
    }
}

void PortableGraymapReader::print(std::ostream &out) const {
    unsigned int width = depth_m > 256? 5: 3;
    for (unsigned int i = 0; i < height_m; ++ i) {
        for (unsigned int j = 0; j < width_m; ++ j) {
            unsigned int idx = getIdx(i, j);
            std::cout << "  " << std::setw(width) << pixels_m[idx];
        }
        std::cout << std::endl;
    }
}