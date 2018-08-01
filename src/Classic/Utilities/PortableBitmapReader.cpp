#include "Utilities/PortableBitmapReader.h"

#include <iostream>
#include <fstream>
#include <sstream>

PortableBitmapReader::PortableBitmapReader(const std::string &input) {

    std::ifstream in(input);
    readHeader(in);
    pixels_m.resize(width_m * height_m);

    if (type_m == ASCII) {
        readImageAscii(in);
    } else {
        readImageBinary(in);
    }
}

std::string PortableBitmapReader::getNextPart(std::istream &in) {
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

void PortableBitmapReader::readHeader(std::istream &in) {
    std::string magicValue = getNextPart(in);

    if (magicValue == "P1") {
        type_m = ASCII;
    } else if (magicValue == "P4") {
        type_m = BINARY;
    } else {
        throw OpalException("PortableBitmapReader::readHeader",
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
}

void PortableBitmapReader::readImageAscii(std::istream &in) {
    unsigned int size = height_m * width_m;
    unsigned int i = 0;
    while (i < size) {
        char c;
        in >> c;

        if (!(c == ' ' ||
              c == '\n' ||
              c == '\t' ||
              c == '\r')) {
            pixels_m[i] = (c == '1');
            ++ i;
        }

    }
}

void PortableBitmapReader::readImageBinary(std::istream &in) {
    static const unsigned int sizeChar = sizeof(char) * 8;

    unsigned int numBytes = width_m / sizeChar;
    if (numBytes * sizeChar != width_m)
        ++ numBytes;

    unsigned int numPixels = 0;
    for (unsigned int i = 0; i < height_m; ++ i) {
        for (unsigned int j = 0; j < numBytes; ++ j) {
            unsigned char c;
            in >> c;

            for (unsigned int j = sizeChar; j > 0 ; -- j) {
                pixels_m[numPixels ++] = ((c >> (j - 1)) & 1);
                if (numPixels % width_m == 0) {
                    break;
                }
            }
        }
    }
}