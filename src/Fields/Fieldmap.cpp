#include "Fields/Fieldmap.h"

#include "Utility/PAssert.h"

/* To be ported from OPAL
#include "Fields/Astra1DElectroStatic.h"
#include "Fields/Astra1DElectroStatic_fast.h"
#include "Fields/Astra1DMagnetoStatic_fast.h"
#include "Fields/Astra1DDynamic_fast.h"
#include "Fields/FM1DDynamic.h"
#include "Fields/FM1DDynamic_fast.h"
#include "Fields/FM1DElectroStatic.h"
#include "Fields/FM1DElectroStatic_fast.h"
#include "Fields/FM1DMagnetoStatic.h"
#include "Fields/FM1DMagnetoStatic_fast.h"
#include "Fields/FM1DProfile1.h"
#include "Fields/FM1DProfile2.h"
#include "Fields/FM2DElectroStatic.h"
#include "Fields/FM3DH5Block.h"
#include "Fields/FM3DH5BlockBase.h"
#include "Fields/FM3DH5Block_nonscale.h"
#include "Fields/FM3DMagnetoStatic.h"
#include "Fields/FM3DMagnetoStaticExtended.h"
#include "Fields/FM3DMagnetoStaticH5Block.h"
#include "Fields/FMDummy.h"
*/

#include "Fields/Astra1DDynamic.h"
#include "Fields/Astra1DMagnetoStatic.h"
#include "Fields/FM2DDynamic.h"
#include "Fields/FM2DMagnetoStatic.h"
#include "Fields/FM3DDynamic.h"
#include "Fields/G4BL2DMagnetoStatic.h"
#include "Fields/G4BL3DGrid.h"

#include "Utilities/GeneralOpalException.h"

#include "H5hut.h"

#include <Kokkos_Core.hpp>

#include <filesystem>

#include <cctype>
#include <cmath>
#include <cstring>
#include <fstream>
#include <ios>
#include <iostream>
#include <new>
#include <sstream>

namespace fs = std::filesystem;

namespace {
    /// Stop with a message that names the map and says what to do. `failure` is the text of
    /// the allocation error, which says how much memory was asked for and where.
    [[noreturn]] void throwOutOfMemory(const std::string& filename, const std::string& failure) {
        throw GeneralOpalException(
                "Fieldmap::readMap",
                "Field map '" + filename + "' does not fit in memory. Loading it failed with:\n  "
                        + failure
                        + "\nEvery MPI rank holds its own copy of the map, and GPU builds keep "
                          "one copy on the CPU and one on the GPU.\nTo make it fit, do one of:\n"
                          "  - use a coarser grid or a smaller box for this map,\n"
                          "  - run fewer MPI ranks per node,\n"
                          "  - run on nodes with more memory.");
    }
}  // namespace

#define REGISTER_PARSE_TYPE(X)            \
    template <>                           \
    struct Fieldmap::TypeParseTraits<X> { \
        static const char* name;          \
    };                                    \
    const char* Fieldmap::TypeParseTraits<X>::name = #X

/**
 * @brief Factory method to obtain a Fieldmap instance.
 *
 * Checks the FieldmapDictionary cache to see if the file is already loaded.
 * If not, it calls readHeader() to determine the map type, instantiates the
 * corresponding subclass, and adds it to the dictionary.
 *
 * @param Filename Absolute path to the fieldmap file.
 * @param fast Not implemented. Previously attempted to load a "fast" version.
 * @param zReverse Load the map mirrored in z with the longitudinal component
 * negated. Only G4beamline cylinder maps support this; the other readers throw
 * when it is set.
 * @return Fieldmap* Pointer to the managed fieldmap instance.
 */
Fieldmap* Fieldmap::getFieldmap(std::string Filename, bool /*fast*/, bool zReverse) {
    std::map<std::string, FieldmapDescription>::iterator position =
            FieldmapDictionary.find(Filename);
    /// Found matching entry?
    if (position != FieldmapDictionary.end()) {
        // The cache is keyed on the filename alone, so the same file cannot be
        // held in two orientations at once.
        if ((*position).second.Map->isZReversed() != zReverse) {
            throw GeneralOpalException(
                    "Fieldmap::getFieldmap()",
                    "Fieldmap \"" + Filename
                            + "\" is already loaded with the opposite ZREVERSE setting. The same "
                              "file cannot be used in both orientations in one run");
        }
        (*position).second.RefCounter++;
        /// Return fieldmap pointer
        return (*position).second.Map;
    } else {
        MapType type;
        std::pair<std::map<std::string, FieldmapDescription>::iterator, bool> position;
        /// Determine fieldmap type
        type = readHeader(Filename);

        if (zReverse && type != TG4BL2DMagnetoStatic) {
            throw GeneralOpalException(
                    "Fieldmap::getFieldmap()",
                    "ZREVERSE is only supported for G4beamline cylinder field maps, not for \""
                            + Filename + "\"");
        }
        /**
         * @brief Based on type:
         * - Create new fieldmap object by parsing the fieldmap file
         * - Create new FieldmapDescription object for this fieldmap
         * - Insert FieldmapDescription into FieldmapDictionary
         * - Return pointer to the fieldmap
         */
        switch (type) {
            case T2DMagnetoStatic:
                position = FieldmapDictionary.insert(
                        std::make_pair(
                                Filename,
                                FieldmapDescription(
                                        T2DMagnetoStatic, new FM2DMagnetoStatic(Filename))));
                return (*position.first).second.Map;
                break;

            case TG4BL2DMagnetoStatic:
                position = FieldmapDictionary.insert(
                        std::make_pair(
                                Filename, FieldmapDescription(
                                                  TG4BL2DMagnetoStatic,
                                                  new G4BL2DMagnetoStatic(Filename, zReverse))));
                return (*position.first).second.Map;
                break;

            case TG4BL3DGrid:
                position = FieldmapDictionary.insert(
                        std::make_pair(
                                Filename,
                                FieldmapDescription(TG4BL3DGrid, new G4BL3DGrid(Filename))));
                return (*position.first).second.Map;
                break;

            case T2DDynamic:
                position = FieldmapDictionary.insert(
                        std::make_pair(
                                Filename,
                                FieldmapDescription(T2DDynamic, new FM2DDynamic(Filename))));
                return (*position.first).second.Map;
                break;

            case T3DDynamic:
                position = FieldmapDictionary.insert(
                        std::make_pair(
                                Filename,
                                FieldmapDescription(T3DDynamic, new FM3DDynamic(Filename))));
                return (*position.first).second.Map;

            case TAstraDynamic:
                // if (fast) {
                //     position = FieldmapDictionary.insert(std::make_pair(
                //         Filename,
                //         FieldmapDescription(TAstraDynamic, new Astra1DDynamic_fast(Filename))));
                // } else {
                position = FieldmapDictionary.insert(
                        std::make_pair(
                                Filename,
                                FieldmapDescription(TAstraDynamic, new Astra1DDynamic(Filename))));
                // }
                return (*position.first).second.Map;
                break;

            case TAstraMagnetoStatic:
                // if (fast) {
                //     position = FieldmapDictionary.insert(
                //             std::make_pair(
                //                 Filename,
                //                 FieldmapDescription(TAstraMagnetoStatic,
                //                 Astra1DMagnetoStatic_fast(Filename))));
                // } else {
                position = FieldmapDictionary.insert(
                        std::make_pair(
                                Filename,
                                FieldmapDescription(
                                        TAstraMagnetoStatic, new Astra1DMagnetoStatic(Filename))));
                // }
                return (*position.first).second.Map;
                break;

            default:
                throw GeneralOpalException(
                        "Fieldmap::getFieldmap()",
                        "Couldn't determine type of fieldmap in file \"" + Filename + "\"");
        }
    }
}

void Fieldmap::clearDictionary() {
    std::map<std::string, FieldmapDescription>::iterator it = FieldmapDictionary.begin();
    for (; it != FieldmapDictionary.end(); ++it) {
        delete it->second.Map;
        it->second.Map = nullptr;
    }
    FieldmapDictionary.clear();
}

/**
 * @brief Determines the fieldmap type
 *
 * @return MapType
 */
MapType Fieldmap::readHeader(std::string Filename) {
    char magicnumber[5] = "    ";
    std::string buffer;
    int lines_read_m = 0;

    // Check for default map(s).
    if (Filename == "1DPROFILE1-DEFAULT") return T1DProfile1;

    if (Filename.empty())
        throw GeneralOpalException("Fieldmap::readHeader()", "No field map file specified");

    if (!fs::exists(Filename))
        throw GeneralOpalException(
                "Fieldmap::readHeader()", "File '" + Filename + "' doesn't exist");

    std::ifstream File(Filename.c_str());
    if (!File.good()) {
        std::cerr << "could not open file " << Filename << std::endl;
        return UNKNOWN;
    }

    getLine(File, lines_read_m, buffer);
    std::istringstream interpreter(buffer, std::istringstream::in);

    interpreter.read(magicnumber, 4);

    if (std::strcmp(magicnumber, "3DDy") == 0) return T3DDynamic;

    if (std::strcmp(magicnumber, "3DMa") == 0) {
        char tmpString[21] = "                    ";
        interpreter.read(tmpString, 20);

        if (std::strcmp(tmpString, "gnetoStatic_Extended") == 0)
            return T3DMagnetoStatic_Extended;
        else
            return T3DMagnetoStatic;
    }

    if (std::strcmp(magicnumber, "3DEl") == 0) return T3DElectroStatic;

    if (std::strcmp(magicnumber, "2DDy") == 0) {
        // char tmpString[14] = "             ";
        // interpreter.read(tmpString, 13);
        return T2DDynamic;
    }

    if (std::strcmp(magicnumber, "2DMa") == 0) {
        // char tmpString[20] = "                   ";
        // interpreter.read(tmpString, 19);
        return T2DMagnetoStatic;
    }

    if (std::strcmp(magicnumber, "2DEl") == 0) {
        // char tmpString[20] = "                   ";
        // interpreter.read(tmpString, 19);
        return T2DElectroStatic;
    }

    if (std::strcmp(magicnumber, "1DDy") == 0) return T1DDynamic;

    if (std::strcmp(magicnumber, "1DMa") == 0) return T1DMagnetoStatic;

    if (std::strcmp(magicnumber, "1DPr") == 0) {
        // char tmpString[7] = "      ";
        // interpreter.read(tmpString, 6);
        // if (strcmp(tmpString, "ofile1") == 0)
        return T1DProfile1;
        // if (strcmp(tmpString, "ofile2") == 0)
        //     return T1DProfile2;
    }

    if (std::strcmp(magicnumber, "1DEl") == 0) return T1DElectroStatic;

    if (std::strcmp(magicnumber, "\211HDF") == 0) {
        h5_err_t h5err [[maybe_unused]] = 0;
        char name[20];
        h5_size_t len_name = sizeof(name);
        h5_prop_t props    = H5CreateFileProp();
        MPI_Comm comm      = ippl::Comm->getCommunicator();
        h5err              = H5SetPropFileMPIOCollective(props, &comm);
        PAssert(h5err != H5_ERR);
        h5_file_t file = H5OpenFile(Filename.c_str(), H5_O_RDONLY, props);
        PAssert(file != (h5_file_t)H5_ERR);
        H5CloseProp(props);

        h5err = H5SetStep(file, 0);
        PAssert(h5err != H5_ERR);

        h5_int64_t num_fields = H5BlockGetNumFields(file);
        PAssert(num_fields != H5_ERR);
        MapType maptype = UNKNOWN;

        for (h5_ssize_t i = 0; i < num_fields; ++i) {
            h5err = H5BlockGetFieldInfo(
                    file, (h5_size_t)i, name, len_name, nullptr, nullptr, nullptr, nullptr);
            PAssert(h5err != H5_ERR);
            // using field name "Bfield" and "Hfield" to distinguish the type
            if (std::strcmp(name, "Bfield") == 0) {
                maptype = T3DMagnetoStaticH5Block;
                break;
            } else if (std::strcmp(name, "Hfield") == 0) {
                maptype = T3DDynamicH5Block;
                break;
            }
        }
        h5err = H5CloseFile(file);
        PAssert(h5err != H5_ERR);
        if (maptype != UNKNOWN) return maptype;
    }

    if (std::strcmp(magicnumber, "Astr") == 0) {
        char tmpString[3] = "  ";
        interpreter.read(tmpString, 2);
        if (std::strcmp(tmpString, "aE") == 0) {
            return TAstraElectroStatic;
        }
        if (std::strcmp(tmpString, "aM") == 0) {
            return TAstraMagnetoStatic;
        }
        if (std::strcmp(tmpString, "aD") == 0) {
            return TAstraDynamic;
        }
    }

    /**
     * G4beamline maps carry no descriptor. They open with an optional `param`
     * line followed by the section keyword that names the grid type, so look
     * for that keyword instead of a magic number.
     */
    if (std::strcmp(magicnumber, "para") == 0 || std::strcmp(magicnumber, "cyli") == 0
        || std::strcmp(magicnumber, "grid") == 0) {
        std::string section = buffer.substr(0, buffer.find_first_of(" \t"));
        while (section == "param" && !File.eof()) {
            getLine(File, lines_read_m, buffer);
            section = buffer.substr(0, buffer.find_first_of(" \t"));
        }

        if (section == "cylinder") {
            return TG4BL2DMagnetoStatic;
        }
        if (section == "grid") {
            return TG4BL3DGrid;
        }
    }

    return UNKNOWN;
}

void Fieldmap::readMap(std::string Filename) {
    std::map<std::string, FieldmapDescription>::iterator position =
            FieldmapDictionary.find(Filename);
    if (position != FieldmapDictionary.end())
        if (!(*position).second.read) {
            // A failed allocation is turned into a message that names the map. freeMap()
            // first releases whatever was allocated before the failure.
            try {
                (*position).second.Map->readMap();
            }
#if KOKKOS_VERSION >= 50200
            // Kokkos 5.2 and later throw this, which is not a std::bad_alloc.
            catch (const Kokkos::Experimental::BadAlloc& e) {
                (*position).second.Map->freeMap();
                throwOutOfMemory(Filename, e.what());
            }
#endif
            catch (const std::bad_alloc& e) {
                (*position).second.Map->freeMap();
                throwOutOfMemory(Filename, e.what());
            }
            (*position).second.read = true;
        }
}

const std::vector<double>& Fieldmap::getFileValues(std::size_t minimum) const {
    if (fileValues_m.size() < minimum) {
        throw GeneralOpalException(
                "Fieldmap::getFileValues", "Field map '" + Filename_m + "': expected at least "
                                                   + std::to_string(minimum)
                                                   + " numbers from the data lines, have "
                                                   + std::to_string(fileValues_m.size()) + ".");
    }
    return fileValues_m;
}

void Fieldmap::freeMap(std::string Filename) {
    std::map<std::string, FieldmapDescription>::iterator position =
            FieldmapDictionary.find(Filename);
    /*
      FIXME: find( ) make problem, crashes
    */
    if (position != FieldmapDictionary.end()) {
        if ((*position).second.RefCounter > 0) {
            (*position).second.RefCounter--;
        }

        if ((*position).second.RefCounter == 0) {
            delete (*position).second.Map;
            (*position).second.Map = nullptr;
            FieldmapDictionary.erase(position);
        }
    }
}

void Fieldmap::getLine(std::ifstream& in, int& lines_read, std::string& buffer) {
    size_t firstof = 0;
    size_t lastof;

    do {
        ++lines_read;
        std::getline(in, buffer);

        size_t comment = buffer.find("#");
        buffer         = buffer.substr(0, comment);

        lastof  = buffer.find_last_of(alpha_numeric);
        firstof = buffer.find_first_of(alpha_numeric);
    } while (in.good() && lastof == std::string::npos);

    if (firstof != std::string::npos) {
        buffer = buffer.substr(firstof, lastof - firstof + 1);
    }
}

bool Fieldmap::interpreteEOF(std::ifstream& in) {
    std::string buffer;
    while (std::getline(in, buffer)) {
        ++lines_read_m;
        size_t comment = buffer.find_first_of("#");
        buffer         = buffer.substr(0, comment);
        size_t lasto   = buffer.find_first_of(alpha_numeric);
        if (lasto != std::string::npos) {
            exceedingValuesWarning();
            return false;
        }
    }
    return true;
}

void Fieldmap::interpretWarning(
        const std::ios_base::iostate& state, const bool& read_all, const std::string& expecting,
        const std::string& found) {
    std::stringstream errormsg;
    errormsg << "THERE SEEMS TO BE SOMETHING WRONG WITH YOUR FIELD MAP '" << Filename_m << "'.\n";
    if (!read_all) {
        errormsg << "Didn't find enough values!" << std::endl;
    } else if (state & std::ios_base::eofbit) {
        errormsg << "Found more values than expected!" << std::endl;
    } else if (state & std::ios_base::failbit) {
        errormsg << "Found wrong type of values!"
                 << "\n"
                 << "expecting: '" << expecting << "' on line " << lines_read_m << ",\n"
                 << "instead found: '" << found << "'." << std::endl;
    }
    throw GeneralOpalException("Fieldmap::interpretWarning()", errormsg.str());
}

void Fieldmap::missingValuesWarning() {
    std::stringstream errormsg;
    errormsg << "THERE SEEMS TO BE SOMETHING WRONG WITH YOUR FIELD MAP '" << Filename_m << "'.\n"
             << "There are only " << lines_read_m - 1 << " lines in the file, expecting more.\n"
             << "Please check the section about field maps in the user manual.";

    throw GeneralOpalException("Fieldmap::missingValuesWarning()", errormsg.str());
}

void Fieldmap::exceedingValuesWarning() {
    std::stringstream errormsg;
    errormsg << "THERE SEEMS TO BE SOMETHING WRONG WITH YOUR FIELD MAP '" << Filename_m << "'.\n"
             << "There are too many lines in the file, expecting only " << lines_read_m
             << " lines.\n"
             << "Please check the section about field maps in the user manual.";

    throw GeneralOpalException("Fieldmap::exceedingValuesWarning()", errormsg.str());
}

void Fieldmap::disableFieldmapWarning() {
    std::stringstream errormsg;
    errormsg << "DISABLING FIELD MAP '" + Filename_m + "' DUE TO PARSING ERRORS.";

    throw GeneralOpalException("Fieldmap::disableFieldmapsWarning()", errormsg.str());
}

void Fieldmap::noFieldmapWarning() {
    std::stringstream errormsg;
    errormsg << "DISABLING FIELD MAP '" << Filename_m << "' SINCE FILE COULDN'T BE FOUND!";

    throw GeneralOpalException("Fieldmap::noFieldmapsWarning()", errormsg.str());
}

std::string Fieldmap::typeset_msg(const std::string& msg, const std::string& title) {
    static std::string frame(
            "* ******************************************************************************\n");
    static unsigned int frame_width = frame.length() - 5;
    static std::string closure(
            "                                                                               *\n");

    std::string return_string("\n" + frame);

    int remaining_length          = msg.length();
    unsigned int current_position = 0;

    unsigned int ii = 0;
    for (; ii < title.length(); ++ii) {
        char c = title[ii];
        c      = std::toupper(c);
        return_string.replace(15 + 2 * ii, 1, " ");
        return_string.replace(16 + 2 * ii, 1, &c, 1);
    }
    return_string.replace(15 + 2 * ii, 1, " ");

    while (remaining_length > 0) {
        size_t eol = msg.find("\n", current_position);
        std::string next_to_process;
        if (eol != std::string::npos) {
            next_to_process = msg.substr(current_position, eol - current_position);
        } else {
            next_to_process = msg.substr(current_position);
            eol             = msg.length();
        }

        if (eol - current_position < frame_width) {
            return_string += "* " + next_to_process + closure.substr(eol - current_position + 2);
        } else {
            unsigned int last_space = next_to_process.rfind(" ", frame_width);
            if (last_space > 0) {
                if (last_space < frame_width) {
                    return_string += "* " + next_to_process.substr(0, last_space)
                                     + closure.substr(last_space + 2);
                } else {
                    return_string += "* " + next_to_process.substr(0, last_space) + " *\n";
                }
                if (next_to_process.length() - last_space + 1 < frame_width) {
                    return_string += "* " + next_to_process.substr(last_space + 1)
                                     + closure.substr(next_to_process.length() - last_space + 1);
                } else {
                    return_string += "* " + next_to_process.substr(last_space + 1) + " *\n";
                }
            } else {
                return_string += "* " + next_to_process + " *\n";
            }
        }

        current_position = eol + 1;
        remaining_length = msg.length() - current_position;
    }
    return_string += frame;

    return return_string;
}

void Fieldmap::getOnaxisEz(std::vector<std::pair<double, double>>& /*onaxis*/) {}

REGISTER_PARSE_TYPE(int);
REGISTER_PARSE_TYPE(unsigned int);
REGISTER_PARSE_TYPE(double);
REGISTER_PARSE_TYPE(std::string);

std::string Fieldmap::alpha_numeric(
        "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789.-+\211");
std::map<std::string, Fieldmap::FieldmapDescription> Fieldmap::FieldmapDictionary =
        std::map<std::string, Fieldmap::FieldmapDescription>();
