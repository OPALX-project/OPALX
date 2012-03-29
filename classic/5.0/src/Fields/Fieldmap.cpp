#include <iostream>
#include <fstream>
#include <ios>
#include "Fields/Fieldmap.hh"
#include "Fields/FM3DDynamic.hh"
#include "Fields/FM3DH5Block.hh"
#include "Fields/FM3DH5Block_nonescale.hh"
#include "Fields/FM3DMagnetoStaticH5Block.hh"
#include "Fields/FM2DDynamic.hh"
#include "Fields/FM2DDynamic_cspline.hh"
#include "Fields/FM2DElectroStatic.hh"
#include "Fields/FM2DElectroStatic_cspline.hh"
#include "Fields/FM2DMagnetoStatic.hh"
#include "Fields/FM2DMagnetoStatic_cspline.hh"
#include "Fields/FM1DDynamic.hh"
#include "Fields/FM1DDynamic_fast.hh"
#include "Fields/Astra1DDynamic.hh"
#include "Fields/Astra1DDynamic_fast.hh"
#include "Fields/FM1DElectroStatic.hh"
#include "Fields/FM1DElectroStatic_fast.hh"
#include "Fields/Astra1DElectroStatic.hh"
#include "Fields/Astra1DElectroStatic_fast.hh"
#include "Fields/FM1DMagnetoStatic.hh"
#include "Fields/FM1DMagnetoStatic_fast.hh"
#include "Fields/Astra1DMagnetoStatic.hh"
#include "Fields/Astra1DMagnetoStatic_fast.hh"
#include "Fields/FM1DProfile1.hh"
#include "Fields/FM1DProfile2.hh"
#include "Fields/FMDummy.hh"
#include "H5hut.h"

#define REGISTER_PARSE_TYPE(X) template <> struct Fieldmap::TypeParseTraits<X> \
    { static const char* name; } ; const char* Fieldmap::TypeParseTraits<X>::name = #X

extern Inform *gmsg;

using namespace std;

Fieldmap *Fieldmap::getFieldmap(std::string Filename, bool fast) {
    map<std::string, FieldmapDescription>::iterator position = FieldmapDictionary.find(Filename);
    if(position != FieldmapDictionary.end()) {
        (*position).second.RefCounter++;
        return (*position).second.Map;
    } else {
        MapType type;
        pair<map<std::string, FieldmapDescription>::iterator, bool> position;
        type = readHeader(Filename);
        switch(type) {
            case T1DDynamic:
                if(fast) {
                    position = FieldmapDictionary.insert(pair<std::string, FieldmapDescription>(Filename, FieldmapDescription(T1DDynamic, new FM1DDynamic_fast(Filename))));
                } else {
                    position = FieldmapDictionary.insert(pair<std::string, FieldmapDescription>(Filename, FieldmapDescription(T1DDynamic, new FM1DDynamic(Filename))));
                }
                return (*position.first).second.Map;
                break;
            case TAstraDynamic:
                if(fast) {
                    position = FieldmapDictionary.insert(pair<std::string, FieldmapDescription>(Filename, FieldmapDescription(TAstraDynamic, new Astra1DDynamic_fast(Filename))));
                } else {
                    position = FieldmapDictionary.insert(pair<std::string, FieldmapDescription>(Filename, FieldmapDescription(TAstraDynamic, new Astra1DDynamic(Filename))));
                }
                return (*position.first).second.Map;
                break;
            case T1DElectroStatic:
                if(fast) {
                    position = FieldmapDictionary.insert(pair<std::string, FieldmapDescription>(Filename, FieldmapDescription(T1DElectroStatic, new FM1DElectroStatic_fast(Filename))));
                } else {
                    position = FieldmapDictionary.insert(pair<std::string, FieldmapDescription>(Filename, FieldmapDescription(T1DElectroStatic, new FM1DElectroStatic(Filename))));
                }
                return (*position.first).second.Map;
                break;
            case TAstraElectroStatic:
                if(fast) {
                    position = FieldmapDictionary.insert(pair<std::string, FieldmapDescription>(Filename, FieldmapDescription(TAstraElectroStatic, new Astra1DElectroStatic_fast(Filename))));
                } else {
                    position = FieldmapDictionary.insert(pair<std::string, FieldmapDescription>(Filename, FieldmapDescription(TAstraElectroStatic, new Astra1DElectroStatic(Filename))));
                }
                return (*position.first).second.Map;
                break;
            case T1DMagnetoStatic:
                if(fast) {
                    position = FieldmapDictionary.insert(pair<std::string, FieldmapDescription>(Filename, FieldmapDescription(T1DMagnetoStatic, new FM1DMagnetoStatic_fast(Filename))));
                } else {
                    position = FieldmapDictionary.insert(pair<std::string, FieldmapDescription>(Filename, FieldmapDescription(T1DMagnetoStatic, new FM1DMagnetoStatic(Filename))));
                }
                return (*position.first).second.Map;
                break;
            case TAstraMagnetoStatic:
                if(fast) {
                    position = FieldmapDictionary.insert(pair<std::string, FieldmapDescription>(Filename, FieldmapDescription(TAstraMagnetoStatic, new Astra1DMagnetoStatic_fast(Filename))));
                } else {
                    position = FieldmapDictionary.insert(pair<std::string, FieldmapDescription>(Filename, FieldmapDescription(TAstraMagnetoStatic, new Astra1DMagnetoStatic(Filename))));
                }
                return (*position.first).second.Map;
                break;
            case T1DProfile1:
                position = FieldmapDictionary.insert(pair<std::string, FieldmapDescription>(Filename, FieldmapDescription(T1DProfile1, new FM1DProfile1(Filename))));
                return (*position.first).second.Map;
                break;
            case T1DProfile2:
                position = FieldmapDictionary.insert(pair<std::string, FieldmapDescription>(Filename, FieldmapDescription(T1DProfile2, new FM1DProfile2(Filename))));
                return (*position.first).second.Map;
                break;
            case T2DDynamic:
                position = FieldmapDictionary.insert(pair<std::string, FieldmapDescription>(Filename, FieldmapDescription(T2DDynamic, new FM2DDynamic(Filename))));
                return (*position.first).second.Map;
                break;
            case T2DDynamic_cspline:
                position = FieldmapDictionary.insert(pair<std::string, FieldmapDescription>(Filename, FieldmapDescription(T2DDynamic_cspline, new FM2DDynamic_cspline(Filename))));
                return (*position.first).second.Map;
                break;
            case T2DElectroStatic:
                position = FieldmapDictionary.insert(pair<std::string, FieldmapDescription>(Filename, FieldmapDescription(T2DElectroStatic, new FM2DElectroStatic(Filename))));
                return (*position.first).second.Map;
                break;
            case T2DElectroStatic_cspline:
                position = FieldmapDictionary.insert(pair<std::string, FieldmapDescription>(Filename, FieldmapDescription(T2DElectroStatic_cspline, new FM2DElectroStatic_cspline(Filename))));
                return (*position.first).second.Map;
                break;
            case T2DMagnetoStatic:
                position = FieldmapDictionary.insert(pair<std::string, FieldmapDescription>(Filename, FieldmapDescription(T2DMagnetoStatic, new FM2DMagnetoStatic(Filename))));
                return (*position.first).second.Map;
                break;
            case T2DMagnetoStatic_cspline:
                position = FieldmapDictionary.insert(pair<std::string, FieldmapDescription>(Filename, FieldmapDescription(T2DMagnetoStatic_cspline, new FM2DMagnetoStatic_cspline(Filename))));
                return (*position.first).second.Map;
                break;
            case T3DDynamic:
                position = FieldmapDictionary.insert(pair<std::string, FieldmapDescription>(Filename, FieldmapDescription(T3DDynamic, new FM3DDynamic(Filename))));
                return (*position.first).second.Map;
                break;
                //        case T3DElectroStatic:
                //            break;
                //        case T3DMagnetoStatic:
                //            break;
            case T3DMagnetoStaticH5Block:
                position = FieldmapDictionary.insert(pair<std::string, FieldmapDescription>(Filename, FieldmapDescription(T3DMagnetoStaticH5Block, new FM3DMagnetoStaticH5Block(Filename))));
                return (*position.first).second.Map;
                break;
            case T3DDynamicH5Block:
                if(fast) {
                    position = FieldmapDictionary.insert(pair<std::string, FieldmapDescription>(Filename, FieldmapDescription(T3DDynamic, new FM3DH5Block_nonescale(Filename))));
                } else {
                    position = FieldmapDictionary.insert(pair<std::string, FieldmapDescription>(Filename, FieldmapDescription(T3DDynamic, new FM3DH5Block(Filename))));
                }
                return (*position.first).second.Map;
                break;
            default:
                position = FieldmapDictionary.insert(pair<std::string, FieldmapDescription>(Filename, FieldmapDescription(UNKNOWN, new FMDummy(Filename))));
                return (*position.first).second.Map;
        }
    }
}

vector<std::string> Fieldmap::getListFieldmapNames() {
    vector<std::string> name_list;
    for(map<std::string, FieldmapDescription>::const_iterator it = FieldmapDictionary.begin(); it != FieldmapDictionary.end(); ++ it) {
        name_list.push_back((*it).first);
    }
    return name_list;
}

void Fieldmap::deleteFieldmap(std::string Filename) {
    map<std::string, FieldmapDescription>::iterator position = FieldmapDictionary.find(Filename);
    /*
      FIXME: find( ) make problem
    0x00000001018f2b40 in std::string::compare ()
    (gdb) where
    #0  0x00000001018f2b40 in std::string::compare ()
    #1  0x0000000100031239 in std::operator< <char, std::char_traits<char>, std::allocator<char> > (__lhs=@0x28, __rhs=@0x7fff5fbfea60) at basic_string.h:2512
    #2  0x0000000100030dd3 in std::less<std::string>::operator() (this=0x100c01d40, __x=@0x28, __y=@0x7fff5fbfea60) at stl_function.h:236
    #3  0x000000010081ad6c in std::_Rb_tree<std::string, std::pair<std::string const, Fieldmap::FieldmapDescription>, std::_Select1st<std::pair<std::string const, Fieldmap::FieldmapDescription> >, std::less<std::string>, std::allocator<std::pair<std::string const, Fieldmap::FieldmapDescription> > >::_M_lower_bound (this=0x100c01d40, __x=0x8, __y=0x1029b80e0, __k=@0x7fff5fbfea60) at stl_tree.h:1080
    #4  0x000000010081a8f5 in std::_Rb_tree<std::string, std::pair<std::string const, Fieldmap::FieldmapDescription>, std::_Select1st<std::pair<std::string const, Fieldmap::FieldmapDescription> >, std::less<std::string>, std::allocator<std::pair<std::string const, Fieldmap::FieldmapDescription> > >::find (this=0x100c01d40, __k=@0x7fff5fbfea60) at stl_tree.h:1526
    #5  0x000000010081a60b in std::map<std::string, Fieldmap::FieldmapDescription, std::less<std::string>, std::allocator<std::pair<std::string const, Fieldmap::FieldmapDescription> > >::find (this=0x100c01d40, __x=@0x7fff5fbfea60) at stl_map.h:737
    #6  0x0000000100817b9a in Fieldmap::deleteFieldmap (Filename={static npos = <optimized out>, _M_dataplus = {<allocator<char>> = {<__gnu_cxx::new_allocator<char>> = {<No data fields>}, <No data fields>}, _M_p = 0x102962608 "CTF3_Ez_ASTRA.opal"}}) at /Users/adelmann/svnwork/opal/classic/5.0/src/Fields/Fieldmap.cpp:167
    #7  0x00000001007204b6 in ~RFCavity (this=0x102962320) at /Users/adelmann/svnwork/opal/classic/5.0/src/AbsBeamline/RFCavity.cpp:121
    #8
        */
    if(position != FieldmapDictionary.end()) {
        if((*position).second.RefCounter > 1) {
            (*position).second.RefCounter--;
        } else {
            delete(*position).second.Map;
            FieldmapDictionary.erase(position);
        }
    }
}

MapType Fieldmap::readHeader(std::string Filename) {
    char magicnumber[5] = "    ";
    std::string buffer;
    int lines_read_m = 0;

    // Check for default map(s).
    if(Filename == "1DPROFILE1-DEFAULT")
        return T1DProfile1;

    ifstream File(Filename.c_str());
    if(!File.good()) {
        cerr << "could not open file " << Filename << endl;
        return UNKNOWN;
    }

    getLine(File, lines_read_m, buffer);
    istringstream interpreter(buffer, istringstream::in);

    interpreter.read(magicnumber, 4);

    if(strcmp(magicnumber, "3DDy") == 0)
        return T3DDynamic;

    if(strcmp(magicnumber, "3DMa") == 0)
        return T3DMagnetoStatic;

    if(strcmp(magicnumber, "3DEl") == 0)
        return T3DElectroStatic;

    if(strcmp(magicnumber, "2DDy") == 0) {
        char tmpString[14] = "             ";
        interpreter.read(tmpString, 13);
        if(strcmp(tmpString, "namic_cspline") == 0) {
            return T2DDynamic_cspline;
        } else {
            return T2DDynamic;
        }
    }

    if(strcmp(magicnumber, "2DMa") == 0) {
        char tmpString[20] = "                   ";
        interpreter.read(tmpString, 19);
        if(strcmp(tmpString, "gnetoStatic_cspline") == 0) {
            return T2DMagnetoStatic_cspline;
        } else {
            return T2DMagnetoStatic;
        }
    }

    if(strcmp(magicnumber, "2DEl") == 0) {
        char tmpString[20] = "                   ";
        interpreter.read(tmpString, 19);
        if(strcmp(tmpString, "ectroStatic_cspline") == 0) {
            return T2DElectroStatic_cspline;
        } else {
            return T2DElectroStatic;
        }
    }

    if(strcmp(magicnumber, "1DDy") == 0)
        return T1DDynamic;

    if(strcmp(magicnumber, "1DMa") == 0)
        return T1DMagnetoStatic;

    if(strcmp(magicnumber, "1DPr") == 0) {
        char tmpString[7] = "      ";
        interpreter.read(tmpString, 6);
        if(strcmp(tmpString, "ofile1") == 0)
            return T1DProfile1;
        if(strcmp(tmpString, "ofile2") == 0)
            return T1DProfile2;
    }

    if(strcmp(magicnumber, "1DEl") == 0)
        return T1DElectroStatic;

    if(strcmp(magicnumber, "\211HDF") == 0) {
        h5_err_t h5err;
        h5_size_t grid_rank;
        h5_size_t grid_dims[3];
        h5_size_t field_dims;
        char name[20];
        h5_size_t len_name = 20;
        h5_int64_t ftype;

        h5_file_t *file = H5OpenFile(Filename.c_str(), H5_O_RDONLY, MPI_COMM_WORLD);
        if(file) {
            h5err = H5SetStep(file, 0);
            if(h5err != H5_SUCCESS)
                ERRORMSG("H5 rc= " << h5err << " in " << __FILE__ << " @ line " << __LINE__ << endl);
            h5_int64_t num_fields = H5BlockGetNumFields(file);
            for(h5_ssize_t i = 0; i < num_fields; ++ i) {
                /*
                  Work around API changes in H5Block
                */
                h5err = H5BlockGetFieldInfo(file, (h5_size_t)i, name, len_name, &grid_rank, grid_dims, &field_dims, &ftype);
                if(h5err != H5_SUCCESS)
                    ERRORMSG("H5 rc= " << h5err << " in " << __FILE__ << " @ line " << __LINE__ << endl);
                if(strcmp(name, "Bfield") == 0) {// using dataset name "Bfield" and "Hfield" to distinguish the type

                    return T3DMagnetoStaticH5Block;
                } else if(strcmp(name, "Hfield") == 0) {
                    return T3DDynamicH5Block;
                }
            }
            h5err = H5CloseFile(file);
            if(h5err != H5_SUCCESS)
                ERRORMSG("H5 rc= " << h5err << " in " << __FILE__ << " @ line " << __LINE__ << endl);
            //return T3DDynamicH5Block;

        }
    }
    if(strcmp(magicnumber, "Astr") == 0) {
        char tmpString[3] = "  ";
        interpreter.read(tmpString, 2);
        if(strcmp(tmpString, "aE") == 0) {
            return TAstraElectroStatic;
        }
        if(strcmp(tmpString, "aM") == 0) {
            return TAstraMagnetoStatic;
        }
        if(strcmp(tmpString, "aD") == 0) {
            return TAstraDynamic;
        }
    }


    return UNKNOWN;
}

void Fieldmap::readMap(std::string Filename) {
    map<std::string, FieldmapDescription>::iterator position = FieldmapDictionary.find(Filename);
    if(position != FieldmapDictionary.end())
        if(!(*position).second.read)
            (*position).second.Map->readMap();
}

void Fieldmap::freeMap(std::string Filename) {
    map<std::string, FieldmapDescription>::iterator position = FieldmapDictionary.find(Filename);
    if(position != FieldmapDictionary.end()) {
        (*position).second.FreeCounter++;
        if((*position).second.FreeCounter == (*position).second.RefCounter) {// Fixme: the initial value of FreeCounter may not be 0 automatically.
            (*position).second.Map->freeMap();
            FieldmapDictionary.erase(position);
        }
    }
}

void Fieldmap::setExitFaceSlope(const double &)
{ }

void Fieldmap::setEdgeConstants(const double &bendAngle, const double &entranceAngle, const double &exitAngle)
{};

void Fieldmap::setFieldGap(const double &)
{};

void Fieldmap::setFieldLength(const double &)
{};

bool Fieldmap::adjustFringeFields()
{ return false; };

void Fieldmap::getLine(ifstream &in, int &lines_read, std::string &buffer) {
    size_t firstof = 0;
    size_t lastof;
    size_t comment;

    do {
        ++ lines_read;
        in.getline(buffer_m, READ_BUFFER_LENGTH);

        buffer = std::string(buffer_m);

        comment = buffer.find("#");
        buffer = buffer.substr(0, comment);

        lastof = buffer.find_last_of(alpha_numeric);
        firstof = buffer.find_first_of(alpha_numeric);
    } while(!in.eof() && lastof == std::string::npos);

    if(firstof != std::string::npos) {
        buffer = buffer.substr(firstof, lastof - firstof + 1);
    }
}

bool Fieldmap::interpreteEOF(ifstream &in) {
    while(!in.eof()) {
        ++lines_read_m;
        in.getline(buffer_m, READ_BUFFER_LENGTH);
        std::string buffer(buffer_m);
        std::string rest;
        size_t comment = buffer.find_first_of("#");
        buffer = buffer.substr(0, comment);
        size_t lasto = buffer.find_first_of(alpha_numeric);
        if(lasto != std::string::npos) {
            exceedingValuesWarning();
            return false;
        }
    }
    return true;
}

void Fieldmap::interpreteWarning(const std::string &error_msg,
                                 const std::string &expecting,
                                 const std::string &found) {
    Inform msg("Fieldmap ");
    std::stringstream errormsg;
    std::stringstream tmpmsg;
    errormsg << "THERE SEEMS TO BE SOMETHING WRONG WITH YOUR FIELD MAP '" << Filename_m << "'.\n"
             << "expecting: '" << expecting << "' on line " << lines_read_m << ",\n"
             << "found instead: '" << found << "'.";
    std::string errormsg_str = typeset_msg(errormsg.str(), "error");
    msg << errormsg_str << "\n"
        << endl;
    if(Ippl::myNode() == 0) {
        ofstream omsg("errormsg.txt", ios_base::app);
        omsg << errormsg_str << endl;
        omsg.close();
    }
}

void Fieldmap::interpreteWarning(const ios_base::iostate &state,
                                 const bool &read_all,
                                 const std::string &expecting,
                                 const std::string &found) {
    std::string error_msg;
    if(!read_all) {
        error_msg = std::string("Didn't find enough values!");
    } else if(state & ios_base::eofbit) {
        error_msg = std::string("Found more values than expected!");
    } else if(state & ios_base::failbit) {
        error_msg = std::string("Found wrong type of values!");
    }
    interpreteWarning(error_msg, expecting, found);
}

void Fieldmap::missingValuesWarning() {
    Inform msg("Fieldmap ");
    stringstream errormsg;
    errormsg << "THERE SEEMS TO BE SOMETHING WRONG WITH YOUR FIELD MAP '" << Filename_m << "'.\n"
             << "There are only " << lines_read_m - 1 << " lines in the file, expecting more.\n"
             << "Please check the section about field maps in the user manual.";
    std::string errormsg_str = typeset_msg(errormsg.str(), "error");

    msg << errormsg_str << "\n"
        << endl;
    if(Ippl::myNode() == 0) {
        ofstream omsg("errormsg.txt", ios_base::app);
        omsg << errormsg_str << endl;
        omsg.close();
    }
}

void Fieldmap::exceedingValuesWarning() {
    Inform msg("Fieldmap ");
    stringstream errormsg;
    errormsg << "THERE SEEMS TO BE SOMETHING WRONG WITH YOUR FIELD MAP '" << Filename_m << "'.\n"
             << "There are too many lines in the file, expecting only " << lines_read_m << " lines.\n"
             << "Please check the section about field maps in the user manual.";
    std::string errormsg_str = typeset_msg(errormsg.str(), "error");
    msg << errormsg_str << "\n"
        << endl;
    if(Ippl::myNode() == 0) {
        ofstream omsg("errormsg.txt", ios_base::app);
        omsg << errormsg_str << endl;
        omsg.close();
    }
}

void Fieldmap::disableFieldmapWarning() {
    Inform msg("Fieldmap ");
    stringstream errormsg;
    errormsg << "DISABLING FIELD MAP '" + Filename_m + "' DUE TO PARSING ERRORS." ;
    std::string errormsg_str = typeset_msg(errormsg.str(), "error");
    msg << errormsg_str << "\n"
        << endl;
    if(Ippl::myNode() == 0) {
        ofstream omsg("errormsg.txt", ios_base::app);
        omsg << errormsg_str << endl;
        omsg.close();
    }
}

void Fieldmap::noFieldmapWarning() {
    Inform msg("Fieldmap ");
    stringstream errormsg;
    errormsg << "DISABLING FIELD MAP '" << Filename_m << "' SINCE FILE COULDN'T BE FOUND!";
    std::string errormsg_str = typeset_msg(errormsg.str(), "error");
    msg << errormsg_str << "\n"
        << endl;
    if(Ippl::myNode() == 0) {
        ofstream omsg("errormsg.txt", ios_base::app);
        omsg << errormsg.str() << endl;
        omsg.close();
    }
}

std::string Fieldmap::typeset_msg(const std::string &msg, const std::string &title) {
    static std::string frame("* ******************************************************************************\n");
    static unsigned int frame_width = frame.length() - 5;
    static std::string closure("                                                                               *\n");

    std::string return_string("\n" + frame);

    int remaining_length = msg.length();
    unsigned int current_position = 0;

    unsigned int ii = 0;
    for(; ii < title.length(); ++ ii) {
        char c = title[ii];
        c = toupper(c);
        return_string.replace(15 + 2 * ii, 1, " ");
        return_string.replace(16 + 2 * ii, 1, &c, 1);
    }
    return_string.replace(15 + 2 * ii, 1, " ");

    while(remaining_length > 0) {
        size_t eol = msg.find("\n", current_position);
        std::string next_to_process;
        if(eol + 1 == msg.length()) {
            break;
        } else if(eol != std::string::npos) {
            next_to_process = msg.substr(current_position, eol - current_position);
        } else {
            next_to_process = msg.substr(current_position);
            eol = msg.length();
        }

        if(eol - current_position < frame_width) {
            return_string += "* " + next_to_process + closure.substr(eol - current_position + 2);
        } else {
            unsigned int last_space = next_to_process.rfind(" ", frame_width);
            if(last_space > 0) {
                if(last_space < frame_width) {
                    return_string += "* " + next_to_process.substr(0, last_space) + closure.substr(last_space + 2);
                } else {
                    return_string += "* " + next_to_process.substr(0, last_space) + " *\n";
                }
                if(next_to_process.length() - last_space + 1 < frame_width) {
                    return_string += "* " + next_to_process.substr(last_space + 1) + closure.substr(next_to_process.length() - last_space + 1);
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

void Fieldmap::getOnaxisEz(vector<pair<double, double> > & onaxis)
{ }


REGISTER_PARSE_TYPE(int);
REGISTER_PARSE_TYPE(double);
REGISTER_PARSE_TYPE(std::string);

std::string Fieldmap::alpha_numeric("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789.-+\211");
map<std::string, Fieldmap::FieldmapDescription> Fieldmap::FieldmapDictionary = map<std::string, Fieldmap::FieldmapDescription>();
char Fieldmap::buffer_m[READ_BUFFER_LENGTH];
