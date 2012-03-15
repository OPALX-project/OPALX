#include <iostream>
#include <fstream>
#include <ios>
#include "Fields/Fieldmap.hh"
#include "Fields/FM3DDynamic.hh"
#include "Fields/FM3DH5Block.hh"
#include "Fields/FM2DDynamic.hh"
#include "Fields/FM2DDynamic_cspline.hh"
#include "Fields/FM2DElectroStatic.hh"
#include "Fields/FM2DElectroStatic_cspline.hh"
#include "Fields/FM2DMagnetoStatic.hh"
#include "Fields/FM2DMagnetoStatic_cspline.hh"
#include "Fields/FM1DDynamic.hh"
#include "Fields/FM1DDynamic_fast.hh"
#include "Fields/Astra1DDynamic.hh"
#include "Fields/FM1DElectroStatic.hh"
#include "Fields/FM1DElectroStatic_fast.hh"
#include "Fields/Astra1DElectroStatic.hh"
#include "Fields/FM1DMagnetoStatic.hh"
#include "Fields/FM1DMagnetoStatic_fast.hh"
#include "Fields/Astra1DMagnetoStatic.hh"
#include "Fields/FM1DProfile1.hh"
#include "Fields/FM1DProfile2.hh"
#include "Fields/FMDummy.hh"

#define REGISTER_PARSE_TYPE(X) template <> struct Fieldmap::TypeParseTraits<X> \
    { static const char* name; } ; const char* Fieldmap::TypeParseTraits<X>::name = #X

extern Inform *gmsg;

using namespace std;

Fieldmap* Fieldmap::getFieldmap(string Filename, bool fast)
{
    map<string,FieldmapDescription>::iterator position = FieldmapDictionary.find(Filename);
    if (position != FieldmapDictionary.end()) {
        (*position).second.RefCounter++;
        return (*position).second.Map;
    } else {
        MapType type;
        bool swap;
        pair<map<string,FieldmapDescription>::iterator, bool> position;
        type = readHeader(Filename);
        switch(type) {
        case T1DDynamic:
            if (fast) {
                position = FieldmapDictionary.insert(pair<string,FieldmapDescription>(Filename,FieldmapDescription(T1DDynamic,new FM1DDynamic_fast(Filename))));
            } else {
                position = FieldmapDictionary.insert(pair<string,FieldmapDescription>(Filename,FieldmapDescription(T1DDynamic,new FM1DDynamic(Filename))));
            }
            return (*position.first).second.Map;
            break;
        case TAstraDynamic:
            position = FieldmapDictionary.insert(pair<string,FieldmapDescription>(Filename,FieldmapDescription(TAstraDynamic,new Astra1DDynamic(Filename))));
            return (*position.first).second.Map;
            break;
        case T1DElectroStatic:
            if (fast) {
                position = FieldmapDictionary.insert(pair<string,FieldmapDescription>(Filename,FieldmapDescription(T1DElectroStatic,new FM1DElectroStatic_fast(Filename))));
            } else {
                position = FieldmapDictionary.insert(pair<string,FieldmapDescription>(Filename,FieldmapDescription(T1DElectroStatic,new FM1DElectroStatic(Filename))));
            }
            return (*position.first).second.Map;
            break;
        case TAstraElectroStatic:
            position = FieldmapDictionary.insert(pair<string,FieldmapDescription>(Filename,FieldmapDescription(TAstraElectroStatic,new Astra1DElectroStatic(Filename))));
            return (*position.first).second.Map;
            break;
        case T1DMagnetoStatic:
            if (fast) {
                position = FieldmapDictionary.insert(pair<string,FieldmapDescription>(Filename,FieldmapDescription(T1DMagnetoStatic,new FM1DMagnetoStatic_fast(Filename))));
            } else {
                position = FieldmapDictionary.insert(pair<string,FieldmapDescription>(Filename,FieldmapDescription(T1DMagnetoStatic,new FM1DMagnetoStatic(Filename))));
            }
            return (*position.first).second.Map;
            break;
        case TAstraMagnetoStatic:
            position = FieldmapDictionary.insert(pair<string,FieldmapDescription>(Filename,FieldmapDescription(TAstraMagnetoStatic,new Astra1DMagnetoStatic(Filename))));
            return (*position.first).second.Map;
            break;
        case T1DProfile1:
            position = FieldmapDictionary.insert(pair<string,FieldmapDescription>(Filename,FieldmapDescription(T1DProfile1, new FM1DProfile1(Filename))));
            return (*position.first).second.Map;
            break;
        case T1DProfile2:
            position = FieldmapDictionary.insert(pair<string,FieldmapDescription>(Filename,FieldmapDescription(T1DProfile2,new FM1DProfile2(Filename))));
            return (*position.first).second.Map;
            break;
        case T2DDynamic:
            position = FieldmapDictionary.insert(pair<string,FieldmapDescription>(Filename,FieldmapDescription(T2DDynamic,new FM2DDynamic(Filename))));
            return (*position.first).second.Map;
            break;
        case T2DDynamic_cspline:
            position = FieldmapDictionary.insert(pair<string,FieldmapDescription>(Filename,FieldmapDescription(T2DDynamic_cspline,new FM2DDynamic_cspline(Filename))));
            return (*position.first).second.Map;
            break;
        case T2DElectroStatic:
            position = FieldmapDictionary.insert(pair<string,FieldmapDescription>(Filename,FieldmapDescription(T2DElectroStatic,new FM2DElectroStatic(Filename))));
            return (*position.first).second.Map;
            break;
        case T2DElectroStatic_cspline:
            position = FieldmapDictionary.insert(pair<string,FieldmapDescription>(Filename,FieldmapDescription(T2DElectroStatic_cspline,new FM2DElectroStatic_cspline(Filename))));
            return (*position.first).second.Map;
            break;
        case T2DMagnetoStatic:
            position = FieldmapDictionary.insert(pair<string,FieldmapDescription>(Filename,FieldmapDescription(T2DMagnetoStatic,new FM2DMagnetoStatic(Filename))));
            return (*position.first).second.Map;
            break;
        case T2DMagnetoStatic_cspline:
            position = FieldmapDictionary.insert(pair<string,FieldmapDescription>(Filename,FieldmapDescription(T2DMagnetoStatic_cspline,new FM2DMagnetoStatic_cspline(Filename))));
            return (*position.first).second.Map;
            break;
        case T3DDynamic:
            position = FieldmapDictionary.insert(pair<string,FieldmapDescription>(Filename,FieldmapDescription(T3DDynamic,new FM3DDynamic(Filename))));
            return (*position.first).second.Map;
            break;
//        case T3DElectroStatic:
//            break;
//        case T3DMagnetoStatic:
//            break;
        case T3DDynamicH5Block:
            position = FieldmapDictionary.insert(pair<string,FieldmapDescription>(Filename,FieldmapDescription(T3DDynamic,new FM3DH5Block(Filename))));
            return (*position.first).second.Map;
            break;
        default:
            position = FieldmapDictionary.insert(pair<string, FieldmapDescription>(Filename,FieldmapDescription(UNKNOWN, new FMDummy(Filename))));
            return (*position.first).second.Map;
        }
    }
}

vector<string> Fieldmap::getListFieldmapNames()
{
    vector<string> name_list;
    for (map<string, FieldmapDescription>::const_iterator it = FieldmapDictionary.begin(); it != FieldmapDictionary.end(); ++ it) {
        name_list.push_back((*it).first);
    }
    return name_list;
}

void Fieldmap::deleteFieldmap(string Filename)
{
    map<string,FieldmapDescription>::iterator position = FieldmapDictionary.find(Filename);
    if (position != FieldmapDictionary.end()) {
        if ((*position).second.RefCounter > 1) {
            (*position).second.RefCounter--;
        } else {
            delete (*position).second.Map;
            FieldmapDictionary.erase(position);
        }
    }
}

MapType Fieldmap::readHeader(string Filename)
{
    char magicnumber[5] = "    ";
    string buffer;
    int lines_read_m = 0;

    ifstream File(Filename.c_str());
    if (!File.good()) {
        cerr << "could not open file " << Filename << endl;
        return UNKNOWN;
    }

    getLine(File, lines_read_m, buffer);
    istringstream interpreter(buffer, istringstream::in);

    interpreter.read(magicnumber,4);

    if (strcmp(magicnumber,"3DDy") == 0)
        return T3DDynamic;

    if (strcmp(magicnumber,"3DMa") == 0)
        return T3DMagnetoStatic;

    if (strcmp(magicnumber,"3DEl") == 0)
        return T3DElectroStatic;

    if (strcmp(magicnumber,"2DDy") == 0) {
        char tmpString[14] = "             ";
        interpreter.read(tmpString, 13);
        if (strcmp(tmpString,"namic_cspline") == 0) {
            return T2DDynamic_cspline;
        } else {
            return T2DDynamic;
        }
    }

    if (strcmp(magicnumber,"2DMa") == 0) {
        char tmpString[20] = "                   ";
        interpreter.read(tmpString, 19);
        if (strcmp(tmpString, "gnetoStatic_cspline") == 0) {
            return T2DMagnetoStatic_cspline;
        } else {
            return T2DMagnetoStatic;
        }
    }

    if (strcmp(magicnumber,"2DEl") == 0) {
        char tmpString[20] = "                   ";
        interpreter.read(tmpString, 19);
        if (strcmp(tmpString, "ectroStatic_cspline") == 0) {
            return T2DElectroStatic_cspline;
        } else {
            return T2DElectroStatic;
        }
    }

    if (strcmp(magicnumber,"1DDy") == 0)
        return T1DDynamic;

    if (strcmp(magicnumber,"1DMa") == 0)
        return T1DMagnetoStatic;

    if (strcmp(magicnumber,"1DPr") == 0) {
        char tmpString[7] = "      ";
        interpreter.read(tmpString,6);
        if (strcmp(tmpString,"ofile1") == 0)
            return T1DProfile1;
        if (strcmp(tmpString,"ofile2") == 0)
            return T1DProfile2;
    }

    if (strcmp(magicnumber,"1DEl") == 0)
        return T1DElectroStatic;

    if (strcmp(magicnumber, "HDF\r") == 0) {
        return T3DDynamicH5Block;
    }

    if (strcmp(magicnumber, "Astr") == 0) {
        char tmpString[3] = "  ";
        interpreter.read(tmpString, 2);
        if (strcmp(tmpString, "aE") == 0) {
            return TAstraElectroStatic;
        }
        if (strcmp(tmpString, "aM") == 0) {
            return TAstraMagnetoStatic;
        }
        if (strcmp(tmpString, "aD") == 0) {
            return TAstraDynamic;
        }
    }


    return UNKNOWN;               
}

void Fieldmap::readMap(string Filename)
{
    map<string, FieldmapDescription>::iterator position = FieldmapDictionary.find(Filename);
    if (position != FieldmapDictionary.end())
        if (!(*position).second.read)
            (*position).second.Map->readMap();
}

void Fieldmap::freeMap(string Filename)
{
    map<string, FieldmapDescription>::iterator position = FieldmapDictionary.find(Filename);
    if (position != FieldmapDictionary.end()) {
        (*position).second.FreeCounter++;
        if ((*position).second.FreeCounter == (*position).second.RefCounter) {
            (*position).second.Map->freeMap();
            FieldmapDictionary.erase(position);
        }
    }
}

void Fieldmap::setExitFaceSlope(const double&)
{ }

void Fieldmap::getLine(ifstream & in, int &lines_read, string & buffer)
{
    size_t firstof = 0;
    size_t lastof;
    size_t comment;

    do {
        ++ lines_read;
        in.getline(buffer_m, READ_BUFFER_LENGTH);
        buffer = string(buffer_m);
        comment = buffer.find("#");
        lastof = buffer.find_last_of(alpha_numeric, comment);
        firstof = buffer.find_first_of(alpha_numeric);
    } while (!in.eof() && lastof == string::npos && firstof < READ_BUFFER_LENGTH && firstof < lastof);
    if (firstof < READ_BUFFER_LENGTH && firstof < lastof) {
        buffer = buffer.substr(firstof, lastof + 1);
    }
}

bool Fieldmap::interpreteEOF(ifstream & in) 
{
    while (!in.eof()) {
        ++lines_read_m;
        in.getline(buffer_m, READ_BUFFER_LENGTH);
        string buffer(buffer_m);
        string rest;
        size_t comment = buffer.find_first_of("#");
        buffer = buffer.substr(0, comment);
        size_t lasto = buffer.find_first_of(alpha_numeric);
        if (lasto != string::npos) {
            exceedingValuesWarning();
            return false;
        }
    }
    return true;
}

void Fieldmap::interpreteWarning(const string & error_msg, 
                                 const string & expecting, 
                                 const string & found) 
{
    Inform msg("Fieldmap ");
    stringstream errormsg;
    stringstream tmpmsg;
    double tmp = lines_read_m;
    errormsg << "THERE SEEMS TO BE SOMETHING WRONG WITH YOUR FIELD MAP '" << Filename_m << "'.\n"
             << "expecting: '" << expecting << "' on line " << lines_read_m << ",\n"
             << "found instead: '" << found << "'.";
    string errormsg_str = typeset_msg(errormsg.str(), "error");
    msg << errormsg_str << "\n"
        << endl;
    if (Ippl::myNode() == 0) {
        ofstream omsg("errormsg.txt", ios_base::app);
        omsg << errormsg_str << endl;
        omsg.close();
    }
}

void Fieldmap::interpreteWarning(const ios_base::iostate & state, 
                                 const bool & read_all, 
                                 const string & expecting,
                                 const string & found)
{
    string error_msg;
    if (!read_all) {
        error_msg = string("Didn't find enough values!");
    } else if (state & ios_base::eofbit) {
        error_msg = string("Found more values than expected!");
    } else if (state & ios_base::failbit) {
        error_msg = string("Found wrong type of values!");
    }
    interpreteWarning(error_msg, expecting, found);
}

void Fieldmap::missingValuesWarning() 
{
    Inform msg("Fieldmap ");
    stringstream errormsg;
    errormsg << "THERE SEEMS TO BE SOMETHING WRONG WITH YOUR FIELD MAP '" << Filename_m << "'.\n"
             << "There are only " << lines_read_m - 1 << " lines in the file, expecting more.\n"
             << "Please check the section about field maps in the user manual.";
    string errormsg_str = typeset_msg(errormsg.str(), "error");

    msg << errormsg_str << "\n"
        << endl;
    if (Ippl::myNode() == 0) {
        ofstream omsg("errormsg.txt", ios_base::app);
        omsg << errormsg_str << endl;
        omsg.close();
    }
}

void Fieldmap::exceedingValuesWarning()
{
    Inform msg("Fieldmap ");
    stringstream errormsg;
    errormsg << "THERE SEEMS TO BE SOMETHING WRONG WITH YOUR FIELD MAP '" << Filename_m << "'.\n"
             << "There are too many lines in the file, expecting only " << lines_read_m << " lines.\n"
             << "Please check the section about field maps in the user manual.";
    string errormsg_str = typeset_msg(errormsg.str(), "error");
    msg << errormsg_str << "\n"
        << endl;
    if (Ippl::myNode() == 0) {
        ofstream omsg("errormsg.txt", ios_base::app);
        omsg << errormsg_str << endl;
        omsg.close();
    }
}

void Fieldmap::disableFieldmapWarning()
{
    Inform msg("Fieldmap ");
    stringstream errormsg;
    errormsg << "DISABLING FIELD MAP '" + Filename_m + "' DUE TO PARSING ERRORS." ;
    string errormsg_str = typeset_msg(errormsg.str(), "error");
    msg << errormsg_str << "\n"
        << endl;
    if (Ippl::myNode() == 0) {
        ofstream omsg("errormsg.txt", ios_base::app);
        omsg << errormsg_str << endl;
        omsg.close();
    }
}

void Fieldmap::noFieldmapWarning()
{
    Inform msg("Fieldmap ");
    stringstream errormsg;
    errormsg << "DISABLING FIELD MAP '" << Filename_m << "' SINCE FILE COULDN'T BE FOUND!";
    string errormsg_str = typeset_msg(errormsg.str(), "error");
    msg << errormsg_str << "\n"
        << endl;
    if (Ippl::myNode() == 0) {
        ofstream omsg("errormsg.txt", ios_base::app);
        omsg << errormsg.str() << endl;
        omsg.close();
    }
}

string Fieldmap::typeset_msg(const string & msg, const string & title)
{
    static string frame("* ******************************************************************************\n");
    static int frame_width = frame.length() - 5;
    static string closure("                                                                               *\n");

    string return_string("\n" + frame);

    int remaining_length = msg.length();
    int current_position = 0;

    int ii = 0;
    for (; ii < title.length(); ++ ii) {
        char c = title[ii];
        c = toupper(c);
        return_string.replace(15 + 2*ii, 1, " ");
        return_string.replace(16 + 2*ii, 1, &c, 1);
    }
    return_string.replace(15 + 2*ii, 1, " ");

    while (remaining_length > 0) {
        int eol = msg.find("\n", current_position);
        string next_to_process;
        if (eol == msg.length() - 1) {
            break;
        } else if (eol != string::npos) {
            next_to_process = msg.substr(current_position, eol-current_position);
        } else {
            next_to_process = msg.substr(current_position);
            eol = msg.length();
        }

        if (eol - current_position < frame_width) {
            return_string += "* " + next_to_process + closure.substr(eol - current_position + 2);
        } else {
            int last_space = next_to_process.rfind(" ", frame_width);
            if (last_space > 0) {
                if (last_space < frame_width) {
                    return_string += "* " + next_to_process.substr(0, last_space) + closure.substr(last_space + 2);
                } else {
                    return_string += "* " + next_to_process.substr(0, last_space) + " *\n";
                }
                if (next_to_process.length() - last_space + 1 < frame_width) {
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

REGISTER_PARSE_TYPE(int);
REGISTER_PARSE_TYPE(double);
REGISTER_PARSE_TYPE(string);

string Fieldmap::alpha_numeric("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789.-+");
map<string,Fieldmap::FieldmapDescription> Fieldmap::FieldmapDictionary = map<string,Fieldmap::FieldmapDescription>();
char Fieldmap::buffer_m[READ_BUFFER_LENGTH];
