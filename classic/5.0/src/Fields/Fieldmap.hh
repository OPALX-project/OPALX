#ifndef CLASSIC_FIELDMAP_HH
#define CLASSIC_FIELDMAP_HH

#define READ_BUFFER_LENGTH 256

#include <string>
#include <map>
#include <vector>
#include "Ippl.h"

using namespace std;


typedef ParticleSpatialLayout<double, 3>::SingleParticlePos_t Vector_t;

enum MapType {
    UNKNOWN = 0,
    T1DDynamic,
    TAstraDynamic,
    T1DElectroStatic,
    TAstraElectroStatic,
    T1DMagnetoStatic,
    TAstraMagnetoStatic,
    T1DProfile1,
    T1DProfile2,
    T2DDynamic,
    T2DDynamic_cspline,
    T2DElectroStatic,
    T2DElectroStatic_cspline,
    T2DMagnetoStatic,
    T2DMagnetoStatic_cspline,
    T3DDynamic,
    T3DElectroStatic,
    T3DMagnetoStatic,
    T3DMagnetoStaticH5Block,
    T3DDynamicH5Block
};

enum SwapType {

    XZ = 0,
    ZX,
    XYZ = 10,
    XZMY,
    XMYMZ,
    XMZY,
    YMXZ,
    MXMYZ,
    MYXZ,
    ZYMX,
    MXYMZ,
    MZYX
};

enum DiffDirection {
    DX = 0,
    DY,
    DZ
};

class Fieldmap {

public:

    static Fieldmap *getFieldmap(string Filename, bool fast = false);
    static vector<string> getListFieldmapNames();
    static void deleteFieldmap(string Filename);
    static MapType readHeader(string Filename);
    static void readMap(string Filename);
    static void freeMap(string Filename);

    static std::string typeset_msg(const std::string &msg, const std::string &title);

    virtual bool getFieldstrength(const Vector_t &R, Vector_t &E, Vector_t &B) const = 0;
    virtual bool getFieldstrength_fdiff(const Vector_t &R, Vector_t &E, Vector_t &B, const DiffDirection &dir) const = 0;
    virtual void getFieldDimensions(double &zBegin, double &zEnd, double &rBegin, double &rEnd) const = 0;
    virtual void swap() = 0;
    virtual void getInfo(Inform *msg) = 0;
    virtual double getFrequency() const = 0;
    virtual void setFrequency(double freq) = 0;
    virtual void setExitFaceSlope(const double &);
    MapType getType() { return Type;}

    virtual void getOnaxisEz(vector<pair<double, double> > & onaxis);

protected:
    Fieldmap(const string &aFilename);
    ~Fieldmap() { ;};
    MapType Type;

    string Filename_m;
    int lines_read_m;

    void getLine(ifstream &in, string &buffer);
    static void getLine(ifstream &in, int &lines_read, string &buffer);
    template<class S>
    bool interpreteLine(ifstream &in, S &value, const bool &file_length_known = true);
    template<class S, class T>
    bool interpreteLine(ifstream &in, S &value1, T &value2, const bool &file_length_known = true);
    template<class S, class T, class U>
    bool interpreteLine(ifstream &in, S &value1, T &value2, U &value3, const bool &file_length_known = true);
    template<class S, class T, class U, class V>
    bool interpreteLine(ifstream &in, S &value1, T &value2, U &value3, V &value4, const bool &file_length_known = true);
    template<class S>
    bool interpreteLine(ifstream &in, S &value1, S &value2, S &value3, S &value4, S &value5, S &value6, const bool &file_length_known = true);

    bool interpreteEOF(ifstream &in);

    void interpreteWarning(const string &error_msg, const string &expecting, const string &found);
    void interpreteWarning(const ios_base::iostate &state,
                           const bool &read_all,
                           const string &error_msg,
                           const string &found);
    void missingValuesWarning();
    void exceedingValuesWarning();

    void disableFieldmapWarning();
    void noFieldmapWarning();

private:
    virtual void readMap() = 0;
    virtual void freeMap() = 0;

    template<typename T>
    struct TypeParseTraits {
        static const char *name;
    };

    static char buffer_m[READ_BUFFER_LENGTH];
    static string alpha_numeric;

    struct FieldmapDescription {
        MapType Type;
        Fieldmap *Map;
        unsigned int RefCounter;
        unsigned int FreeCounter;
        bool read;
        FieldmapDescription(MapType aType, Fieldmap *aMap) {
            Type = aType;
            Map = aMap;
            RefCounter = 1;
	    FreeCounter = 0;// fixme: chuan add this because in multipacting simulation FreeCounter can not be properly initialized.
            read = false;
        }
        //     increaseCounter
    };

    static map<string, FieldmapDescription> FieldmapDictionary;

};

#endif
