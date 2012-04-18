#ifndef CLASSIC_FIELDMAP_HH
#define CLASSIC_FIELDMAP_HH

#define READ_BUFFER_LENGTH 256

#include <string>
#include <map>
#include <vector>
#include "Algorithms/Vektor.h"

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

    static Fieldmap *getFieldmap(std::string Filename, bool fast = false);
    static std::vector<std::string> getListFieldmapNames();
    static void deleteFieldmap(std::string Filename);
    static MapType readHeader(std::string Filename);
    static void readMap(std::string Filename);
    static void freeMap(std::string Filename);

    static std::string typeset_msg(const std::string &msg, const std::string &title);

    virtual bool getFieldstrength(const Vector_t &R, Vector_t &E, Vector_t &B) const = 0;
    virtual bool getFieldstrength_fdiff(const Vector_t &R, Vector_t &E, Vector_t &B, const DiffDirection &dir) const = 0;
    virtual void getFieldDimensions(double &zBegin, double &zEnd, double &rBegin, double &rEnd) const = 0;
    virtual void getFieldDimensions(double &xIni, double &xFinal, double &yIni, double &yFinal, double &zIni, double &zFinal) const = 0;
    virtual void swap() = 0;
    virtual void getInfo(Inform *msg) = 0;
    virtual double getFrequency() const = 0;
    virtual void setFrequency(double freq) = 0;
    virtual void setExitFaceSlope(const double &);
    virtual void setEdgeConstants(const double &bendAngle, const double &entranceAngle, const double &exitAngle);
    virtual void setFieldGap(const double &);
    virtual void setFieldLength(const double &);
    virtual bool adjustFringeFields();

    MapType getType() { return Type;}

    virtual void getOnaxisEz(std::vector<std::pair<double, double> > & onaxis);

protected:
    Fieldmap(const std::string &aFilename);
    virtual ~Fieldmap() { ;};
    MapType Type;

    std::string Filename_m;
    int lines_read_m;

    void getLine(std::ifstream &in, std::string &buffer);
    static void getLine(std::ifstream &in, int &lines_read, std::string &buffer);
    template<class S>
    bool interpreteLine(std::ifstream &in, S &value, const bool &file_length_known = true);
    template<class S, class T>
    bool interpreteLine(std::ifstream &in, S &value1, T &value2, const bool &file_length_known = true);
    template<class S, class T, class U>
    bool interpreteLine(std::ifstream &in, S &value1, T &value2, U &value3, const bool &file_length_known = true);
    template<class S, class T, class U, class V>
    bool interpreteLine(std::ifstream &in, S &value1, T &value2, U &value3, V &value4, const bool &file_length_known = true);
    template<class S>
    bool interpreteLine(std::ifstream &in, S &value1, S &value2, S &value3, S &value4, S &value5, S &value6, const bool &file_length_known = true);

    bool interpreteEOF(std::ifstream &in);

    void interpreteWarning(const std::string &error_msg, const std::string &expecting, const std::string &found);
    void interpreteWarning(const std::ios_base::iostate &state,
                           const bool &read_all,
                           const std::string &error_msg,
                           const std::string &found);
    void missingValuesWarning();
    void exceedingValuesWarning();

    void disableFieldmapWarning();
    void noFieldmapWarning();

public:
    virtual void readMap() = 0;
    virtual void freeMap() = 0;
private:
    template<typename T>
    struct TypeParseTraits {
        static const char *name;
    };

    static char buffer_m[READ_BUFFER_LENGTH];
    static std::string alpha_numeric;

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

    static std::map<std::string, FieldmapDescription> FieldmapDictionary;

};

#endif
