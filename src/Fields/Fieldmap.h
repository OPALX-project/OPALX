#ifndef CLASSIC_FIELDMAP_HH
#define CLASSIC_FIELDMAP_HH

#define READ_BUFFER_LENGTH 256

#include <map>
#include <string>
#include <vector>

#include "OPALTypes.h"

#include "Utilities/GSLSpline.h"

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
    T3DMagnetoStatic_Extended,
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

enum DiffDirection { DX = 0, DY, DZ };

/**
 * @class Fieldmap
 * @brief Abstract base class for all field maps. It acts as a factory for
 * creating specific field map types based on the file content and manages a
 * cache of loaded field maps.
 */
class Fieldmap {
public:
    /* ========================================================================== */
    /* ======================== Dictionary Functions ============================ */
    /**
     * @brief Get a field map instance.
     * Use this factory method to obtain a field map. It checks the cache
     * (FieldmapDictionary) first.
     * @param Filename The path to the field map file.
     * @param fast If true, load a "fast" version if available
     * (implementation dependent).
     * @return Pointer to the Fieldmap instance.
     */
    static Fieldmap* getFieldmap(std::string Filename, bool fast = false);

    /**
     * @brief Get a list of all loaded field map names.
     * @return Vector of filenames.
     */
    static std::vector<std::string> getListFieldmapNames();

    /**
     * @brief Delete a specific field map from the cache and memory.
     * @param Filename The filename of the map to delete.
     */
    static void deleteFieldmap(std::string Filename);

    /**
     * @brief Clear the entire field map cache.
     */
    static void clearDictionary();

    /**
     * @brief Read the header of a field map file to determine its type.
     * @param Filename The file to check.
     * @return The MapType determined from the file header.
     */
    static MapType readHeader(std::string Filename);

    /**
     * @brief Trigger the actual reading of the field map data.
     * @param Filename The file to read.
     */
    static void readMap(std::string Filename);

    /**
     * @brief Decrease reference count or delete field map if unused.
     * @param Filename The file to free.
     */
    static void freeMap(std::string Filename);

    static std::string typeset_msg(const std::string& msg, const std::string& title);

    /* ========================================================================== */
    /* =========================== Field Functions=============================== */
    /**
     * @brief Apply the FM to all the particles
     *
     * @param pc Particle container
     */
    virtual void applyField(std::shared_ptr<ParticleContainer_t> pc) = 0;

    /**
     * @brief Get the field strength at a given point.
     *
     * @param R Position [m] relative to the field map origin.
     * @param E Output Electric field [MV/m].
     * @param B Output Magnetic field [T].
     * @return true if R is outside of the field map, false otherwise.
     */
    virtual bool getFieldstrength(
        const Vector_t<double, 3>& R, Vector_t<double, 3>& E, Vector_t<double, 3>& B
    ) const = 0;

    /**
     * @brief Get the field derivative with respect to a direction.
     *
     * @param R Position [m].
     * @param E Output derivative of Electric field.
     * @param B Output derivative of Magnetic field.
     * @param dir Direction of derivative (DX, DY, DZ).
     * @return true if R is outside, false otherwise.
     */
    virtual bool getFieldDerivative(
        const Vector_t<double, 3>& R, Vector_t<double, 3>& E, Vector_t<double, 3>& B,
        const DiffDirection& dir
    ) const = 0;

    /**
     * @brief Get the longitudinal dimensions of the field.
     * @param zBegin Output start of field [m].
     * @param zEnd Output end of field [m].
     */
    virtual void getFieldDimensions(double& zBegin, double& zEnd) const = 0;

    /**
     * @brief Get the full 3D bounding box of the field.
     *
     * @param xIni Output minimum x [m].
     * @param xFinal Output maximum x [m].
     * @param yIni Output minimum y [m].
     * @param yFinal Output maximum y [m].
     * @param zIni Output minimum z [m].
     * @param zFinal Output maximum z [m].
     */
    virtual void getFieldDimensions(
        double& xIni, double& xFinal, double& yIni, double& yFinal, double& zIni, double& zFinal
    ) const = 0;

    /// @brief Swap coordinates (implementation dependent).
    virtual void swap() = 0;

    /// @brief Print info about the field map.
    virtual void getInfo(Inform* msg) = 0;

    /**
     * @brief Get the frequency.
     * @return Frequency [MHz].
     */
    virtual double getFrequency() const = 0;

    /**
     * @brief Set the frequency.
     * @param freq Frequency [MHz].
     */
    virtual void setFrequency(double freq) = 0;

    virtual void setEdgeConstants(
        const double& bendAngle, const double& entranceAngle, const double& exitAngle
    );

    virtual void setFieldLength(const double&);

    virtual void get1DProfile1EngeCoeffs(
        std::vector<double>& engeCoeffsEntry, std::vector<double>& engeCoeffsExit
    );

    virtual void get1DProfile1EntranceParam(
        double& entranceParameter1, double& entranceParameter2, double& entranceParameter3
    );

    virtual void get1DProfile1ExitParam(
        double& exitParameter1, double& exitParameter2, double& exitParameter3
    );

    virtual double getFieldGap();

    virtual void setFieldGap(double gap);

    MapType getType() { return Type; }

    virtual void getOnaxisEz(std::vector<std::pair<double, double>>& onaxis);

    /// @brief Check if a point is inside the field map.
    virtual bool isInside(const Vector_t<double, 3>& /*r*/) const = 0;

    /**
     * @brief Pure virtual method to read the map data.
     * Called by the public static readMap().
     */
    virtual void readMap() = 0;

    /**
     * @brief Pure virtual method to free the map data.
     */
    virtual void freeMap() = 0;

protected:
    Fieldmap() = delete;

    Fieldmap(const std::string& aFilename)
        : Filename_m(aFilename), lines_read_m(0), normalize_m(true) {};

    virtual ~Fieldmap() { ; };
    MapType Type;

    std::string Filename_m;
    int lines_read_m;

    bool normalize_m;
    void getLine(std::ifstream& in, std::string& buffer) { getLine(in, lines_read_m, buffer); }

    static void getLine(std::ifstream& in, int& lines_read, std::string& buffer);

    template <class S>
    bool interpretLine(std::ifstream& in, S& value, const bool& file_length_known = true);

    template <class S, class T>
    bool interpretLine(
        std::ifstream& in, S& value1, T& value2, const bool& file_length_known = true
    );

    template <class S, class T, class U>
    bool interpretLine(
        std::ifstream& in, S& value1, T& value2, U& value3, const bool& file_length_known = true
    );

    template <class S, class T, class U, class V>
    bool interpretLine(
        std::ifstream& in, S& value1, T& value2, U& value3, V& value4,
        const bool& file_length_known = true
    );

    template <class S>
    bool interpretLine(
        std::ifstream& in, S& value1, S& value2, S& value3, S& value4, S& value5, S& value6,
        const bool& file_length_known = true
    );

    bool interpreteEOF(std::ifstream& in);

    void interpretWarning(
        const std::ios_base::iostate& state, const bool& read_all, const std::string& error_msg,
        const std::string& found
    );

    void missingValuesWarning();

    void exceedingValuesWarning();

    void disableFieldmapWarning();

    void noFieldmapWarning();

    void lowResolutionWarning(double squareError, double maxError);

    void checkMap(
        unsigned int accuracy, std::pair<double, double> fieldDimensions, double deltaZ,
        const std::vector<double>& fourierCoefficients, gsl_spline* splineCoefficients,
        gsl_interp_accel* splineAccelerator
    );

    void checkMap(
        unsigned int accuracy, double length, const std::vector<double>& zSampling,
        const std::vector<double>& fourierCoefficients, gsl_spline* splineCoefficients,
        gsl_interp_accel* splineAccelerator
    );

    void write3DField(
        unsigned int nx, unsigned int ny, unsigned int nz, const std::pair<double, double>& xrange,
        const std::pair<double, double>& yrange, const std::pair<double, double>& zrange,
        const std::vector<Vector_t<double, 3>>& ef, const std::vector<Vector_t<double, 3>>& bf
    );

private:
    template <typename T>
    struct TypeParseTraits {
        static const char* name;
    };

    static char buffer_m[READ_BUFFER_LENGTH];
    static std::string alpha_numeric;

    struct FieldmapDescription {
        MapType Type;
        Fieldmap* Map;
        unsigned int RefCounter;
        bool read;
        FieldmapDescription(MapType aType, Fieldmap* aMap)
            : Type(aType), Map(aMap), RefCounter(1), read(false) {}
    };

    static std::map<std::string, FieldmapDescription> FieldmapDictionary;
};

#endif
