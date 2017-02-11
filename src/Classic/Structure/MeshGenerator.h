#include "AbsBeamline/ElementBase.h"

struct MeshData {
    std::vector<Vector_t> vertices_m;
    std::vector<Vektor<unsigned int, 3> > triangles_m;
    std::vector<std::pair<Vector_t, Vector_t> > decorations_m;
    int type_m;
};


class MeshGenerator {
public:
    MeshGenerator();

    void add(const ElementBase &element);

    void write(const std::string &fname);
private:
    enum MeshType {OTHER = 0
                 , DIPOLE
                 , QUADRUPOLE
                 , SEXTUPOLE
                 , OCTUPOLE
                 , SOLENOID
                 , RFCAVITY
    };

    static MeshData getCylinder(double length,
                                double minor,
                                double major,
                                double formFactor,
                                const unsigned int numSegments = 36);

    static MeshData getBox(double length,
                           double width,
                           double height,
                           double formFactor);


    std::vector<MeshData> elements_m;
};