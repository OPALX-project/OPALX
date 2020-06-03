#include "FieldWriter.h"

#include <iomanip>
#include <fstream>

#include <boost/filesystem.hpp>
#include <boost/format.hpp>

#include "Algorithms/PBunchDefs.h"
#include "Utilities/Util.h"
#include "AbstractObjects/OpalData.h"

template<typename FieldType>
void FieldWriter::dumpField(FieldType& field, std::string name,
                            std::string unit, long long step,
                            FieldType* image)
{
    constexpr bool isVectorField = std::is_same<VField_t, FieldType>::value;
    std::string type = (isVectorField) ? "field" : "scalar";

    INFOMSG("*** START DUMPING " + Util::toUpper(name) + " FIELD ***" << endl);

    boost::filesystem::path file("data");
    boost::format filename("%1%-%2%-%|3$05|.dat");
    std::string basename = OpalData::getInstance()->getInputBasename();
    filename % basename % (name + std::string("_") + type) % step;
    file /= filename.str();

    std::ofstream fout(file.string(), std::ios::out);
    fout.precision(9);

    fout << "# " << name << " " << type << " data on grid" << std::endl
         << "#"
         << std::setw(4)  << "i"
         << std::setw(5)  << "j"
         << std::setw(5)  << "k"
         << std::setw(17) << "x [m]"
         << std::setw(17) << "y [m]"
         << std::setw(17) << "z [m]";
    if (isVectorField) {
         fout << std::setw(10) << name << "x [" << unit << "]"
              << std::setw(10) << name << "y [" << unit << "]"
              << std::setw(10) << name << "z [" << unit << "]";
    } else {
        fout << std::setw(13) << name << " [" << unit << "]";
    }

    if (image) {
        fout << std::setw(13) << name << " image [" << unit << "]";
    }

    fout << std::endl;

    Vector_t origin = field.get_mesh().get_origin();
    Vector_t spacing(field.get_mesh().get_meshSpacing(0),
                     field.get_mesh().get_meshSpacing(1),
                     field.get_mesh().get_meshSpacing(2));

    NDIndex<3> localIdx = field.getLayout().getLocalNDIndex();
    for (int x = localIdx[0].first(); x <= localIdx[0].last(); x++) {
        for (int y = localIdx[1].first(); y <= localIdx[1].last(); y++) {
            for (int z = localIdx[2].first(); z <= localIdx[2].last(); z++) {
                fout << std::setw(5) << x + 1
                     << std::setw(5) << y + 1
                     << std::setw(5) << z + 1
                     << std::setw(17) << origin(0) + x * spacing(0)
                     << std::setw(17) << origin(1) + y * spacing(1)
                     << std::setw(17) << origin(2) + z * spacing(2);
                if (isVectorField) {
                    Vector_t vfield = field[x][y][z].get();
                    fout << std::setw(17) << vfield[0]
                         << std::setw(17) << vfield[1]
                         << std::setw(17) << vfield[2];
                } else {
                    fout << std::setw(17) << field[x][y][z].get();
                }

                if (image) {
                    fout << std::setw(17) << image[x][y][z].get();
                }
                fout << std::endl;
            }
        }
    }
    fout.close();
    INFOMSG("*** FINISHED DUMPING " + Util::toUpper(name) + " FIELD ***" << endl);
}