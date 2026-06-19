#include "PartBunch/Diagnostics/FieldDiagnostics.hpp"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>

#include "AbstractObjects/OpalData.h"
#include "Utilities/Util.h"

extern Inform* gmsg;

template <>
void FieldDiagnostics<double, 3>::dumpVectorField(
        const std::string& what, const VField_t<double, 3>* field) const {
    /*
      what == ef
     */

    Inform m("FieldDiagnostics::dumpVectorField");

    if (ippl::Comm->size() > 1 || callCounter_m < 2) {
        m << level5 << "Skipping vector field dump for multiple ranks or first call." << endl;
        return;
    }

    /* Save the files in the output directory of the simulation. The file
     * name of vector fields is
     *
     * 'basename'-'name'_field-'******'.dat
     *
     * and of scalar fields
     *
     * 'basename'-'name'_scalar-'******'.dat
     *
     * with
     *   'basename': OPAL input file name (*.in)
     *   'name':     field name (input argument of function)
     *   '******':   callCounter_m padded with zeros to 6 digits
     */

    std::string dirname = "data/";

    std::string type;
    std::string unit;

    if (Util::toUpper(what) == "EF") {
        type = "vector";
        unit = "";
    }

    auto localIdx = field->getOwned();
    auto mesh_mp  = &(field->get_mesh());
    auto spacing  = mesh_mp->getMeshSpacing();
    auto origin   = mesh_mp->getOrigin();
    int nghost    = field->getNghost();  // ghosts are excluded in getLocalNDIndex()

    auto fieldV      = field->getView();
    auto field_hostV = field->getHostMirror();
    Kokkos::deep_copy(field_hostV, fieldV);

    std::filesystem::path file(dirname);
    std::string basename = OpalData::getInstance()->getInputBasename();
    std::ostringstream oss;
    oss << basename << "-" << (what + std::string("_") + type) << "-" << std::setfill('0')
        << std::setw(6) << callCounter_m << ".dat";
    std::string filename = oss.str();
    file /= filename;
    std::ofstream fout(file.string(), std::ios::out);

    fout << std::setprecision(9);

    fout << "# " << Util::toUpper(what) << " " << type << " data on grid" << std::endl
         << "# origin= " << std::fixed << origin << " h= " << std::fixed << spacing << std::endl
         << "#" << std::setw(4) << "i" << std::setw(5) << "j" << std::setw(5) << "k"
         << std::setw(17) << "x [m]" << std::setw(17) << "y [m]" << std::setw(17) << "z [m]";

    fout << std::setw(10) << what << "x [" << unit << "]" << std::setw(10) << what << "y [" << unit
         << "]" << std::setw(10) << what << "z [" << unit << "]";

    fout << std::endl;

    for (int i = localIdx[0].first() + nghost; i <= localIdx[0].last() + nghost; i++) {
        for (int j = localIdx[1].first() + nghost; j <= localIdx[1].last() + nghost; j++) {
            for (int k = localIdx[2].first() + nghost; k <= localIdx[2].last() + nghost; k++) {
                // define the physical points (cell-centered)
                double x = (i - nghost) * spacing[0] + origin[0];
                double y = (j - nghost) * spacing[1] + origin[1];
                double z = (k - nghost) * spacing[2] + origin[2];

                fout << std::setw(5) << i << std::setw(5) << j << std::setw(5) << k << std::setw(17)
                     << x << std::setw(17) << y << std::setw(17) << z << std::scientific << "\t"
                     << field_hostV(i, j, k)[0] << "\t" << field_hostV(i, j, k)[1] << "\t"
                     << field_hostV(i, j, k)[2] << std::endl;
            }
        }
    }
    fout.close();
    m << level5 << "*** FINISHED DUMPING " + Util::toUpper(what) + " FIELD *** to " << file.string()
      << endl;
}

template <>
void FieldDiagnostics<double, 3>::dumpScalarField(
        const std::string& what, const Field_t<3>* rho, const Field_t<3>* phi,
        const std::string& solverType) const {
    /*
      what == phi | rho
     */

    Inform m("FieldDiagnostics::dumpScalarField");
    m << level5 << "Dumping scalar field: " << what << endl;

    if (ippl::Comm->size() > 1 /*|| callCounter_m<2*/) {
        m << level5 << "Skipping scalar field dump for multiple ranks or first call." << endl;
        return;
    }

    /* Save the files in the output directory of the simulation. The file
     * name of vector fields is
     *
     * 'basename'-'name'_field-'******'.dat
     *
     * and of scalar fields
     *
     * 'basename'-'name'_scalar-'******'.dat
     *
     * with
     *   'basename': OPAL input file name (*.in)
     *   'name':     field name (input argument of function)
     *   '******':   callCounter_m padded with zeros to 6 digits
     */

    // Needs to be empty...?
    std::string dirname = "data/";

    std::string type;
    std::string unit;
    bool isVectorField = false;

    if (Util::toUpper(what) == "RHO") {
        type = "scalar";
        unit = "Cb/m^3";
    } else if (Util::toUpper(what) == "PHI") {
        type = "scalar";
        unit = "V";
    }

    const Field_t<3>* field = (solverType == "CG" && Util::toUpper(what) == "PHI")
                                      ? phi
                                      : rho;  // both rho and phi are in the same variable (in
                                             // place computation)

    // auto localIdx = field->getOwned();
    ippl::NDIndex<3> localIdx = field->getLayout().getLocalNDIndex();
    int nghost = field->getNghost();  // ghosts are excluded in getLocalNDIndex(), but we still need
                                      // to shift indices
    auto mesh_mp = &(field->get_mesh());
    auto spacing = mesh_mp->getMeshSpacing();
    auto origin  = mesh_mp->getOrigin();

    Field_t<3>::view_type fieldV             = field->getView();
    Field_t<3>::host_mirror_type field_hostV = field->getHostMirror();
    Kokkos::deep_copy(field_hostV, fieldV);

    std::filesystem::path file(dirname);
    std::string basename = OpalData::getInstance()->getInputBasename();
    std::ostringstream oss;
    oss << basename << "-" << (what + std::string("_") + type) << "-" << std::setfill('0')
        << std::setw(6) << callCounter_m << ".dat";
    std::string filename = oss.str();
    file /= filename;
    std::ofstream fout(file.string(), std::ios::out);

    fout << std::setprecision(9);

    fout << "# " << Util::toUpper(what) << " " << type << " data on grid" << std::endl
         << "# origin= " << std::fixed << origin << " h= " << std::fixed << spacing
         << " nghosts=" << nghost << std::endl
         << "#" << std::setw(4) << "i" << std::setw(5) << "j" << std::setw(5) << "k"
         << std::setw(17) << "x [m]" << std::setw(17) << "y [m]" << std::setw(17) << "z [m]";

    if (isVectorField) {
        fout << std::setw(10) << what << "x [" << unit << "]" << std::setw(10) << what << "y ["
             << unit << "]" << std::setw(10) << what << "z [" << unit << "]";
    } else {
        fout << std::setw(13) << what << " [" << unit << "]";
    }

    fout << std::endl;

    if (Util::toUpper(what) == "RHO") {
        for (int i = localIdx[0].first() + nghost; i <= localIdx[0].last() + nghost; i++) {
            for (int j = localIdx[1].first() + nghost; j <= localIdx[1].last() + nghost; j++) {
                for (int k = localIdx[2].first() + nghost; k <= localIdx[2].last() + nghost; k++) {
                    // define the physical points (cell-centered)
                    double x = (i - nghost) * spacing[0] + origin[0];
                    double y = (j - nghost) * spacing[1] + origin[1];
                    double z = (k - nghost) * spacing[2] + origin[2];

                    fout << std::setw(5) << i << std::setw(5) << j << std::setw(5) << k
                         << std::setw(17) << x << std::setw(17) << y << std::setw(17) << z
                         << std::scientific << "\t" << field_hostV(i, j, k) << std::endl;
                }
            }
        }
    } else {
        for (int i = localIdx[0].first() + nghost; i <= localIdx[0].last() + nghost; i++) {
            for (int j = localIdx[1].first() + nghost; j <= localIdx[1].last() + nghost; j++) {
                for (int k = localIdx[2].first() + nghost; k <= localIdx[2].last() + nghost; k++) {
                    // define the physical points (cell-centered)
                    double x = (i - nghost) * spacing[0] + origin[0];
                    double y = (j - nghost) * spacing[1] + origin[1];
                    double z = (k - nghost) * spacing[2] + origin[2];

                    // "+ 1" matches OPAL indexing in the output
                    fout << std::setw(5) << i << std::setw(5) << j << std::setw(5) << k
                         << std::setw(17) << x << std::setw(17) << y << std::setw(17) << z
                         << std::scientific << "\t" << field_hostV(i, j, k) << std::endl;
                }
            }
        }
    }
    fout.close();
    m << level5 << "*** FINISHED DUMPING " + Util::toUpper(what) + " FIELD *** to " << file.string()
      << endl;
}

template class FieldDiagnostics<double, 3>;
