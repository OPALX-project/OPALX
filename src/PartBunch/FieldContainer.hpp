#ifndef OPAL_FIELD_CONTAINER_H
#define OPAL_FIELD_CONTAINER_H

#include <memory>
#include <string>

#include "Manager/BaseManager.h"

template <unsigned Dim>
using Mesh_t = ippl::UniformCartesian<double, Dim>;

template <typename T, unsigned Dim>
using PLayout_t = typename ippl::ParticleSpatialLayout<T, Dim, Mesh_t<Dim>>;

template <unsigned Dim>
using Centering_t = typename Mesh_t<Dim>::DefaultCentering;

template <unsigned Dim>
using FieldLayout_t = ippl::FieldLayout<Dim>;

template <typename T, unsigned Dim, class... ViewArgs>
using Field = ippl::Field<T, Dim, Mesh_t<Dim>, Centering_t<Dim>, ViewArgs...>;

template <unsigned Dim, class... ViewArgs>
using Field_t = Field<double, Dim, ViewArgs...>;

template <typename T, unsigned Dim>
using Vector_t = ippl::Vector<T, Dim>;

template <typename T, unsigned Dim, class... ViewArgs>
using VField_t = Field<Vector_t<T, Dim>, Dim, ViewArgs...>;

// Define the FieldsContainer class
template <typename T, unsigned Dim = 3>
class FieldContainer {
public:
    FieldContainer(
            Vector_t<T, Dim>& hr, Vector_t<T, Dim>& rmin, Vector_t<T, Dim>& rmax,
            std::array<bool, Dim> decomp, ippl::NDIndex<Dim> domain, Vector_t<T, Dim> origin,
            bool isAllPeriodic)
        : hr_m(hr),
          rmin_m(rmin),
          rmax_m(rmax),
          decomp_m(decomp),
          mesh_m(domain, hr, origin),
          fl_m(MPI_COMM_WORLD, domain, decomp, isAllPeriodic) {}

    ~FieldContainer() {}

private:
    Vector_t<double, Dim> hr_m;
    Vector_t<double, Dim> rmin_m;
    Vector_t<double, Dim> rmax_m;
    std::array<bool, Dim> decomp_m;
    VField_t<T, Dim> E_m;
    Field_t<Dim> rho_m;
    Field<T, Dim> phi_m;
    std::shared_ptr<VField_t<T, Dim>> Etmp_m;
    std::shared_ptr<VField_t<T, Dim>> Btmp_m;
    std::shared_ptr<VField_t<T, Dim>> flippedZSlabField_m;
    Mesh_t<Dim> mesh_m;
    FieldLayout_t<Dim> fl_m;

public:
    VField_t<T, Dim>& getE() { return E_m; }
    void setE(VField_t<T, Dim>& E) { E_m = E; }

    Field_t<Dim>& getRho() { return rho_m; }
    void setRho(Field_t<Dim>& rho) { rho_m = rho; }

    Field<T, Dim>& getPhi() { return phi_m; }
    void setPhi(Field<T, Dim>& phi) { phi_m = phi; }

    Vector_t<double, Dim>& getHr() { return hr_m; }
    void setHr(const Vector_t<double, Dim>& hr) { hr_m = hr; }

    Vector_t<double, Dim>& getRMin() { return rmin_m; }
    void setRMin(const Vector_t<double, Dim>& rmin) { rmin_m = rmin; }

    Vector_t<double, Dim>& getRMax() { return rmax_m; }
    void setRMax(const Vector_t<double, Dim>& rmax) { rmax_m = rmax; }

    std::array<bool, Dim> getDecomp() { return decomp_m; }
    void setDecomp(std::array<bool, Dim> decomp) { decomp_m = decomp; }

    Mesh_t<Dim>& getMesh() { return mesh_m; }
    void setMesh(Mesh_t<Dim>& mesh) { mesh_m = mesh; }

    FieldLayout_t<Dim>& getFL() { return fl_m; }
    void setFL(std::shared_ptr<FieldLayout_t<Dim>>& fl) { fl_m = fl; }

    std::shared_ptr<VField_t<T, Dim>> getTempEField() { return Etmp_m; }
    void setTempEField(std::shared_ptr<VField_t<T, Dim>> Etmp) { Etmp_m = Etmp; }

    std::shared_ptr<VField_t<T, Dim>> getTempBField() { return Btmp_m; }
    void setTempBField(std::shared_ptr<VField_t<T, Dim>> Btmp) { Btmp_m = Btmp; }

    std::shared_ptr<VField_t<T, Dim>> getFlippedZSlabField() { return flippedZSlabField_m; }
    void resetFlippedZSlabField() { flippedZSlabField_m.reset(); }

    void initializeTemporaryFields() {
        if (!Etmp_m) {
            Etmp_m = std::make_shared<VField_t<T, Dim>>();
        }
        Etmp_m->initialize(mesh_m, fl_m);

        if (!Btmp_m) {
            Btmp_m = std::make_shared<VField_t<T, Dim>>();
        }
        Btmp_m->initialize(mesh_m, fl_m);
    }

    void updateFieldLayoutsAfterLayoutChange(const std::string& stype_m = "") {
        E_m.updateLayout(fl_m);
        rho_m.updateLayout(fl_m);

        if (stype_m == "CG") {
            phi_m.updateLayout(fl_m);
            phi_m = 0.0;
            phi_m.setFieldBC(phi_m.getFieldBC());
        }

        if (Etmp_m) {
            Etmp_m->updateLayout(fl_m);
        } else {
            Etmp_m = std::make_shared<VField_t<T, Dim>>();
            Etmp_m->initialize(mesh_m, fl_m);
        }

        if (Btmp_m) {
            Btmp_m->updateLayout(fl_m);
        } else {
            Btmp_m = std::make_shared<VField_t<T, Dim>>();
            Btmp_m->initialize(mesh_m, fl_m);
        }

        resetFlippedZSlabField();
    }

    std::shared_ptr<VField_t<T, Dim>> getOrCreateFlippedZSlabField(
            const VField_t<T, Dim>& src) {
        auto& layout        = src.getLayout();
        auto& mesh          = src.get_mesh();
        const int srcNghost = src.getNghost();

        const bool needsInit = !flippedZSlabField_m
                               || &flippedZSlabField_m->getLayout() != &layout
                               || flippedZSlabField_m->getNghost() != srcNghost;
        if (!flippedZSlabField_m) {
            flippedZSlabField_m = std::make_shared<VField_t<T, Dim>>();
        }
        if (needsInit) {
            flippedZSlabField_m->initialize(mesh, layout, srcNghost);
        }

        return flippedZSlabField_m;
    }

    void initializeFields(const std::string& stype_m = "") {
        Inform m("FieldContainer::initializeFields");
        m << level3 << "Mesh spacing = " << mesh_m.getMeshSpacing() << endl;
        m << level3 << "Origin       = " << mesh_m.getOrigin() << endl;
        m << level3 << "FL           = " << fl_m << endl;

        E_m.initialize(mesh_m, fl_m);
        rho_m.initialize(mesh_m, fl_m);
        m << level3 << "E_m, rho_m field initialized." << endl;
        if (stype_m == "CG") {
            phi_m.initialize(mesh_m, fl_m);
            m << level3 << "Phi field initialized for " << stype_m << endl;
        }
        initializeTemporaryFields();
        resetFlippedZSlabField();
    }
};

#endif
