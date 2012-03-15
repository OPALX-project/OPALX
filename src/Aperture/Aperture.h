#ifndef OPAL_Aperture_HH
#define OPAL_Aperture_HH


#include "Algorithms/DefaultVisitor.h"
#include "Beamlines/FlaggedElmPtr.h"
#include "FixedAlgebra/FMatrix.h"
#include "FixedAlgebra/FVector.h"
#include "AbstractObjects/Table.h"
#include "Beamlines/TBeamline.h"
#include "Tables/Twiss.h"
#include <cmath>
#include <vector>
#include <string>

#include <iostream>
using namespace std;

class Aperture: public DefaultVisitor, public Table {
protected:

    //Properties of the optic element
    struct Data {
        double kq, k0s, e1, e2, l, phi, courb;
        string type;
    } data;

    //data needed for the aperture calculation
    double co, K_beta, delta_p, n1;

    double tolerance_x, tolerance_y, offsetx, offsety;

    //vector of points of the beam screen
    struct coord {
        double x;
        double y;
    };
    vector<coord> ShapeBeamScreen;
    vector<coord> ShapeHalo;

public:

    class A_row: public FlaggedElmPtr {
        friend class Aperture;
    public:


        A_row() {};
        A_row(const A_row &a);

        A_row(ElementBase *elem, int occur): FlaggedElmPtr(elem) {
            setCounter(occur);
        }
        A_row(const FlaggedElmPtr &rhs, int order):
            FlaggedElmPtr(rhs), Interpol(order)
        {}
        inline double getBeta_x(int ind) ;
        inline double getBeta_y(int ind) ;
        inline double getDisp_x(int ind) ;
        inline double getDisp_y(int ind) ;
        inline double getDisp_x_prim(int ind) ;
        inline double getDisp_y_prim(int ind) ;
        inline double getApert(int ind);
        inline string getType_elm() ;
        inline double getOrb();

    private:
        string Type_elm;
        double Orb;
        struct pt_interpol {
            double Beta_x, Beta_y;
            double Disp_x, Disp_x_prim;
            double Disp_y, Disp_y_prim;
            double delta_x, delta_y;
            double apert;
        };
        //Vector of pt_interpol of size "order" containing the different values of Beta
        //and Disp... inside the optic element pointed by  p.
        vector<pt_interpol> Interpol;


    };
    typedef TBeamline<A_row> A_Tline;

    Aperture();
    Aperture(const string &name, Aperture *parent);
    virtual ~Aperture() {};

    //return the maximum of the betax(y) function in the elm
    double getBETXMAX(const A_row &, int i1 = 0, int i2 = 0) const;
    double getBETYMAX(const A_row &, int i1 = 0, int i2 = 0) const;
    //return the minimum of the aperture function inside the elm
    double getAPERTMIN(const A_row &row, int i1 = 0, int i2 = 0) const;
    //return the shape of the beam screen or the halo
    vector<Aperture::coord> getShape(vector<double> vec);
    //return the minimum apert of a slice if an element
    void calcul_Apert(A_row &a, int slice, Twiss *tp);

    //return the lenght of the table
    virtual double getLength();

    virtual double getCell(const PlaceRep &place, const string &colName);

    //return the beamline used
    virtual const Beamline *getLine() const;
    //comments in Twiss class
    virtual CellArray getDefault() const;
    virtual std::vector<double> getColumn(const RangeRep &rng, const string
                                          &colName);
    virtual std::vector<double> getRow(const PlaceRep &pos, const
                                       std::vector<string> &cols);
    virtual bool isDependent(const string &name) const;
    virtual bool matches(Table *rhs) const;
    virtual Expressions::PtrToScalar<double> makeColumnExpression(const string &colname) const;
    virtual Object *clone(const string &name);
    virtual void makeTFS(std::ostream &os, const CellArray &cells)const;
    virtual void printTable(std::ostream &, const CellArray &)const;
    virtual void fill();
    void tfsBody(std::ostream &os, const CellArray &cells) const;

    A_Tline::const_iterator begin() const;
    A_Tline::iterator begin();
    A_Tline::const_iterator end() const;
    A_Tline::iterator end();

    const A_row &getCurrent() const;

    //functions in charge of catching the properties of optic elements
    virtual void visitMultipole(const Multipole &);
    virtual void visitRBend(const RBend &);
    virtual void visitSBend(const SBend &);
    virtual void visitCyclotron(const Cyclotron &);
    virtual void applyDefault(const ElementBase &);

    //call of the aperture command
    void execute();
    //execution of the aperture command
    void run();

    void calcul(Twiss::TLine::iterator i, A_row &a, int order, Twiss *tp);
    void PrtTfsApert(A_Tline::iterator n, int order, ofstream
                     &outFile);


protected:

    enum {
        TABLE,
        BEAM,
        NSLICE,
        STATIC,
        DATA,
        DEFAULTAPERTURE,
        FILE,
        SIZE
    };

    A_Tline itsTable;
    A_Tline::iterator n;
    const Beam *beam;
    Twiss *tp;
    Twiss::TLine::iterator i;

    //Transformation matrix for optic elements
    FMatrix<double, 6, 6> Transf_mat;

    //vector alpha,beta,gamma
    FVector<double, 3> Euler_x;
    FVector<double, 3> Euler_y;

    //Variable for dispersion
    FVector<double, 2> Dispx;
    FVector<double, 2> Dispy;

private:
    A_row &findRow(const PlaceRep &row);
    mutable A_Tline::const_iterator current;
    //line name
    string itsLine;
};



#endif //OPAL_Aperture_HH
