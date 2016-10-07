#ifndef PARTBUNCHBASE__H
#define PARTBUNCHBASE__H

// dimension of our positions
const unsigned Dim = 3;

// some typedefs
typedef ParticleSpatialLayout<double,Dim>::SingleParticlePos_t Vector_t;
typedef ParticleSpatialLayout<double,Dim> playout_t;
typedef UniformCartesian<Dim,double> Mesh_t;
typedef Cell                                       Center_t;
typedef CenteredFieldLayout<Dim, Mesh_t, Center_t> FieldLayout_t;
typedef Field<double, Dim, Mesh_t, Center_t>       Field_t;
typedef Field<Vector_t, Dim, Mesh_t, Center_t>     VField_t;
typedef IntCIC IntrplCIC_t;

const double pi = acos(-1.0);




class PartBunchBase {
    
public:
    
    virtual ~PartBunchBase() {}
    
    virtual void myUpdate() = 0;
    virtual void create(int) = 0;
    virtual void gatherStatistics() = 0;
    
    virtual size_t getLocalNum() const = 0;
    virtual Vector_t getRMin() = 0;
    virtual Vector_t getRMax() = 0;
    virtual Vector_t getHr() = 0;

#if 0    
    virtual const Mesh_t& getMesh() const = 0;
    virtual Mesh_t& getMesh() = 0;
    virtual FieldLayout_t& getFieldLayout() = 0;
#endif
    
    virtual double scatter() = 0;
    virtual void initFields() = 0;
    virtual void gatherCIC() = 0;
    
    virtual Vector_t& getR(int i) = 0;
    virtual double& getQM(int i) = 0;
    virtual Vector_t& getP(int i) = 0;
    virtual Vector_t& getE(int i) = 0;
    virtual Vector_t& getB(int i) = 0;
};

#endif