#ifndef CLASSIC_FIELDMAP_HH
#define CLASSIC_FIELDMAP_HH


#include <string>
#include <map>
#include "Ippl.h"

typedef ParticleSpatialLayout<double,3>::SingleParticlePos_t Vector_t;

using namespace std;

enum MapType {
  UNKNOWN = 0,
  T1DDynamic,
  T1DElectroStatic,
  T1DMagnetoStatic,
  T1DMagnetoStaticEnge,
  T2DDynamic,
  T2DElectroStatic,
  T2DMagnetoStatic,
  T3DDynamic,
  T3DElectroStatic,
  T3DMagnetoStatic
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

class Fieldmap
{

 public:

  static Fieldmap* getFieldmap(string Filename, bool fast=false);
  static void deleteFieldmap(string Filename);
  static MapType readHeader(string Filename);
  static void readMap(string Filename);
  static void freeMap(string Filename);

  virtual bool getFieldstrength(const Vector_t &R, Vector_t &E, Vector_t &B) const = 0;
  virtual void getFieldDimensions(double &zBegin, double &zEnd, double &rBegin, double &rEnd) const = 0;
  virtual void swap() = 0;
  virtual void rescale(double factor) = 0;
  virtual void getInfo(Inform *msg) = 0;
  virtual double getFrequency() const = 0;
  virtual void setFrequency(double freq) = 0;
  MapType getType() { return Type;}

protected:
  Fieldmap(){ ;};
  ~Fieldmap(){ ;};
  MapType Type;

 private:
  virtual void readMap() = 0;
  virtual void freeMap() = 0;

  struct FieldmapDescription{
    MapType Type;
    Fieldmap* Map;
    unsigned int RefCounter;
    unsigned int FreeCounter;
    bool read;
    FieldmapDescription(MapType aType, Fieldmap* aMap)
    {
      Type = aType;
      Map = aMap;
      RefCounter = 1;
      read = false;
    }
//     increaseCounter
  };

  static map<string,FieldmapDescription> FieldmapDictionary;
};


#endif
