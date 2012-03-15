#include <iostream>
#include <fstream>
#include <ios>
#include "Fields/Fieldmap.hh"
#include "Fields/FM2DElectroStatic.hh"
#include "Fields/FM2DMagnetoStatic.hh"
#include "Fields/FM2DDynamic.hh"
#include "Fields/FM1DDynamic.hh"
#include "Fields/FM1DDynamic_fast.hh"

using namespace std;


Fieldmap* Fieldmap::getFieldmap(string Filename, bool fast)
{
  map<string,FieldmapDescription>::iterator position = FieldmapDictionary.find(Filename);
  if (position != FieldmapDictionary.end())
    {
      (*position).second.RefCounter++;
      return (*position).second.Map;
    }
  else
    {
      MapType type;
      bool swap;
      pair<map<string,FieldmapDescription>::iterator, bool> position;
      type = readHeader(Filename);
      switch(type)
        {
        case T1DDynamic:
          if (fast)
              position = FieldmapDictionary.insert(pair<string,FieldmapDescription>(Filename,FieldmapDescription(T1DDynamic,new FM1DDynamic_fast(Filename))));
          else
              position = FieldmapDictionary.insert(pair<string,FieldmapDescription>(Filename,FieldmapDescription(T1DDynamic,new FM1DDynamic(Filename))));
          return (*position.first).second.Map;
          break;
        case T1DElectroStatic:
//           position = FieldmapDictionary.insert(pair<string,FieldmapDescription>(Filename,FieldmapDescription(T1DElectroStatic,new FM1DElectroStatic(Filename))));
//           return (*position.first).second.Map;
          break;
        case T1DMagnetoStatic:
//           position = FieldmapDictionary.insert(pair<string,FieldmapDescription>(Filename,FieldmapDescription(T1DMagnetoStatic,new FM1DMagnetoStatic(Filename))));
//           return (*position.first).second.Map;
          break;
        case T2DDynamic:
          position = FieldmapDictionary.insert(pair<string,FieldmapDescription>(Filename,FieldmapDescription(T2DDynamic,new FM2DDynamic(Filename))));
          return (*position.first).second.Map;
          break;
        case T2DElectroStatic:
          position = FieldmapDictionary.insert(pair<string,FieldmapDescription>(Filename,FieldmapDescription(T2DElectroStatic,new FM2DElectroStatic(Filename))));
          return (*position.first).second.Map;
          break;
        case T2DMagnetoStatic:
          position = FieldmapDictionary.insert(pair<string,FieldmapDescription>(Filename,FieldmapDescription(T2DMagnetoStatic,new FM2DMagnetoStatic(Filename))));
          return (*position.first).second.Map;
          break;
        case T3DDynamic:
//           position = FieldmapDictionary.insert(pair<string,FieldmapDescription>(Filename,FieldmapDescription(T3DDynamic,new FM3DDynamic(Filename))));
//           return (*position.first).second.Map;
          break;
        case T3DElectroStatic:
//           position = FieldmapDictionary.insert(pair<string,FieldmapDescription>(Filename,FieldmapDescription(T3DElectroStatic,new FM3DElectroStatic(Filename))));
//           return (*position.first).second.Map;
          break;
        case T3DMagnetoStatic:
//           position = FieldmapDictionary.insert(pair<string,FieldmapDescription>(Filename,FieldmapDescription(T3DMagnetoStatic,new FM3DMagnetoStatic(Filename))));
//           return (*position.first).second.Map;
          break;
        default:
          cerr << "could not determine type of fieldmap" << endl;
          return NULL;
        }
    }
}

void Fieldmap::deleteFieldmap(string Filename)
{
  map<string,FieldmapDescription>::iterator position = FieldmapDictionary.find(Filename);
  if (position != FieldmapDictionary.end())
    {
      if ((*position).second.RefCounter > 1)
        (*position).second.RefCounter--;
      else
        {
          delete (*position).second.Map;
          FieldmapDictionary.erase(position);
        }
    }
}

MapType Fieldmap::readHeader(string Filename)
{
  int tmpInt;
  string tmpString;
  ifstream File(Filename.c_str());
  if (!File.good())
    {
      cerr << "could not open file " << Filename << endl;
      return UNKNOWN;
    }
  
//   File >> tmpInt;
//   File.close();
  
//   return MapType(tmpInt);
  File >> tmpString;
  File.close();

  if (tmpString == "2DMagnetoStatic")
    return T2DMagnetoStatic;
  else if (tmpString == "2DElectroStatic")
    return T2DElectroStatic;
  else if (tmpString == "2DDynamic")
    return T2DDynamic;
  else if (tmpString == "1DDynamic")
    return T1DDynamic;
  else
    {
      cerr << "unknow fieldmap type" << endl;
      return UNKNOWN;
    }
}

map<string,Fieldmap::FieldmapDescription> Fieldmap::FieldmapDictionary = map<string,Fieldmap::FieldmapDescription>();
