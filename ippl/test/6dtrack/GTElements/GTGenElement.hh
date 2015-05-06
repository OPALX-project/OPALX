#ifndef GENELEM_HH
#define GENELEM_HH
class GenElem {

public:
  GenElem()
  { 
  };

  ~GenElem() 
  { 
  };

  void basicInit()
  {
  };

  void setBMultipoleField(BMultipoleField field) {
    field_m = field;
  }

private:
   BMultipoleField field_m;
};
#endif
