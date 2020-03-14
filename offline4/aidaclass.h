#ifndef AIDACLASS_H
#define AIDACLASS_H
#include "TObject.h"

class AIDAClass {
 public:
  //! default constructor
  AIDAClass(){}
  virtual ~AIDAClass(){}
  unsigned long long T;        // Calibrated time
  unsigned long long Tfast;
  double E;    // Energy
  double EX;
  double EY;
  double x,y,z;// number of pixel for AIDA, or number of tube for BELEN
  int nx, ny, nz;// multiplicity of hits in x and y strips, and dssd planes
  unsigned char ID;    // Detector type (BigRips, Aida ion, AIDA beta, BELEN, Clovers)
  //** other stuff pending to define **//
  ClassDef(AIDAClass,1);
};
#endif // AIDACLASS_H
