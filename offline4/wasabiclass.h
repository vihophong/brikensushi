#ifndef WASABICLASS_H
#define WASABICLASS_H
#include "TObject.h"

#define kmaxhits 10

class WASABIClass {
 public:
  //! default constructor
  WASABIClass(){}
  virtual ~WASABIClass(){}
  Long64_t timestamp;        // Calibrated time
  Double_t xcog[kmaxhits]; // x position determined by central of gravity method
  Double_t ycog[kmaxhits]; // y position determined by central of gravity method
  Double_t xmax[kmaxhits]; // x position determined by position with max E method
  Double_t ymax[kmaxhits]; // y position determined by position with max E method
  Double_t ex[kmaxhits];
  Double_t ey[kmaxhits];
  Int_t z[kmaxhits]; // z position
  Double_t blxmax[kmaxhits];
  Double_t blymax[kmaxhits];
  Int_t nx[kmaxhits];
  Int_t ny[kmaxhits];

  Int_t zmax;
  Int_t nhits;
  Int_t nz;

  void Clear(){
      timestamp=0;
      memset(xcog,0,sizeof(xcog));
      memset(ycog,0,sizeof(ycog));
      memset(xmax,0,sizeof(xmax));
      memset(ymax,0,sizeof(ymax));
      memset(ex,0,sizeof(ex));
      memset(ey,0,sizeof(ey));
      memset(z,0,sizeof(z));
      memset(blxmax,0,sizeof(blxmax));
      memset(blymax,0,sizeof(blymax));
      memset(nx,0,sizeof(nx));
      memset(ny,0,sizeof(ny));
      nhits=0;
      nz=0;
      zmax=-1;
  }

  ClassDef(WASABIClass,1);
};
#endif // WASABICLASS_H
