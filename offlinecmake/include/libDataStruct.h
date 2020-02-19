#ifndef __libDataStruct_H
#define __libDataStruct_H

#include "TObject.h"
#include "TROOT.h"
#include <vector>
#include <iostream>

#define TDC_N_CHANNEL 64
#define TDC_MAX_MULT 3

using namespace std;
/*!
  Container for the full beam, tof, beta and pid information
*/
class TDCHit : public TObject {
public:
  //! default constructor
  TDCHit(){
      Clear();
  }
  virtual ~TDCHit(){}
  void Clear(){
      ch = -9999;
      t = -9999;
  }
    UShort_t ch ;//channel number
    Int_t t;//baseline
  /// \cond CLASSIMP
  ClassDef(TDCHit,1);
  /// \endcond
};

class NIGIRIHit : public TObject {
public:
  //! default constructor
  NIGIRIHit(){
    Clear();
  }
  virtual ~NIGIRIHit(){}


  virtual void Copy(NIGIRIHit& obj){
        obj.ch=ch;
        obj.finets=finets;
        obj.cshort=cshort;
        obj.clong=clong;
        obj.baseline=baseline;
        obj.nsample=nsample;
        obj.pulse=pulse;
  }
  void Clear(){
        ch = -1;
        finets  = 0;
        cshort = 0;
        clong  = 0;
        baseline = 0;
        nsample = 0;
        pulse.clear();
  }
  void Print(){
      cout<<"ch = "<<ch<<endl;
      cout<<"finets = "<<finets<<endl;
      cout<<"clong = "<<clong<<endl;
      cout<<"pulse 0 = "<<pulse[0]<<endl;
  }
    Short_t ch ;//channel number
    Double_t finets;//finets
    Double_t cshort;//charge short
    Double_t clong;//charge long
    Double_t baseline;//baseline
    Short_t nsample;
    std::vector<UShort_t> pulse;//pulse

  /// \cond CLASSIMP
  ClassDef(NIGIRIHit,1);
  /// \endcond
};


class NIGIRI : public TObject {
public:
  //! default constructor
  NIGIRI(){      
      Clear();
  }
  virtual ~NIGIRI(){}
  void Clear(){
      losttrigger = 0;
      overrange = 0 ;
      ts = 0;
      evt_type = -1;
      evt = 0;
      b = -1;
      fmult = 0;
      for (size_t idx=0;idx<fhits.size();idx++){
          delete fhits[idx];
      }      
      fhits.clear();
  }

  virtual void Copy(NIGIRI& obj){
      for (vector<NIGIRIHit*>::iterator hitin_it=fhits.begin(); hitin_it!=fhits.end(); hitin_it++){
          NIGIRIHit* clonehit = new NIGIRIHit;
          NIGIRIHit* originhit = *hitin_it;
          originhit->Copy(*clonehit);
          obj.AddHit(clonehit);
      }
      obj.fmult=fmult;
      obj.evt_type=evt_type;
      obj.overrange=overrange;
      obj.losttrigger=losttrigger;
      obj.evt=evt;
      obj.b=b;
      obj.ts=ts;
  }
  void Print(){
      cout<<"ts = "<<ts<<endl;
      cout<<"evt = "<<evt<<endl;
      cout<<"fhits size = "<<fhits.size()<<endl;

  }
  Int_t GetMult(){return fmult;}
  NIGIRIHit* GetHit(unsigned short n){return fhits.at(n);}
  void AddHit(NIGIRIHit* hit){
	fmult++;
    fhits.push_back(hit);
  }
  Int_t fmult;
  //! common stuff
  Char_t evt_type;
  Char_t overrange;
  Char_t losttrigger;//trigger lost flag (=1: lost)
  Int_t evt;//evt number
  Short_t b;//board number
  ULong64_t ts;//timestamp
  std::vector<NIGIRIHit*> fhits;
  /// \cond CLASSIMP
  ClassDef(NIGIRI,1);
  /// \endcond
};


#endif


