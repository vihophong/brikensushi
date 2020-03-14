#ifndef __libDataStruct_H
#define __libDataStruct_H

#include "TObject.h"
#include "TROOT.h"
#include <vector>
#include <iostream>

#include <TChain.h>
#include <TFile.h>
#include <fstream>

#define V1740_HDR 6
#define V1740_N_CH 64
#define V1740_N_BOARD 4
#define V1740_PACKET_0 100
#define V1740_PACKET_1 101
#define V1740_PACKET_2 102
#define V1740_PACKET_3 103
#define V1740_N_MAX_CH 64
#define NSBL 8
#define N_REFESH_WF 1000
#define N_MAX_WF_LENGTH 500

#define COUNTS_TO_DISPLAY 10000
#define COUNTS_TO_CAL 100
#define MAX_N_SAMPLE 300

#define MAX_MAP_LENGTH 50


#define TDC_N_CHANNEL 64
#define TDC_MAX_MULT 3

// default constant
#define CH_THRESHOLD 50
//DSSDs constant
#define N_DSSD 4
#define N_STRIP_X 16
#define N_STRIP_Y 16

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

class WASABIStruct : public TObject {
public:
  //! default constructor
  WASABIStruct(){
    Clear();
  }
  virtual ~WASABIStruct(){}

  //! Clear the information
  virtual void Clear(){
      fid=0.;
  }

  //! copy cluster
  virtual void Copy(WASABIStruct& obj){

  }

  //! Printing information
  void Print(Option_t *option = "") const {

    return;
  }

protected:
  // Declaration of leaf types
  unsigned char fid;

  Int_t           evt[V1740_N_CH*V1740_N_BOARD];
  Double_t        dgtz_e[V1740_N_CH*V1740_N_BOARD];
  Double_t        dgtz_bl[V1740_N_CH*V1740_N_BOARD];
  Int_t           dgtz_ch[V1740_N_CH*V1740_N_BOARD];
  UShort_t        dgtz_nsample[V1740_N_CH*V1740_N_BOARD];
  ULong64_t       dgtz_ts[V1740_N_CH*V1740_N_BOARD];
  UShort_t        dgtz_waveform[V1740_N_CH*V1740_N_BOARD][300];
  UShort_t        dgtz_sample[V1740_N_CH*V1740_N_BOARD][300];

  Double_t           dssd_adc[N_DSSD*(N_STRIP_X+N_STRIP_Y)];
  Double_t           dssd_e[N_DSSD*(N_STRIP_X+N_STRIP_Y)];

  Double_t           dssd_thr[N_DSSD*(N_STRIP_X+N_STRIP_Y)];
  Double_t           dssd_cal0[N_DSSD*(N_STRIP_X+N_STRIP_Y)];
  Double_t           dssd_cal1[N_DSSD*(N_STRIP_X+N_STRIP_Y)];

  Int_t           striptoch[N_DSSD*(N_STRIP_X+N_STRIP_Y)];
  Bool_t          flag_mapping;

  /// \cond CLASSIMP
  ClassDef(WASABIStruct,1);
  /// \endcond
};

#endif


