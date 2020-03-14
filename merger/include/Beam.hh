#ifndef __BEAM_H
#define __BEAM_H
#include <iostream>
#include <vector>
#include <cstdlib>
#include <math.h>

#include "TObject.h"

#define kMaxBeamInfo 1
#define kMaxRIPSInfo 1
#define kMaxTOFInfo 1



//using namespace std;
/*!
  Container for the full beam, tof, beta and pid information
*/

struct TreeData {
ULong64_t ts;
ULong64_t sts;
Double_t tof;
Double_t zet;
Double_t aoq;
Double_t f5x;
Double_t f11x;
Double_t f11y;
Double_t f11dt;
Double_t beta;
};

class Beam : public TObject {
public:
  //! default constructor  
  Beam(){
    Clear();
  }
  virtual ~Beam(){}
  //! Clear all information
  void Clear(Option_t *option = ""){
    fts = 0;
    memset(faoq,0.,sizeof(faoq));
    memset(fzet,0.,sizeof(fzet));
    memset(ftof,0.,sizeof(ftof));
    memset(fbeta,0.,sizeof(fbeta));
    memset(fdelta,0.,sizeof(fdelta));
    memset(fbrho,0,sizeof(fbrho));
  }
  virtual void Copy(Beam& obj){
      obj.SetTimestamp(fts);
      for (Int_t i=0;i<kMaxBeamInfo;i++){
          obj.SetAQ(i,faoq[i]);
          obj.SetZ(i,fzet[i]);
          obj.SetBeta(i,fbeta[i]);
	  obj.SetBrho(i,fbrho[i]);
      }
      for (Int_t i=0;i<kMaxTOFInfo;i++){
          obj.SetTOF(i,ftof[i]);
      }
      for (Int_t i=0;i<kMaxRIPSInfo;i++){
          obj.SetDelta(i,fdelta[i]);
      }
  }

  //! Set the timestamp
  void SetTimestamp(unsigned long long ts){fts = ts;}

  //! Set the A/Q ratio
  void SetAQ(unsigned short j, double aoq){
    if( j>kMaxBeamInfo-1) return;
    faoq[j] = aoq;
  }

  //! Set the Z number
  void SetZ(unsigned short j, double zet){
    if( j>kMaxBeamInfo-1) return;
    fzet[j] = zet;
  }
  //! Set both A/Q and Z
  void SetAQZ(unsigned short j, double aoq, double zet){
    if( j>kMaxBeamInfo-1) return;
    faoq[j] = aoq;
    fzet[j] = zet;
  }
  //! Set the time-of-flight
  void SetTOF(unsigned short j, double tof){
    if( j>kMaxTOFInfo-1) return;
    ftof[j] = tof;
  }
  //! Set the beta
  void SetBrho(unsigned short j, double brho){
    if( j>kMaxBeamInfo-1) return;
    fbrho[j] = brho;
  }
  //! Set the beta
  void SetBeta(unsigned short j, double beta){
    if( j>kMaxBeamInfo-1) return;
    fbeta[j] = beta;
  }

  //! Set the delta
  void SetDelta(unsigned short j, double delta){
    if( j>kMaxRIPSInfo-1) return;
    fdelta[j] = delta;
  }

  //! Get timestamp
  unsigned long long GetTimestamp(){return fts;}

  //! Get the A/Q ratio
  double GetAQ(unsigned short j){
    if( j>kMaxBeamInfo-1) return sqrt(-1.);
    return faoq[j];
  }
  //! Get the Z number
  double GetZ(unsigned short j){
    if( j>kMaxBeamInfo-1) return sqrt(-1.);
    return fzet[j];
  }
  //! Get the time-of-flight
  double GetTOF(unsigned short j){
    if( j>2) return sqrt(-1.);
    return ftof[j];
  }
  //! Get Brho
  double GetBrho(unsigned short j){
    if( j>2) return sqrt(-1.);
    return fbrho[j];
  }
  //! Get beta
  double GetBeta(unsigned short j){
    if( j>2) return sqrt(-1.);
    return fbeta[j];
  }
  //! Get Delta
  double GetDelta(unsigned short j){
    if( j>2) return sqrt(-1.);
    return fdelta[j];
  }

protected:
  unsigned long long fts;
  //! A/Q for 3-5, 5-7, 3-7
  double faoq[kMaxBeamInfo];
  //! Z for 3-5, 5-7, 3-7
  double fzet[kMaxBeamInfo];

  //! brho for 3-5, 5-7, 3-7
  double fbrho[kMaxBeamInfo];

  //! time-of-flight for 3-5, 5-7, 3-7
  double ftof[kMaxTOFInfo];
  //! beta for 3-7, 8-11, 7-8
  double fbeta[kMaxBeamInfo];
  //! delta momentum 3-5, 5-7, 3-7
  double fdelta[kMaxRIPSInfo];

  /// \cond CLASSIMP
  ClassDef(Beam,1);
  /// \endcond
};


#endif
