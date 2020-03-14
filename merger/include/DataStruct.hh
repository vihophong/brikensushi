#ifndef DATASTRUCT_H
#define DATASTRUCT_H

#include "TTree.h"
const Int_t kMaxGamma = 500;
const Int_t kMaxNeutron = 500;

typedef struct {
    ULong64_t evt;
    ULong64_t ts; 	 //timestamp in ns
    Double_t t,x,y,ex,ey,ion_x,ion_y,ion_ex,ion_ey,zet,aoq,beta,deltaxy;//ts diff in ns
    Short_t z,ion_z,multx,multy,multz,ndecay,isbump;

    Int_t gc_hit;
    Double_t gc_E[kMaxGamma];
    Double_t gc_T[kMaxGamma];//gamma time in ns
    Int_t gc_ch[kMaxGamma];

    Int_t neu_hit;
    Double_t neu_E[kMaxNeutron];
    Double_t neu_T[kMaxNeutron];//moderation time in us
    Int_t neu_ch[kMaxNeutron];
    Double_t neu_x[kMaxNeutron];
    Double_t neu_y[kMaxNeutron];
    Int_t neub_hit;
    Double_t neub_E[kMaxNeutron];
    Double_t neub_T[kMaxNeutron];//moderation time in us
    Int_t neub_ch[kMaxNeutron];
    Double_t neub_x[kMaxNeutron];
    Double_t neub_y[kMaxNeutron];
} datatype;

typedef struct {
    Double_t yso_t,yso_e,yso_x,yso_y;
    Double_t zet,aoq,beta;
    Double_t F11L_T,F11L_E,F11R_T,F11R_E,F7_T,veto_T,veto_E;

    Int_t gc_hit;
    Double_t gc_E[kMaxGamma];
    Double_t gc_T[kMaxGamma];//gamma time in ns
    Double_t gc_Tslew[kMaxGamma];//gamma time in ns
    Int_t gc_ch[kMaxGamma];

    Int_t gc1_hit;
    Double_t gc1_E[kMaxGamma];
    Double_t gc1_T[kMaxGamma];//gamma time in ns
    Double_t gc1_Tslew[kMaxGamma];//gamma time in ns
    Int_t gc1_ch[kMaxGamma];

    Int_t gc2_hit;
    Double_t gc2_E[kMaxGamma];
    Double_t gc2_T[kMaxGamma];//gamma time in ns
    Double_t gc2_Tslew[kMaxGamma];//gamma time in ns
    Int_t gc2_ch[kMaxGamma];

    Int_t ab1_hit;
    Double_t ab1_E[kMaxGamma];
    Double_t ab1_T[kMaxGamma];//gamma time in ns
    Double_t ab1_Tslew[kMaxGamma];//gamma time in ns
    Int_t ab1_ch[kMaxGamma];//first hit channel
    Short_t ab1_mult[kMaxGamma];//multiplicity

    Int_t ab2_hit;
    Double_t ab2_E[kMaxGamma];
    Double_t ab2_T[kMaxGamma];//gamma time in ns
    Double_t ab2_Tslew[kMaxGamma];//gamma time in ns
    Int_t ab2_ch[kMaxGamma];//first hit channel
    Short_t ab2_mult[kMaxGamma];//multiplicity

    Int_t neu_hit;
    Double_t neu_E[kMaxNeutron];
    Double_t neu_T[kMaxNeutron];//moderation time in us
    Int_t neu_ch[kMaxNeutron];

} datatypeisomer;


typedef struct{
    Double_t gc_E;
    Long64_t gc_T;//gamma time in ns
    Double_t gc_Tslew;//gamma time in ns
    Int_t gc_ch;
} gammahit;

typedef struct{
    Double_t ab_E;
    Long64_t ab_T;//gamma time in ns
    Double_t ab_Tslew;//gamma time in ns
    Int_t ab_ch;//first hit channel
    Short_t ab_mult[4];//multiplicity
} gammaab;

typedef struct{
    int correntrybrips;
    int correntryf11r;
    int correntryf11l;
    int correntrydEtop;
    int correntrydEbot;
    int correntryvetodown;
    Double_t yso_t,yso_e,yso_x,yso_y;
    std::vector<gammahit*> gammagc1_vector;
    std::vector<gammahit*> gammagc2_vector;
    std::vector<gammaab*> gammaab1_vector;
    std::vector<gammaab*> gammaab2_vector;
}ImplantCorrelationVector;

class wasabiHit : public TObject{
public:
    wasabiHit(){}
    virtual ~wasabiHit(){}
    void Clear(Option_t *option = ""){
      ts = -9999;
      x = -9999;
      y = -9999;
      ex = -9999;
      ey = -9999;
      blxmax = -9999;
      blymax = -9999;
      nx = -9999;
      ny = -9999;
      nz = -9999;
      zmax = -9999;
      z = -9999;
    }
    virtual void Copy(wasabiHit& obj){
        obj.ts=ts;
        obj.x=x;
        obj.y=y;
        obj.z=z;
        obj.ex=ex;
        obj.ey=ey;
        obj.blxmax=blxmax;
        obj.blymax=blymax;
        obj.nx=nx;
        obj.ny=ny;
        obj.nz=nz;
        obj.zmax=zmax;
        obj.z=z;
    }
    Long64_t ts;
    Double_t x;
    Double_t y;
    Double_t ex;
    Double_t ey;
    Double_t blxmax;
    Double_t blymax;
    Short_t nx;
    Short_t ny;
    Short_t nz;
    Short_t zmax;
    Int_t z;
    /// \cond CLASSIMP
    ClassDef(wasabiHit,1);
    /// \endcond
} ;

#endif
