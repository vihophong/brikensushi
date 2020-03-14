#ifndef MERGER_H
#define MERGER_H 1

#include "BELEN.hh"
#include "Clover.hh"
#include "Beam.hh"

#include "Wasabi.hh"

#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH1F.h"
#include "TCut.h"
#include "TCutG.h"
#include "TLatex.h"
#include "TSpectrum.h"
#include "TF1.h"
#include "fstream"
#include "DataStruct.hh"


//#define POS_CORR_CENTER_OF_GRAVITY 1



#define MaxNRI 1000

class Merger
{
public:
    Merger();
    virtual ~Merger();
    void SetWasabiFile(char* wasabifile){finputWasabi = wasabifile;}
    void SetBigripsFile(char* bigripsfile){finputBigrips = bigripsfile;}
    void SetBrikenFile(char* brikenfile){finputBriken = brikenfile;}
    void Init();

    void ReadPID(char* pidfile, Int_t ncutpts=20);
    void ReadWasabi(Int_t opt=2);
    void ReadBigrips();
    void ReadBRIKEN(unsigned int startN=0, unsigned int stopN=0,unsigned int startG=0, unsigned int stopG=0,unsigned int startA=0, unsigned int stopA=0);

    void BookIonBetaTree();
    void BookImplantTree();

    Int_t GetNri(){return nri;}
    TTree* GetTreeRI(Int_t i){if (i<0) return ftreeallRI; else return ftreeRI[i];}
    TTree* GetTreeImpRI(Int_t i){if (i<0) return ftreeimplantAll;return ftreeimplantRI[i];}
    TCutG* GetCUTRI(Int_t i){return cutg[i];}

    void ResetImplantData();

    void MergeImplant();

    void MergeIonBeta();

    void MakeF11NeutronVeto();

    TH1F* getProbeHisto1D(Int_t ihisto){return h1dProbe[ihisto];}


    void copyWasabiHit(wasabiHit* source,wasabiHit* destination){
        destination->ts=source->ts;
        destination->x=source->x;
        destination->y=source->x;
        destination->ex=source->ex;
        destination->ey=source->ey;
        destination->blxmax=source->blxmax;
        destination->blymax=source->blymax;
        destination->nx=source->nx;
        destination->ny=source->ny;
        destination->nx=source->nx;
        destination->zmax=source->zmax;
        destination->z=source->z;
    }

    //!stuff for addback
    void CopyAddbackData(gammaab* ab_src,gammaab* ab_des)
    {
        ab_des->ab_ch=ab_src->ab_ch;
        ab_des->ab_E=ab_src->ab_E;
        ab_des->ab_T=ab_src->ab_T;
        ab_des->ab_Tslew=ab_src->ab_Tslew;
        //for (Int_t i=0;i<4;i++) ab_des->ab_mult[i]=ab_src->ab_mult[i];
        memcpy(ab_des->ab_mult,ab_src->ab_mult,sizeof(Short_t)*4);
    }
    void DoAddback();

    //!stuff for slew correction
    void ReadSlewCorr();

protected:

    //! file name to be read
    char* finputWasabi;
    char* finputBigrips;
    char* finputBriken;

    //! file to be read
    TFile* fBigripsFile;
    TFile* fBrikenFile;

    //! tree to be read
    TTree* ftrBigrips;
    TTree* ftrNeutron;
    TTree* ftrGamma;
    TTree* ftrAnc;

    //! enclosing data (tree,file, filename) of Wasabi
    Wasabi * fWasabiBeta;
    Wasabi * fWasabiIon;

    //! number of enries in each tree
    Long64_t fnentriesWasabiBeta;
    Long64_t fnentriesWasabiIon;
    Long64_t fnentriesBigrips;
    Long64_t fnentriesNeutron;
    Long64_t fnentriesGamma;
    Long64_t fnentriesAnc;
    Long64_t fnentriesF11Anc;

    //! data read from stream
    TreeData* fbigrips;
    CloverHit* fclover;
    BELENHit* fneutron;
    BELENHit* fanc;

    //! time maps
    std::multimap < Long64_t,  wasabiHit*> fwasabiIonMap;
    std::multimap < Long64_t,  wasabiHit*>::iterator fwasabiIonMap_it;

    std::multimap < Long64_t,  std::pair<ImplantCorrelationVector*,wasabiHit*>> fwasabiImplantMap;
    std::multimap < Long64_t,  std::pair<ImplantCorrelationVector*,wasabiHit*>>::iterator fwasabiImplantMap_it;


    std::multimap < Long64_t,  wasabiHit*> fwasabiBetaMap;
    std::multimap < Long64_t,  wasabiHit*>::iterator fwasabiBetaMap_it;


    std::multimap < Long64_t, unsigned int> fbigripsMap;
    std::multimap < Long64_t, unsigned int>::iterator fbigripsMap_it;
    std::multimap < Long64_t, BELENHit*> fhe3Map;
    std::multimap < Long64_t, BELENHit*>::iterator fhe3Map_it;
    std::multimap < Long64_t, unsigned int> fcloverMap;
    std::multimap < Long64_t, unsigned int>::iterator fcloverMap_it;

    std::multimap < Long64_t, gammaab*> faddbackclover1Map;
    std::multimap < Long64_t, gammaab*> faddbackclover2Map;
    std::multimap < Long64_t, gammaab*>::iterator faddbackclover1Map_it;
    std::multimap < Long64_t, gammaab*>::iterator faddbackclover2Map_it;

    std::multimap < Long64_t, unsigned int> fancMap;
    std::multimap < Long64_t, unsigned int>::iterator fancMap_it;
    std::multimap < Long64_t, BELENHit*> fdtpulserMap;
    std::multimap < Long64_t, BELENHit*>::iterator fdtpulserMap_it;

    std::multimap < Long64_t, unsigned int> fcorrimpMap;
    std::multimap < Long64_t, unsigned int>::iterator fcorrimpMap_it;

    std::multimap < Long64_t, unsigned int> ff11ancMap;
    std::multimap < Long64_t, unsigned int> ff11ancMap_it;

    std::multimap < Long64_t, unsigned int> fF11LMap;
    std::multimap < Long64_t, unsigned int> fF11RMap;
    std::multimap < Long64_t, unsigned int> fF11LRMap;
    std::multimap < Long64_t, unsigned int>::iterator fF11MapL_it;
    std::multimap < Long64_t, unsigned int>::iterator fF11MapR_it;
    std::multimap < Long64_t, unsigned int>::iterator fF11MapLR_it;

    std::multimap < Long64_t, unsigned int> fVetoTopMap;
    std::multimap < Long64_t, unsigned int> fVetoBotMap;
    std::multimap < Long64_t, unsigned int> fVetoDownMap;
    std::multimap < Long64_t, unsigned int> fvetoMap;
    std::multimap < Long64_t, unsigned int>::iterator fVetoTopMap_it;
    std::multimap < Long64_t, unsigned int>::iterator fVetoBotMap_it;
    std::multimap < Long64_t, unsigned int>::iterator fVetoDownMap_it;
    std::multimap < Long64_t, unsigned int>::iterator fvetoMap_it;

    std::multimap < Long64_t, unsigned int> fdETopMap;
    std::multimap < Long64_t, unsigned int> fdEBotMap;
    std::multimap < Long64_t, unsigned int>::iterator fdETopMap_it;
    std::multimap < Long64_t, unsigned int>::iterator fdEBotMap_it;

    //! time windows
    Long64_t fTW_IonBetalow;
    Long64_t fTW_IonBetaup;

    Long64_t fTW_IonPidlow;
    Long64_t fTW_IonPidup;
    Long64_t fTW_IonF11LRlow;
    Long64_t fTW_IonF11LRup;

    Long64_t fTW_IondElow;
    Long64_t fTW_IondEup;
    Long64_t fTW_IonDownVetolow;
    Long64_t fTW_IonDownVetoup;


    Long64_t fTW_PIDGammalow;
    Long64_t fTW_PIDGammaup;


    //! pulser stuff
    Long64_t ftsbeginpulser;
    Long64_t ftsendpulser;
    TH1F* fhpulser[141];
    TH1F* fhpulserall[141];
    Double_t ftotaltimepulser;
    TH2F* fh2deadtime;
    TH1F* fh1deadtime;
    TH2F* fh2deadtime2;
    TH1F* fh1deadtime2;
    TH2F* fh2deadtime3;
    TH1F* fh1deadtime3;
    TH2F* fh2deadtime4;
    TH1F* fh1deadtime4;
    TH1F* fh1dtpulser;

    //! Spatial windows
    Long64_t fsquareSpatialWindowX;
    Long64_t fsquareSpatialWindowY;

    //! stuff for PID separation
    Int_t nri;
    Int_t nbinszet,nbinsaoq;
    Double_t zetrange[2];
    Double_t aoqrange[2];
    Int_t enablepid[MaxNRI];
    Int_t enablepid2[MaxNRI];
    TString nameri[MaxNRI];
    TString latexnametri[MaxNRI];
    Double_t parmsri[MaxNRI][7];
    Double_t halflife[MaxNRI];
    TCutG* cutg[MaxNRI];
    TLatex* pidtag[MaxNRI];

    //! Oututs
    TH1F* h1dProbe[100];

    wasabiHit* fdataOutWasabiBeta;
    wasabiHit* fdataOutWasabiIon;

    datatypeisomer implantdata;

    TTree* ftreeallRI;
    TTree* ftreeRI[MaxNRI];
    TTree* ftreeimplantAll;
    TTree* ftreeimplantRI[MaxNRI];

    //!stuff for slew correction
    Bool_t isslewcorr;
    Double_t a[8];
    Double_t b[8];
    Double_t c[8];
    Double_t d[8];
    //! gamma calibration provided by JJ
    Double_t fsep[8];
    Double_t flow_offset[8];
    Double_t flow_gain[8];
    Double_t flow_se[8];
    Double_t fhigh_offset[8];
    Double_t fhigh_gain[8];
    Double_t fhigh_se[8];
    Double_t fcgainold[8];
    Double_t fcoffsetold[8];

};

#endif
