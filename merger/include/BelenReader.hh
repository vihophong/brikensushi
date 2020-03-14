#ifndef BELENREADER_H
#define BELENREADER_H

#include <vector>
#include "TTree.h"
#include "TClonesArray.h"
#include "BELEN.hh"
#include "Clover.hh"
#include "BRIKENTreeData.hh"

#include <fstream>
#include "TRandom.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TH2F.h"
#include "TEllipse.h"
#include "TColor.h"
#include "TROOT.h"
#include "TLatex.h"

#include <string>

#define MaxID 1000
#define MaxIndex1 30
#define MaxIndex2 10
#define belenClockResolution 20 //in ns

class BelenReader
{
public:
    BelenReader();
    virtual ~BelenReader();

    //! set geometry mapping files
    void SetGeoMapping(char* mapf){fmappingfile = mapf;}
    //! Get geometry mapping (currently being called by Init)
    void GetGeoMapping();
    //! initialization
    void Init(char* belenfile);
    //! book tree if we want to store files
    void BookTree(TTree* treeNeutron, TTree *treeGamma, TTree *treeAnc, BELENHit* neutron,CloverHit* gamma, BELENHit* anc);

    //! book tree if we want to store files
    void BookYSOTree(TTree* treeYSOion, TTree *treeYSObeta);

    //! close file reader
    void CloseReader(){finfile->Close();}

    //! get next event (any: gamma or neutron or ancillary)
    bool GetNextEvent();
    //! get next neutron event
    bool GetNextNeutronEvent();
     //! get next gamma event
    bool GetNextGammaEvent();
    //! get next ancillary event
    bool GetNextAncEvent();

    //! get next YSO event
    bool GetNextYSOEvent();


    //! Get total number of events from the input file
    int GetNEvents(){return fnentries;}
    //! Get total number of YSO events from the input file
    int GetNYSOEvents(){return fnentriesYSO;}

    //! Get current processing entry
    int GetCurrentEvent(){return fcurentry;}

    //! Get current processing entry
    int GetCurrentYSOEvent(){return fcurentryYSO;}

    //! Get current processing neutron entry
    int GetCurrentNeutronEvent(){return fBLNeuEntry;}
    //! Get current processing gamma entry
    int GetCurrentGammaEvent(){return fBLGamEntry;}
    //! Get current processing ancillary entry
    int GetCurrentAncEvent(){return fBLAncEntry;}

    //! Get current processing YSO ion event entry
    int GetCurrentYSOIonEvent(){return fBLYSOEntryIon;}
    //! Get current processing YSO ion event entry
    int GetCurrentYSOBetaEvent(){return fBLYSOEntryBeta;}

    //! neutron hits
    BELENHit* GetNeutron(){return flocalNeutron;}
    //! gamma hits
    CloverHit* GetGamma(){return flocalGamma;}
    //! ancillary hist
    BELENHit* GetAnc(){return flocalAnc;}

    //! ancillary hits
    std::vector<BELENHit*> GetAncAIDAPL(){return flocalAncAIDAPL;}
    std::vector<BELENHit*> GetAncUpstreamPL(){return flocalAncUpstreamPL;}
    std::vector<BELENHit*> GetAncdE(){return flocalAncdE;}
    std::vector<BELENHit*> GetAncF11PL(){return flocalAncF11PL;}

    //! Get event type, return 1 if neutron, 2 if gamma and 3 if ancillary!
    Double_t GetEnergy(){return fE;}
    ULong64_t GetTimestamp(){return fT;}
    UShort_t GetId(){return fId;}
    UShort_t GetType(){return ftype;}
    UShort_t GetIndex1(){return fIndex1;}
    UShort_t GetIndex2(){return fIndex2;}
    UShort_t GetInfoFlag(){return fInfoFlag;}
    std::string GetName(){return fName;}


    //! clear Ancillary hits and free memory!
    void ClearAncHits();

    //! special function from MC simulation!
    void PerturbateHe3(UShort_t He3Id);
    void PerturbateClover(UShort_t Index1,UShort_t Index2);

protected:
    //! mapping file
    char* fmappingfile;
    //! number of entries in belen file
    unsigned long long fnentries;

    //! number of entries in belen file for YSO
    unsigned long long fnentriesYSO;

    //! current entry in belen file
    unsigned long long fcurentry;

    //! current YSO entry in belen file
    unsigned long long fcurentryYSO;

    //! current neutron entry
    unsigned long long fBLNeuEntry;
    //! current gamma entry
    unsigned long long fBLGamEntry;
    //! current ancillary  entry
    unsigned long long fBLAncEntry;


    //! current YSO  entry
    unsigned long long fBLYSOEntryIon;
    unsigned long long fBLYSOEntryBeta;

    //! current time stamp of gamma
    unsigned long long fBLtsGamma;
    //! current time stamp of neutron
    unsigned long long fBLtsNeutron;
    //! current time stamp of ancillary
    unsigned long long fBLtsAnc;

    //! current YSO  entry
    unsigned long long fBLtsYSO;

    //! data structure
    Double_t        fE;
    ULong64_t       fT;
    UShort_t        fId;
    UShort_t        ftype;
    UShort_t        fIndex1;
    UShort_t        fIndex2;
    UShort_t        fInfoFlag;
    std::string          fName;
    Double_t        fposX;
    Double_t        fposY;
    Double_t        fposZ;


    Double_t fHe3Ecal[MaxID][2];

    Double_t fHe3Id2posX[MaxID];
    Double_t fHe3Id2posY[MaxID];
    Double_t fHe3Id2posZ[MaxID];
    Double_t fHe3Id2diameter[MaxID];
    UShort_t fHe3Id2ring[MaxID];
    Double_t fHe3Id2length[MaxID];

    Double_t fCrystalId2posX[MaxIndex1][MaxIndex2];
    Double_t fCrystalId2posY[MaxIndex1][MaxIndex2];
    Double_t fCrystalId2posZ[MaxIndex1][MaxIndex2];

    //! neutron hits
    BELENHit* flocalNeutron;
    //! gamma hits
    CloverHit* flocalGamma;
    //! ancillary hits
    std::vector<BELENHit*> flocalAncUpstreamPL;
    std::vector<BELENHit*> flocalAncAIDAPL;
    std::vector<BELENHit*> flocalAncdE;
    std::vector<BELENHit*> flocalAncF11PL;
    //! ancillary single hit
    BELENHit* flocalAnc;

    //! input files
    TFile* finfile;

    //! a tree of input root files data
    BrikenTreeData* ftreedataNeuron;
    BrikenTreeData* ftreedataGamma;
    BrikenTreeData* ftreedataAnc;

    YSOData* ftreedataYSO;

    //! input tree
    TTree* ftree;

    //! input tree YSO
    TTree* ftreeYSO;

    //!  output tree for Neutron
    TTree* fmtrNeutron;
    //!  output tree for Gamma
    TTree* fmtrGamma;
    //!  output tree for Anc
    TTree* fmtrAnc;

    //!  output tree for YSO
    TTree* fmtrYSObeta;
    TTree* fmtrYSOion;

    TRandom rr;

    bool fflag_filldata;

};

#endif // BELENREADER_H
