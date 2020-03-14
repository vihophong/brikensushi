//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Feb 24 18:27:28 2020 by ROOT version 6.14/00
// from TTree wbeta/wbeta
// found on file: beamrun53_wrun4646to4679.root
//////////////////////////////////////////////////////////

#ifndef Wasabi_h
#define Wasabi_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>


#include "Wasabidefs.hh"

// Header file for the classes stored in the TTree if any.

class Wasabi {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Long64_t        timestamp;
   Int_t           nhits;
   Double_t        xcog[MAX_N_HITS];   //[nhits]
   Double_t        ycog[MAX_N_HITS];   //[nhits]
   Double_t        xmax[MAX_N_HITS];   //[nhits]
   Double_t        ymax[MAX_N_HITS];   //[nhits]
   Double_t        ex[MAX_N_HITS];   //[nhits]
   Double_t        ey[MAX_N_HITS];   //[nhits]
   Int_t           z[MAX_N_HITS];   //[nhits]
   Double_t        blxmax[MAX_N_HITS];   //[nhits]
   Double_t        blymax[MAX_N_HITS];   //[nhits]
   Int_t           nx[MAX_N_HITS];   //[nhits]
   Int_t           ny[MAX_N_HITS];   //[nhits]
   Int_t           zmax;
   Int_t           nz;
   Double_t        dssd_adc[MAX_N_STRIPS];
   Double_t        dssd_e[MAX_N_STRIPS];
   Double_t        dssd_bl[MAX_N_STRIPS];
   Int_t           dssd_ch[MAX_N_STRIPS];

   // List of branches
   TBranch        *b_timestamp;   //!
   TBranch        *b_nhits;   //!
   TBranch        *b_xcog;   //!
   TBranch        *b_ycog;   //!
   TBranch        *b_xmax;   //!
   TBranch        *b_ymax;   //!
   TBranch        *b_ex;   //!
   TBranch        *b_ey;   //!
   TBranch        *b_z;   //!
   TBranch        *b_blxmax;   //!
   TBranch        *b_blymax;   //!
   TBranch        *b_nx;   //!
   TBranch        *b_ny;   //!
   TBranch        *b_zmax;   //!
   TBranch        *b_nz;   //!
   TBranch        *b_dssd_adc;   //!
   TBranch        *b_dssd_e;   //!
   TBranch        *b_dssd_bl;   //!
   TBranch        *b_dssd_ch;   //!

   Wasabi(char* inputFileName,Bool_t isBeta, TTree *tree=0);
   virtual ~Wasabi();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   Long64_t    GetEntries();
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

