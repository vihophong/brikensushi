//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Feb 22 20:03:54 2020 by ROOT version 6.14/00
// from TTree beta/beta
// found on file: outrootfiles/run4089.root
//////////////////////////////////////////////////////////

#ifndef fitchannel_h
#define fitchannel_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class fitchannel {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           x;
   Int_t           y;
   Int_t           z;
   Int_t           nz;
   Double_t        ex;
   Double_t        ey;
   Double_t        adcx;
   Double_t        adcy;

   // List of branches
   TBranch        *b_x;   //!
   TBranch        *b_y;   //!
   TBranch        *b_z;   //!
   TBranch        *b_nz;   //!
   TBranch        *b_ex;   //!
   TBranch        *b_ey;   //!
   TBranch        *b_adcx;   //!
   TBranch        *b_adcy;   //!

   fitchannel(TTree *tree=0);
   virtual ~fitchannel();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(Int_t dssd, Double_t ecutminx,Double_t ecutminy, Double_t ecutmax, Int_t npoint, char* outname, Int_t startch, Int_t stopch);
   void CheckCalib(char *calibfile, char *outfilename);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef fitchannel_cxx
fitchannel::fitchannel(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("outrootfiles/allruns.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("outrootfiles/allruns.root");
      }
      f->GetObject("beta",tree);

   }
   Init(tree);
}

fitchannel::~fitchannel()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t fitchannel::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t fitchannel::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void fitchannel::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("x", &x, &b_x);
   fChain->SetBranchAddress("y", &y, &b_y);
   fChain->SetBranchAddress("z", &z, &b_z);
   fChain->SetBranchAddress("nz", &nz, &b_nz);
   fChain->SetBranchAddress("ex", &ex, &b_ex);
   fChain->SetBranchAddress("ey", &ey, &b_ey);
   fChain->SetBranchAddress("adcx", &adcx, &b_adcx);
   fChain->SetBranchAddress("adcy", &adcy, &b_adcy);
   Notify();
}

Bool_t fitchannel::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void fitchannel::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t fitchannel::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef fitchannel_cxx
