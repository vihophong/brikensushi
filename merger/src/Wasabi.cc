#include "Wasabi.hh"

Wasabi::Wasabi(char *inputFileName,Bool_t isBeta, TTree *tree) : fChain(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(inputFileName);
      if (!f || !f->IsOpen()) {
         f = new TFile(inputFileName);
      }
      if (isBeta)
        f->GetObject("wbeta",tree);
      else{
        f->GetObject("wion",tree);
      }
   }
   Init(tree);
}

Wasabi::~Wasabi()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Wasabi::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

Long64_t Wasabi::GetEntries()
{
    if (!fChain) return 0;
    return fChain->GetEntries();
}

Long64_t Wasabi::LoadTree(Long64_t entry)
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

void Wasabi::Init(TTree *tree)
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

   fChain->SetBranchAddress("timestamp", &timestamp, &b_timestamp);
   fChain->SetBranchAddress("nhits", &nhits, &b_nhits);
   fChain->SetBranchAddress("xcog", xcog, &b_xcog);
   fChain->SetBranchAddress("ycog", ycog, &b_ycog);
   fChain->SetBranchAddress("xmax", xmax, &b_xmax);
   fChain->SetBranchAddress("ymax", ymax, &b_ymax);
   fChain->SetBranchAddress("ex", ex, &b_ex);
   fChain->SetBranchAddress("ey", ey, &b_ey);
   fChain->SetBranchAddress("z", z, &b_z);
   fChain->SetBranchAddress("blxmax", blxmax, &b_blxmax);
   fChain->SetBranchAddress("blymax", blymax, &b_blymax);
   fChain->SetBranchAddress("nx", nx, &b_nx);
   fChain->SetBranchAddress("ny", ny, &b_ny);
   fChain->SetBranchAddress("zmax", &zmax, &b_zmax);
   fChain->SetBranchAddress("nz", &nz, &b_nz);
   fChain->SetBranchAddress("dssd_adc", dssd_adc, &b_dssd_adc);
   fChain->SetBranchAddress("dssd_e", dssd_e, &b_dssd_e);
   fChain->SetBranchAddress("dssd_bl", dssd_bl, &b_dssd_bl);
   fChain->SetBranchAddress("dssd_ch", dssd_ch, &b_dssd_ch);
   Notify();
}

Bool_t Wasabi::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.
   return kTRUE;
}

void Wasabi::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Wasabi::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
