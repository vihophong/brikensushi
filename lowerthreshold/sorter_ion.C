#define tree_cxx
#include "tree.cpp"
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <map>
void sorter_ion(char* infile, char* outfile)
{

//   In a ROOT session, you can do:
//      Root > .L tree.C
//      Root > tree t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

  //TFile* fin= TFile::Open(infile);
  char tmpchar[1000];
  //sprintf(tmpchar,"group1");
  //tree* inbeta=new tree(infile,tmpchar);
  //inbeta->SetID(5);
  sprintf(tmpchar,"group2");    
  tree* inion=new tree(infile,tmpchar);
  inion->SetID(4);


  //! set configuration file
  inion->ReadConfigTable("dssdconfighighe.txt");
  //inbeta->ReadConfigTable("dssdconfiglowe.txt");

   Long64_t nentriesion = inion->fChain->GetEntries();
   //Long64_t nentriesbeta = inbeta->fChain->GetEntries();

   cout<<nentriesion<<endl;

   //! output
   TFile* fout = new TFile(outfile,"RECREATE");


   TH2F* hthrch=new TH2F("h2","h2",128,0,128,1000,0,1000);

   Long64_t nbytes = 0, nb = 0;

//   //! ion entries
//   for (Long64_t jentry=0; jentry<nentriesion;jentry++) {
//       Long64_t ientry = inion->LoadTree(jentry);
//       if (ientry < 0) break;
//       nb = inion->fChain->GetEntry(jentry);   nbytes += nb;
//       //cout<<inion->dgtz_ts[0]<<endl;
//   }
//   cout<<"finished ion entries"<<endl;

   TH2F* hthrchvsevt[128];
   Int_t ncntevtch[128];
   for (Int_t i=0;i<128;i++) {
       hthrchvsevt[i]=new TH2F(Form("hch%d",i),Form("hch%d",i),5000,0,nentriesion/20,1000,0,1000);
       ncntevtch[i]=0;
   }

   //! ion entries
   for (Long64_t jentry=0; jentry<nentriesion;jentry++) {
       Long64_t ientry = inion->LoadTree(jentry);
       if (ientry < 0) break;
       nb = inion->fChain->GetEntry(jentry);   nbytes += nb;
       if(jentry%5000 == 0)
            cout << (Double_t) jentry/(Double_t) nentriesion * 100  << " % proceeded \r "<<flush;
       inion->Reconstruction();

       //! get nz
       Int_t nzz[N_DSSD];
       for (int ii=0;ii<N_DSSD;ii++) nzz[ii]=0;
       Int_t nz=0;
       for (unsigned short i=0;i<inion->GetNRecoData();i++){
           if (inion->GetRecoData(i)->EX>0&&inion->GetRecoData(i)->EY>0){
             nz++;
             nzz[(Int_t)inion->GetRecoData(i)->z]++;
           }
       }
       for (Int_t i=0;i<128;i++) {//loop through all channels
           Int_t zz=i/(N_STRIP_X+N_STRIP_Y);
           Int_t n_over_thr=0;
           if (nzz[zz]==0){//select event with no signal
               for (Int_t j=i-2;j<i+3;j++) {//look if neighboring channel has signals
                   if (j>0&&j<128)
                       if (inion->dssd_adc[j]>inion->dssd_thr[j]) n_over_thr++;
               }
               if (n_over_thr==0){
                     hthrch->Fill(i,inion->dgtz_waveform[inion->striptoch[i]][0]);
                     hthrchvsevt[i]->Fill(ncntevtch[i],inion->dgtz_waveform[inion->striptoch[i]][0]);
                     ncntevtch[i]++;
               }
           }
       }

       inion->ClearRecoData();
   }
   cout<<"finished ion entries"<<endl;
   for (Int_t i=0;i<128;i++) hthrchvsevt[i]->Write();
   hthrch->Write();


   fout->Close();
}
