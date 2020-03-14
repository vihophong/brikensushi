#define tree_cxx
#include "tree.cpp"
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <map>

typedef struct {
  Int_t x,y,z,nz;
  Double_t ex,ey,adcx,adcy;
} datatype2;


void sorter(char* infile, char* outfile)
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

    char tmpchar[1000];
    sprintf(tmpchar,"group1");
    tree* inbeta=new tree(infile,tmpchar);
    inbeta->SetID(5);
    sprintf(tmpchar,"group2");
    tree* inion=new tree(infile,tmpchar);
    inion->SetID(4);


    //! set configuration file
    inbeta->ReadConfigTable("dssdconfiglowe.txt");
    inion->ReadConfigTable("dssdconfighighe.txt");

   Long64_t nentriesbeta = inbeta->fChain->GetEntries();
   Long64_t nentriesion = inion->fChain->GetEntries();

   cout<<nentriesbeta<<"-"<<nentriesion<<endl;
   //! output
   //AIDAClass* wasabi=new AIDAClass;
   TFile* fout = new TFile(outfile,"RECREATE");

   Double_t ex[N_DSSD];
   Double_t ey[N_DSSD];
   Double_t x[N_DSSD];
   Double_t y[N_DSSD];
   Int_t zi[N_DSSD];
   for (Int_t i=0;i<N_DSSD;i++) zi[i]=i;
   Int_t nz;

   TTree* tree = new TTree("beta","beta");
   tree->Branch("ex",ex,Form("ex[%d]/D",N_DSSD));
   tree->Branch("ey",ey,Form("ey[%d]/D",N_DSSD));
   tree->Branch("x",x,Form("x[%d]/D",N_DSSD));
   tree->Branch("y",y,Form("y[%d]/D",N_DSSD));
   tree->Branch("zi",zi,Form("zi[%d]/I",N_DSSD));
   tree->Branch("nz",&nz,"nz/I");

   TTree* treeion = new TTree("ion","ion");
   treeion->Branch("ex",ex,Form("ex[%d]/D",N_DSSD));
   treeion->Branch("ey",ey,Form("ey[%d]/D",N_DSSD));
   treeion->Branch("x",x,Form("x[%d]/D",N_DSSD));
   treeion->Branch("y",y,Form("y[%d]/D",N_DSSD));
   treeion->Branch("zi",zi,Form("zi[%d]/I",N_DSSD));
   treeion->Branch("nz",&nz,"nz/I");

   Long64_t nbytes = 0, nb = 0;
   //! ion entries
   for (Long64_t jentry=0; jentry<1;jentry++) {
       Long64_t ientry = inion->LoadTree(jentry);
       if (ientry < 0) break;
       nb = inion->fChain->GetEntry(jentry);   nbytes += nb;
       inion->Reconstruction();
       if(jentry%5000 == 0)
            cout << (Double_t) jentry/(Double_t) nentriesion * 100  << " % ion entries proceeded \r "<<flush;
       memset(ex,0,sizeof(ex));
       memset(ey,0,sizeof(ey));
       memset(x,0,sizeof(x));
       memset(y,0,sizeof(y));
       nz=0;
       for (unsigned short i=0;i<inion->GetNRecoData();i++){
           if (inion->GetRecoData(i)->EX>0&&inion->GetRecoData(i)->EY>0){
               AIDAClass* data=inion->GetRecoData(i);
               if (data->nx==1&&data->ny==1){//select single pixel events
                   ex[(Int_t)data->z]=data->EX;
                   ey[(Int_t)data->z]=data->EY;
                   x[(Int_t)data->z]=data->x;
                   y[(Int_t)data->z]=data->y;
                   nz=data->nz;
               }
           }
       }
       treeion->Fill();
       inion->ClearRecoData();
   }

   //! beta entries
   for (Long64_t jentry=0; jentry<nentriesbeta;jentry++) {
       Long64_t ientry = inbeta->LoadTree(jentry);
       if (ientry < 0) break;
       nb = inbeta->fChain->GetEntry(jentry);   nbytes += nb;
       inbeta->Reconstruction();
       if(jentry%5000 == 0)
            cout << (Double_t) jentry/(Double_t) nentriesbeta * 100  << " % beta entriesproceeded \r "<<flush;
       memset(ex,0,sizeof(ex));
       memset(ey,0,sizeof(ey));
       memset(x,0,sizeof(x));
       memset(y,0,sizeof(y));
       nz=0;
       for (unsigned short i=0;i<inbeta->GetNRecoData();i++){
           if (inbeta->GetRecoData(i)->EX>0&&inbeta->GetRecoData(i)->EY>0){
               AIDAClass* data=inbeta->GetRecoData(i);
               if (data->nx==1&&data->ny==1){//select single pixel events
                   ex[(Int_t)data->z]=data->EX;
                   ey[(Int_t)data->z]=data->EY;
                   x[(Int_t)data->z]=data->x;
                   y[(Int_t)data->z]=data->y;
                   nz=data->nz;
               }
           }
       }
       tree->Fill();
       inbeta->ClearRecoData();
   }
   tree->Write();
   treeion->Write();
   fout->Close();
}
