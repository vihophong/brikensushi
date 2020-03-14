#define tree_cxx
#include "tree.cpp"
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <map>

typedef struct {
  Int_t x,y,z;
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

  TFile* fin= TFile::Open(infile);
  char tmpchar[1000];
  sprintf(tmpchar,"group1");
  tree* inbeta=new tree(fin,tmpchar);
  inbeta->SetID(5);
  sprintf(tmpchar,"group2");    
  tree* inion=new tree(fin,tmpchar);
  inion->SetID(4);


  //! set configuration file
  inbeta->ReadConfigTable("dssdconfiglowe.txt");

   Long64_t nentriesbeta = inbeta->fChain->GetEntriesFast();
   Long64_t nentriesion = inion->fChain->GetEntriesFast();

   cout<<nentriesbeta<<"-"<<nentriesion<<endl;
   //! output
   //AIDAClass* wasabi=new AIDAClass;
   AIDAClass* wasabi=new AIDAClass;
   datatype2* treedata=new datatype2;
   TFile* fout = new TFile(outfile,"RECREATE");
   TTree* tree = new TTree("aida","aida");
   //tree->Branch("aida",&wasabi);
   tree->Branch("x",&treedata->x,"x/I");
   tree->Branch("y",&treedata->y,"y/I");
   tree->Branch("z",&treedata->z,"z/I");
   tree->Branch("ex",&treedata->ex,"ex/D");
   tree->Branch("ey",&treedata->ey,"ey/D");
   tree->Branch("adcx",&treedata->adcx,"adcx/D");
   tree->Branch("adcy",&treedata->adcy,"adcy/D");


   std::multimap <unsigned long long,AIDAClass*> datamap; //! sort by timestamp
   std::multimap <unsigned long long,AIDAClass*>::iterator datamap_it; //! sort by timestamp

   Long64_t nbytes = 0, nb = 0;

   //! ion entries
   for (Long64_t jentry=0; jentry<nentriesion;jentry++) {
       Long64_t ientry = inion->LoadTree(jentry);
       if (ientry < 0) break;
       nb = inion->fChain->GetEntry(jentry);   nbytes += nb;
       inion->Reconstruction();
       //! get nz
       Int_t nz[N_DSSD];
       for (int ii=0;ii<N_DSSD;ii++) nz[ii]=0;
       for (unsigned short i=0;i<inion->GetNRecoData();i++){
           if (inion->GetRecoData(i)->EX>0&&inion->GetRecoData(i)->EY>0){
               nz[(Int_t)inion->GetRecoData(i)->z]++;
           }
       }
       for (unsigned short i=0;i<inion->GetNRecoData();i++){
           if (inion->GetRecoData(i)->EX>0&&inion->GetRecoData(i)->EY>0){
               AIDAClass* data=new AIDAClass;
               syncrecodata(data,inion->GetRecoData(i));
               data->nz=nz[(Int_t)data->z];
               datamap.insert(std::make_pair(data->T,data));
           }
       }
       inion->ClearRecoData();
   }

   //! beta entries
   for (Long64_t jentry=0; jentry<nentriesbeta;jentry++) {
       Long64_t ientry = inbeta->LoadTree(jentry);
       if (ientry < 0) break;
       nb = inbeta->fChain->GetEntry(jentry);   nbytes += nb;
       inbeta->Reconstruction();

       //! get nz
       Int_t nz[N_DSSD];
       for (int ii=0;ii<N_DSSD;ii++) nz[ii]=0;
       for (unsigned short i=0;i<inbeta->GetNRecoData();i++){
           if (inbeta->GetRecoData(i)->EX>0&&inbeta->GetRecoData(i)->EY>0){
               nz[(Int_t)inbeta->GetRecoData(i)->z]++;
           }
       }
       for (unsigned short i=0;i<inbeta->GetNRecoData();i++){
           if (inbeta->GetRecoData(i)->EX>0&&inbeta->GetRecoData(i)->EY>0){
               AIDAClass* data=new AIDAClass;
               syncrecodata(data,inbeta->GetRecoData(i));
               data->nz=nz[(Int_t)data->z];
               datamap.insert(std::make_pair(data->T,data));
           }
       }
       inbeta->ClearRecoData();
   }

   //! fill time-ordered data to tree
   for(datamap_it=datamap.begin();datamap_it!=datamap.end();datamap_it++){
       AIDAClass* hit=(AIDAClass*)datamap_it->second;
       syncrecodata(wasabi,hit);
       if (wasabi->nx==1&&wasabi->ny==1){
	 treedata->x=wasabi->x;
	 treedata->y=wasabi->y;
	 treedata->z=wasabi->z;
	 treedata->ex=wasabi->EX;
	 treedata->ey=wasabi->EY;
	 treedata->adcx=wasabi->EX;
	 treedata->adcy=wasabi->EY;
	 tree->Fill();
       }
   }
   tree->Write();
   fout->Close();
}
