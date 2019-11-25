#define tree_cxx
#include "tree.cpp"
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <map>

//Fitting function from TPC GET shaper
Double_t ShaperF_GET(Double_t *x, Double_t *p) {

   Double_t semiGaus = p[1]*20*TMath::Exp(-3.0*(x[0] - p[2])/p[3])*sin((x[0] - p[2])/p[3])*pow((x[0] - p[2])/p[3], p[4]);

   return (x[0] >= p[2]) ? p[0] + x[0]*p[5] + semiGaus : p[0] + x[0]*p[5];
}
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



   Int_t npulse=0;

   TH1F* hsigma=new TH1F("sig","sig",200,0,100);
   TH1F* hmaxe=new TH1F("hmaxe","hmaxe",1000,0,2000);
   TH1F* hfite=new TH1F("hfite","hfite",1000,0,2000);

   TH2F* hevse=new TH2F("evse","evse",500,0,1000,500,0,1000);
   //! beta entries
   TCanvas* c1=new TCanvas("c1","c1",900,700);
   for (Long64_t jentry=0; jentry<nentriesbeta;jentry++) {
       Long64_t ientry = inbeta->LoadTree(jentry);
       if (ientry < 0) break;
       nb = inbeta->fChain->GetEntry(jentry);   nbytes += nb;
       inbeta->Reconstruction();

       //! get nz
       //Int_t nz[N_DSSD];
       //for (int ii=0;ii<N_DSSD;ii++) nz[ii]=0;
       Int_t nz=0;
       for (unsigned short i=0;i<inbeta->GetNRecoData();i++){
           if (inbeta->GetRecoData(i)->EX>0&&inbeta->GetRecoData(i)->EY>0){
             nz++;
           }
       }

       //if (nz==4) {//select light particles events
           for (Int_t i=57;i<58;i++) {//loop through all channels
               if (inbeta->dgtz_e[i]>200) {
                   TH1F* h1=new TH1F("h1","h1",120,0,120);
                   for (Int_t j=0;j<inbeta->dgtz_nsample[i];j++) {
                       h1->Fill(j+1,inbeta->dgtz_waveform[i][j]);
                   }


                   if (h1->GetMaximumBin()<110){
		     Int_t min=0;
		     Int_t max=120;
		     if (h1->GetMaximumBin()-30>0) min=h1->GetMaximumBin()-31;
		     if (h1->GetMaximumBin()+30<=120) max=h1->GetMaximumBin()+29;

		     TF1* fWaveform = new TF1("fWaveform", ShaperF_GET, min, max, 6);
		     fWaveform->SetParNames("offset", "amplitude", "peakBeginAt", "sigma", "power", "p2");
                     fWaveform->SetParameters(450, 500, 25, 33., 3, 0.2);
                     fWaveform->SetParLimits(3, 20,70);
		     fWaveform->SetParLimits(1, 0, 10000);


                     h1->Fit(fWaveform,"LQER+","goff");

//                     for (Int_t i=0;i<6;i++){
//                         cout<<fWaveform->GetParameter(i)<<"\t";
//                     }
//                     cout<<endl;
		     //cout<<h1->GetFunction("fWaveform")->GetParameter(1)<<"\t"<<h1->GetFunction("fWaveform")->GetMaximum()<<"\t"<<h1->GetFunction("fWaveform")->GetChisquare()/h1->GetFunction("fWaveform")->GetNDF()<<endl;
                     //cout<<h1->GetFunction("fWaveform")->GetMaximum()<<"\t"<<inbeta->dgtz_bl[i]<<"\t"<<h1->GetFunction("fWaveform")->GetMaximum()-inbeta->dgtz_bl[i]<<"\t"<<inbeta->dgtz_e[i]<<endl;
		     hevse->Fill(inbeta->dgtz_e[i],h1->GetFunction("fWaveform")->GetMaximum()-inbeta->dgtz_bl[i]);

                     if (i==57) {
                         hmaxe->Fill(inbeta->dgtz_e[i]);
                         hfite->Fill(h1->GetFunction("fWaveform")->GetMaximum()-inbeta->dgtz_bl[i]);
                     }

                     h1->Draw();
                     c1->Update();
                     c1->WaitPrimitive();
                     gSystem->Sleep(2000);
		     delete fWaveform;
                   }
                   delete h1;

               }
           }
       //}
       inbeta->ClearRecoData();
       //if (jentry>1000) break;
   }
   hevse->Write();

   hmaxe->Write();
   hfite->Write();

   hsigma->Write();

   fout->Close();
}
