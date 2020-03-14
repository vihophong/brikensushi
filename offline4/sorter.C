#define tree_cxx
#include "tree.cpp"
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <map>

//#define KEEP_ORIGINAL_TREE

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
   AIDAClass* aidadata=new AIDAClass;
   TFile* fout = new TFile(outfile,"RECREATE");

   //! IFIC merger format
   TTree* treeaida = new TTree("aida","aida");
   treeaida->Branch("aida",&aidadata);

   //! Our format
   TTree* treebeta = new TTree("wbeta","wbeta");
   WASABIClass *wasabidatabeta=inbeta->GetWASABIData();
   treebeta->Branch("timestamp",&wasabidatabeta->timestamp,"timestamp/L");
   treebeta->Branch("nhits",&wasabidatabeta->nhits,"nhits/I");
   treebeta->Branch("xcog",wasabidatabeta->xcog,"xcog[nhits]/D");
   treebeta->Branch("ycog",wasabidatabeta->ycog,"ycog[nhits]/D");
   treebeta->Branch("xmax",wasabidatabeta->xmax,"xmax[nhits]/D");
   treebeta->Branch("ymax",wasabidatabeta->ymax,"ymax[nhits]/D");
   treebeta->Branch("ex",wasabidatabeta->ex,"ex[nhits]/D");
   treebeta->Branch("ey",wasabidatabeta->ey,"ey[nhits]/D");
   treebeta->Branch("z",wasabidatabeta->z,"z[nhits]/I");
   treebeta->Branch("blxmax",wasabidatabeta->blxmax,"blxmax[nhits]/D");
   treebeta->Branch("blymax",wasabidatabeta->blxmax,"blymax[nhits]/D");
   treebeta->Branch("nx",wasabidatabeta->nx,"nx[nhits]/I");
   treebeta->Branch("ny",wasabidatabeta->ny,"ny[nhits]/I");
   treebeta->Branch("zmax",&wasabidatabeta->zmax,"zmax/I");
   treebeta->Branch("nz",&wasabidatabeta->nz,"nz/I");
   treebeta->Branch("dssd_adc",inbeta->dssd_adc,Form("dssd_adc[%d]/D",N_DSSD*(N_STRIP_X+N_STRIP_Y)));
   treebeta->Branch("dssd_e",inbeta->dssd_e,Form("dssd_e[%d]/D",N_DSSD*(N_STRIP_X+N_STRIP_Y)));
   treebeta->Branch("dssd_bl",inbeta->dssd_bl,Form("dssd_bl[%d]/D",N_DSSD*(N_STRIP_X+N_STRIP_Y)));
   treebeta->Branch("dssd_ch",inbeta->dssd_ch,Form("dssd_ch[%d]/I",N_DSSD*(N_STRIP_X+N_STRIP_Y)));

#ifdef KEEP_ORIGINAL_TREE
   treebeta->Branch("evt",&inbeta->evt,Form("evt[%d]/I",V1740_N_CH*V1740_N_BOARD_G1));
   treebeta->Branch("dgtz_e",inbeta->dgtz_e,Form("dgtz_e[%d]/D",V1740_N_CH*V1740_N_BOARD_G1));
   treebeta->Branch("dgtz_bl",inbeta->dgtz_bl,Form("dgtz_bl[%d]/D",V1740_N_CH*V1740_N_BOARD_G1));
   treebeta->Branch("dgtz_ch",inbeta->dgtz_ch,Form("dgtz_ch[%d]/I",V1740_N_CH*V1740_N_BOARD_G1));
   treebeta->Branch("dgtz_nsample",inbeta->dgtz_nsample,Form("dgtz_nsample[%d]/s",V1740_N_CH*V1740_N_BOARD_G1));
   treebeta->Branch("dgtz_ts",inbeta->dgtz_ts,Form("dgtz_ts[%d]/l",V1740_N_CH*V1740_N_BOARD_G1));
   treebeta->Branch("dgtz_waveform",inbeta->dgtz_waveform,Form("dgtz_waveform[%d][%d]/s",V1740_N_CH*V1740_N_BOARD_G1,MAX_N_SAMPLE));
   treebeta->Branch("dgtz_sample",inbeta->dgtz_sample,Form("dgtz_sample[%d][%d]/s",V1740_N_CH*V1740_N_BOARD_G1,MAX_N_SAMPLE));
#endif

   TTree* treeion = new TTree("wion","wion");
   WASABIClass *wasabidataion=inion->GetWASABIData();
   treeion->Branch("timestamp",&wasabidataion->timestamp,"timestamp/L");
   treeion->Branch("nhits",&wasabidataion->nhits,"nhits/I");
   treeion->Branch("xcog",wasabidataion->xcog,"xcog[nhits]/D");
   treeion->Branch("ycog",wasabidataion->ycog,"ycog[nhits]/D");
   treeion->Branch("xmax",wasabidataion->xmax,"xmax[nhits]/D");
   treeion->Branch("ymax",wasabidataion->ymax,"ymax[nhits]/D");
   treeion->Branch("ex",wasabidataion->ex,"ex[nhits]/D");
   treeion->Branch("ey",wasabidataion->ey,"ey[nhits]/D");
   treeion->Branch("z",wasabidataion->z,"z[nhits]/I");
   treeion->Branch("blxmax",wasabidataion->blxmax,"blxmax[nhits]/D");
   treeion->Branch("blymax",wasabidataion->blxmax,"blymax[nhits]/D");
   treeion->Branch("nx",wasabidataion->nx,"nx[nhits]/I");
   treeion->Branch("ny",wasabidataion->ny,"ny[nhits]/I");
   treeion->Branch("zmax",&wasabidataion->zmax,"zmax/I");
   treeion->Branch("nz",&wasabidataion->nz,"nz/I");
   treeion->Branch("dssd_adc",inion->dssd_adc,Form("dssd_adc[%d]/D",N_DSSD*(N_STRIP_X+N_STRIP_Y)));
   treeion->Branch("dssd_e",inion->dssd_e,Form("dssd_e[%d]/D",N_DSSD*(N_STRIP_X+N_STRIP_Y)));
   treeion->Branch("dssd_bl",inion->dssd_bl,Form("dssd_bl[%d]/D",N_DSSD*(N_STRIP_X+N_STRIP_Y)));
   treeion->Branch("dssd_ch",inion->dssd_ch,Form("dssd_ch[%d]/I",N_DSSD*(N_STRIP_X+N_STRIP_Y)));

#ifdef KEEP_ORIGINAL_TREE
   treeion->Branch("evt",&inion->evt,Form("evt[%d]/I",V1740_N_CH*V1740_N_BOARD_G1));
   treeion->Branch("dgtz_e",inion->dgtz_e,Form("dgtz_e[%d]/D",V1740_N_CH*V1740_N_BOARD_G1));
   treeion->Branch("dgtz_bl",inion->dgtz_bl,Form("dgtz_bl[%d]/D",V1740_N_CH*V1740_N_BOARD_G1));
   treeion->Branch("dgtz_ch",inion->dgtz_ch,Form("dgtz_ch[%d]/I",V1740_N_CH*V1740_N_BOARD_G1));
   treeion->Branch("dgtz_nsample",inion->dgtz_nsample,Form("dgtz_nsample[%d]/s",V1740_N_CH*V1740_N_BOARD_G1));
   treeion->Branch("dgtz_ts",inion->dgtz_ts,Form("dgtz_ts[%d]/l",V1740_N_CH*V1740_N_BOARD_G1));
   treeion->Branch("dgtz_waveform",inion->dgtz_waveform,Form("dgtz_waveform[%d][%d]/s",V1740_N_CH*V1740_N_BOARD_G1,MAX_N_SAMPLE));
   treeion->Branch("dgtz_sample",inion->dgtz_sample,Form("dgtz_sample[%d][%d]/s",V1740_N_CH*V1740_N_BOARD_G1,MAX_N_SAMPLE));
#endif


   std::multimap <unsigned long long,AIDAClass*> datamap; //! sort by timestamp
   std::multimap <unsigned long long,AIDAClass*>::iterator datamap_it; //! sort by timestamp

   Long64_t nbytes = 0, nb = 0;

   //! ion entries
   for (Long64_t jentry=0; jentry<nentriesion;jentry++) {
       Long64_t ientry = inion->LoadTree(jentry);
       if (ientry < 0) break;
       nb = inion->fChain->GetEntry(jentry);   nbytes += nb;
       inion->ClearRecoData();
       inion->Reconstruction();

       //! stuffs for making aidalike data
       for (unsigned short i=0;i<inion->GetNRecoData();i++){
           if (inion->GetRecoData(i)->EX>0&&inion->GetRecoData(i)->EY>0){
             if (inion->GetWASABIData()->zmax==(Int_t)inion->GetRecoData(i)->z){
                   AIDAClass* data=new AIDAClass;
                   syncrecodata(data,inion->GetRecoData(i));
		   //data->nz=nz[(Int_t)data->z];
                   data->nz=inion->GetWASABIData()->nz;
                   datamap.insert(std::make_pair(data->T,data));
	     }
           }
       }
       treeion->Fill();
   }

   //! beta entries
   for (Long64_t jentry=0; jentry<nentriesbeta;jentry++) {
       Long64_t ientry = inbeta->LoadTree(jentry);
       if (ientry < 0) break;
       nb = inbeta->fChain->GetEntry(jentry);   nbytes += nb;
       inbeta->ClearRecoData();
       inbeta->Reconstruction();

       //! stuffs for making aidalike data
       for (unsigned short i=0;i<inbeta->GetNRecoData();i++){
           if (inbeta->GetRecoData(i)->EX>0&&inbeta->GetRecoData(i)->EY>0){
               AIDAClass* data=new AIDAClass;
               syncrecodata(data,inbeta->GetRecoData(i));
               //data->nz=nz[(Int_t)data->z];
               data->nz=inbeta->GetWASABIData()->nz;
               datamap.insert(std::make_pair(data->T,data));
           }
       }
       treebeta->Fill();
   }

   //! fill time-ordered data to tree
   for(datamap_it=datamap.begin();datamap_it!=datamap.end();datamap_it++){
       AIDAClass* hit=(AIDAClass*)datamap_it->second;       
       syncrecodata(aidadata,hit);
       treeaida->Fill();
   }
   treebeta->Write();
   treeion->Write();
   treeaida->Write();
   fout->Close();
}
