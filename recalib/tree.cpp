//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Oct  5 09:52:51 2017 by ROOT version 5.34/36
// from TTree tree/tree
// found on file: ../test.root
//////////////////////////////////////////////////////////
#include "tree.h"
// default constant
tree::tree(TFile* infilef,char* treename,TTree *tree) : fChain(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
     TFile *f = infilef;     
     f->GetObject(treename,tree);

   }
   Init(tree);

   for (int i=0;i<N_DSSD*(N_STRIP_X+N_STRIP_Y);i++){
       dssd_thr[i]=CH_THRESHOLD;
       dssd_cal0[i]=0.;
       dssd_cal1[i]=1.;
   }
   //! default mapping
   Int_t strip=0;
   for (int z=0;z<N_DSSD;z++){
       for (int xy=0;xy<N_STRIP_X+N_STRIP_Y;xy++){
           strip=z*(N_STRIP_X+N_STRIP_Y)+xy;
           if (xy<N_STRIP_X){
               striptoch[strip]=z*N_STRIP_X+xy;
           }else{
	      if ((xy-N_STRIP_X)<8){// correct mapping
		striptoch[strip]=z*N_STRIP_Y+7-xy+N_STRIP_X+V1740_N_CH;
	      }else{
		striptoch[strip]=z*N_STRIP_Y+xy-N_STRIP_X+V1740_N_CH;
	      }
              //striptoch[strip]=z*N_STRIP_Y+xy-N_STRIP_X+V1740_N_CH;
           }
	   //cout<<strip<<"-"<<striptoch[strip]<<endl;
       }
   }
   ClearRecoData();

   fid=0;
   flag_mapping=false;
}

tree::~tree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

void tree::ReadConfigTable(char* configfile)
{

  ifstream inpf(configfile);
  if (inpf.fail()){
    cout<<"No Configuration table is given"<<endl;
    return;
  }

  cout<<"Start reading configuration table: "<<configfile<<endl;
  Int_t index,dssd_index,strip_index,dgtz_channel;
  Double_t threshold,cal0,cal1;
  Int_t mm=0;

  while (inpf.good()){
    inpf>>index>>dssd_index>>strip_index>>dgtz_channel>>threshold>>cal0>>cal1;
    cout<<index<<"-"<<dssd_index<<"-"<<strip_index<<"-"<<dgtz_channel<<"-"<<threshold<<"-"<<cal0<<"-"<<cal1<<endl;
    dssd_cal0[index]=cal0;
    dssd_cal1[index]=cal1;
    striptoch[index]=dgtz_channel;
    dssd_thr[index]=threshold;
    mm++;
  }
  cout<<"Read "<<mm<<" line"<<endl;
  inpf.close();
}

void tree::ClearRecoData(){
    memset(dssd_adc,0,sizeof(dssd_adc));
    memset(dssd_e,0,sizeof(dssd_e));
    nrecodata=0;
    //maybe not nescessary
    for (size_t idx=0;idx<recodata.size();idx++){
        delete recodata[idx];
    }
    recodata.clear();
}

Int_t tree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t tree::LoadTree(Long64_t entry)
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

void tree::Init(TTree *tree)
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

   fChain->SetBranchAddress("evt", evt, &b_evt);
   fChain->SetBranchAddress("dgtz_e", dgtz_e, &b_dgtz_e);
   fChain->SetBranchAddress("dgtz_bl", dgtz_bl, &b_dgtz_bl);
   fChain->SetBranchAddress("dgtz_ch", dgtz_ch, &b_dgtz_ch);
   fChain->SetBranchAddress("dgtz_nsample", dgtz_nsample, &b_dgtz_nsample);
   fChain->SetBranchAddress("dgtz_ts", dgtz_ts, &b_dgtz_ts);
   fChain->SetBranchAddress("dgtz_waveform", dgtz_waveform, &b_dgtz_waveform);
   fChain->SetBranchAddress("dgtz_sample", dgtz_sample, &b_dgtz_sample);
   Notify();
}

Bool_t tree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void tree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t tree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}


void tree::GetPosMax()
{
    for (Int_t z=0;z<N_DSSD;z++){
        AIDAClass* data=new AIDAClass;
        Double_t maxex=0;
        Double_t maxey=0;
        Int_t xmaxex=-9999;
        Int_t ymaxey=-9999;
        Double_t ex=0;
        Double_t ey=0;
        Int_t nx=0;
        Int_t ny=0;
        Int_t xstrminch=z*(N_STRIP_X+N_STRIP_Y);
        Int_t xstrmaxch=z*(N_STRIP_X+N_STRIP_Y)+N_STRIP_X;
        Int_t ystrminch=z*(N_STRIP_X+N_STRIP_Y)+N_STRIP_X;
        Int_t ystrmaxch=z*(N_STRIP_X+N_STRIP_Y)+N_STRIP_X+N_STRIP_Y;
        //cout<<xstrminch<<"-"<<xstrmaxch<<"-"<<ystrminch<<"-"<<ystrmaxch<<endl;

        for (Int_t i=xstrminch;i<xstrmaxch;i++){
            if (dssd_e[i]>dssd_thr[i]){
                if (dssd_e[i]>maxex){
                  maxex=dssd_e[i];
                  xmaxex=i-xstrminch;
                }
                nx++;
                ex+=dssd_e[i];
            }
        }
        for (Int_t i=ystrminch;i<ystrmaxch;i++){
            if (dssd_e[i]>dssd_thr[i]){
                if (dssd_e[i]>maxey){
                  maxey=dssd_e[i];
                  ymaxey=i-ystrminch;
                }
                ny++;
                ey+=dssd_e[i];
            }
        }

        data->T=dgtz_ts[0]*10;
        data->Tfast=0;
        data->x=xmaxex;
        data->y=ymaxey;
        data->EX=ex;
        data->EY=ey;
        data->E=(ex+ey)/2;
        data->z=z;
        data->nx=nx;
        data->ny=ny;
        data->nz=0;
        data->ID=fid;
        AddRecoData(data);
    }
}
void tree::Reconstruction(){
    //if (!flag_mapping) return;
    //! mapping and calibration
    for (Int_t i=0;i<N_DSSD*(N_STRIP_X+N_STRIP_Y);i++){
        dssd_adc[i]=dgtz_e[striptoch[i]];
        dssd_e[i]=dssd_adc[i]*dssd_cal1[i]+dssd_cal0[i];
    }
    GetPosMax();
}
