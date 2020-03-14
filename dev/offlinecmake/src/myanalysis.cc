
#include <iostream>
#include <pmonitor/pmonitor.h>
#include "myanalysis.h"
#include "sorter.h"

#include "dpp.h"

#include <TH1.h>
#include <TH2.h>

#include <TFile.h>
#include <TTree.h>

#define V1740_HDR 6
#define V1740_N_CH 64
#define V1740_N_BOARD 4
#define V1740_PACKET_0 100

#define DGTZ_CLK_RES 10
#define V1730_EVENT_TYPE 1
#define V1740_EVENT_TYPE 2
#define TDC_EVENT_TYPE 3

#define V1740_N_BOARD_G1 2
#define V1740_N_BOARD_G2 2

//baseline for zero suppresion
#define NSBL 8
//baseline for making average
#define NSBL_AVE 16

//criteria for baseline determination
#define MIN_VARIATION_BASELINE 100

//max waveform lengh
#define N_MAX_WF_LENGTH 500
#define N_WF_LENGTH 500

#define LED_THRESHOLD 100
#define CFD_DELAY 1
#define CFD_FRACTION 0.5


#define MAX_MAP_LENGTH_G1 50
#define MAX_MAP_LENGTH_G2 50
#define TS_WIN_LOW_G1 100
#define TS_WIN_HIGH_G1 100
#define TS_WIN_LOW_G2 100
#define TS_WIN_HIGH_G2 100

#define MAX_N_HISTOGRAMS 100


//! maps for sorter
std::multimap <ULong64_t,NIGIRI*> datamap_lupo; //! sort by timestamp
std::multimap <ULong64_t,NIGIRI*> datamap_dgtz0; //! sort by timestamp
std::multimap <ULong64_t,NIGIRI*> datamap_dgtz1; //! sort by timestamp
std::multimap <ULong64_t,NIGIRI*> datamap_dgtz2; //! sort by timestamp
std::multimap <ULong64_t,NIGIRI*> datamap_dgtz3; //! sort by timestamp
std::multimap<ULong64_t,NIGIRI*>::iterator it_datamap_lupo;
std::multimap<ULong64_t,NIGIRI*>::iterator it_datamap_dgtz0;
std::multimap<ULong64_t,NIGIRI*>::iterator it_datamap_dgtz1;
std::multimap<ULong64_t,NIGIRI*>::iterator it_datamap_dgtz2;
std::multimap<ULong64_t,NIGIRI*>::iterator it_datamap_dgtz3;

Int_t ntotalg1;
Int_t ntotalg2;
Int_t ntotalNoCorrg1;
Int_t ntotalNoCorrg2;
Int_t ntriggerg1[V1740_N_BOARD_G1];
Int_t ntriggerg2[V1740_N_BOARD_G2];


sorter* sorterg1;
sorter* sorterg2;


int init_done = 0;


//! trees
TFile* fout;
//!group 1 data container
TTree* tree;

TH1F *h1;
TH1F *h2;


NIGIRI* data;

using namespace std;
int pinit()
{
  if (init_done) return 1;
  init_done = 1;
  // initialize data
  data = new NIGIRI;

  ntotalg1=0;
  ntotalg2=0;
  ntotalNoCorrg1=0;
  ntotalNoCorrg2=0;

  return 0;
}

// open and close tree, histograms
int phsave(const char* filename)
{
    if (filename==0){//close
        cout<<"closed, with remaining map size processed= "<<datamap_dgtz0.size()<<endl;
        mergedata(true);
        fout->cd();
        TH1F* oh1= (TH1F*)fout->Get("h1");
        oh1->Write();
        TH1F* oh2= (TH1F*)fout->Get("h2");
        oh2->Write();
        sorterg1->writeHistos();
        fout->Close();
    }else{
        fout = new TFile(filename,"RECREATE");
        h1 = new TH1F ( "h1","TS corr group 2", 500, -TS_WIN_LOW_G1, TS_WIN_HIGH_G1);
        h2 = new TH1F ( "h2","TS corr group 1", 500, -TS_WIN_LOW_G2, TS_WIN_HIGH_G2);
        ntotalg1=0;
        ntotalg2=0;
        ntotalNoCorrg1=0;
        ntotalNoCorrg2=0;
        sorterg1=new sorter();
        sorterg2=new sorter();
        sorterg1->initHistos();
        return 0;
    }
}

// process events
int process_event (Event * e)
{
  Packet* p1740[V1740_N_BOARD];

  for (int bb=0;bb<4;bb++) {
      p1740[bb] = e->getPacket(V1740_PACKET_0+bb);

      if (p1740[bb]){
          int* temp;
          int* gg;
          gg=(int*) p1740[bb]->getIntArray(temp);

          int k=V1740_HDR+V1740_N_CH;

          data->Clear();
          data->evt_type = V1740_EVENT_TYPE;
          data->b = bb;//for sorter
          data->evt = gg[2]+1;//this event start from 0
          data->overrange = (Char_t) gg[1];//intepret as channel(group) mask
          UInt_t tslsb = (UInt_t)gg[5];
          UInt_t tsmsb = (UInt_t)gg[4];
          data->ts = (((ULong64_t)tsmsb<<32)&0xFFFF00000000)|(ULong64_t)tslsb;//resolution is 16 ns!
          for (int i=0;i<V1740_N_CH;i++){
            //! header
            NIGIRIHit* chdata=new NIGIRIHit;
            chdata->ch = i;//for sorter

            int nsample = gg[i+V1740_HDR];
            chdata->nsample = nsample;
            UShort_t WaveLine[nsample];

            int ispl = 0;
            for (int j=0;j<nsample/2+nsample%2;j++){
              if (ispl<nsample) {
                  WaveLine[ispl]=gg[k]&0xFFFF;
                  chdata->pulse.push_back(gg[k]&0xFFFF);

              }
              ispl++;
              if (ispl<nsample) {
                  WaveLine[ispl]=(gg[k]>>16)&0xFFFF;
                  chdata->pulse.push_back((gg[k]>>16)&0xFFFF);
              }
              ispl++;
              k++;
            }
            if (nsample>NSBL){
                dpp *oj=new dpp(nsample,WaveLine);
                oj->baselineMean(NSBL_AVE,MIN_VARIATION_BASELINE);
                oj->makeCFD(LED_THRESHOLD,CFD_DELAY,CFD_FRACTION);
                double timeData=oj->cfdFast();
                chdata->cshort = 0;
                chdata->clong = oj->maxAdcPos(N_WF_LENGTH)-oj->bL;
                chdata->baseline = oj->bL;
                chdata->finets = timeData;
                delete oj;
            }
            data->AddHit(chdata);
          }//loop all channels


          mergedata(false,data);
          //data->Print();
          delete p1740[bb];
      }//if there is packet data

  }//loop through all posible packet data
  return 0;
}

int mergedata(bool flagend,NIGIRI*data){
    //! digitizer group 1
    if (datamap_dgtz0.size()>MAX_MAP_LENGTH_G1||flagend){
      for(it_datamap_dgtz0=datamap_dgtz0.begin();it_datamap_dgtz0!=datamap_dgtz0.end();it_datamap_dgtz0++){
        Long64_t ts=(Long64_t)it_datamap_dgtz0->first;
        NIGIRI* hit0=(NIGIRI*)it_datamap_dgtz0->second;
        Long64_t ts1 = ts - TS_WIN_LOW_G1;
        Long64_t ts2 = ts + TS_WIN_HIGH_G1;
        Long64_t corrts = 0;
        it_datamap_dgtz1 = datamap_dgtz1.lower_bound(ts1);
        NIGIRI* hit1=0;
        while(it_datamap_dgtz1!=datamap_dgtz1.end()&&it_datamap_dgtz1->first<ts2){
            corrts = (Long64_t) it_datamap_dgtz1->first;
            hit1=(NIGIRI*)it_datamap_dgtz1->second;
            if (!flagend) {
                std::multimap<ULong64_t,NIGIRI*>::iterator it_datamap_dgtz1p;
                it_datamap_dgtz1p=it_datamap_dgtz1;
                it_datamap_dgtz1p++;
                for (std::multimap<ULong64_t,NIGIRI*>::iterator it_datamaptmp=datamap_dgtz1.begin();it_datamaptmp!=it_datamap_dgtz1p;it_datamaptmp++){
                    NIGIRI* hittmp=it_datamaptmp->second;
                    hittmp->Clear();
                    delete hittmp;
                }
                datamap_dgtz1.erase(datamap_dgtz1.begin(),it_datamap_dgtz1p);
            }
            break;//only find first ts match
        }
        //! fill data here
        if (corrts!=0)
        {
            h1->Fill(ts-corrts);
            ntotalNoCorrg1++;
            sorterg1->setCurrentDgtzHits(hit0,hit1);
        }else{
            sorterg1->setCurrentDgtzHits(hit0,0);
        }
        sorterg1->doSortHighGain();
        ntotalg1++;


        //! progress report
        if (ntotalg1%1000==0&&ntotalg1>0) {
            cout<<ntotalg1<<" events in group one + "<<ntotalg2<<" events in group two are processed \r"<<endl;
            cout<<ntotalNoCorrg1<<" correlated events in group one + "<<ntotalNoCorrg2<<" correlated events in group two are processed \r"<<endl;
        }

        if (!flagend) {
            hit0->Clear();
            delete hit0;
            datamap_dgtz0.erase(datamap_dgtz0.begin(),++it_datamap_dgtz0);
        break;
        }
      }
    }

    //! digitizer group 2
    if (datamap_dgtz2.size()>MAX_MAP_LENGTH_G2||flagend){
      for(it_datamap_dgtz2=datamap_dgtz2.begin();it_datamap_dgtz2!=datamap_dgtz2.end();it_datamap_dgtz2++){
        Long64_t ts=(Long64_t)it_datamap_dgtz2->first;
        NIGIRI* hit0=(NIGIRI*)it_datamap_dgtz2->second;
        Long64_t ts1 = ts - TS_WIN_LOW_G2;
        Long64_t ts2 = ts + TS_WIN_HIGH_G2;
        Long64_t corrts = 0;
        it_datamap_dgtz3 = datamap_dgtz3.lower_bound(ts1);
        NIGIRI* hit1=0;
        while(it_datamap_dgtz3!=datamap_dgtz3.end()&&it_datamap_dgtz3->first<ts2){
            corrts = (Long64_t) it_datamap_dgtz3->first;
            hit1=(NIGIRI*)it_datamap_dgtz3->second;

            if (!flagend) {
                std::multimap<ULong64_t,NIGIRI*>::iterator it_datamap_dgtz3p;
                it_datamap_dgtz3p=it_datamap_dgtz3;
                it_datamap_dgtz3p++;
                for (std::multimap<ULong64_t,NIGIRI*>::iterator it_datamaptmp=datamap_dgtz3.begin();it_datamaptmp!=it_datamap_dgtz3p;it_datamaptmp++){
                    NIGIRI* hittmp=it_datamaptmp->second;
                    hittmp->Clear();
                    delete hittmp;
                }
                datamap_dgtz3.erase(datamap_dgtz3.begin(),it_datamap_dgtz3p);
            }
            break;//only find first ts match
        }
        //! fill data here
        if (corrts!=0)
        {
            h2->Fill(ts-corrts);
            ntotalNoCorrg2++;
            sorterg2->setCurrentDgtzHits(hit0,hit1);
        }else{
            sorterg2->setCurrentDgtzHits(hit0,0);
        }
        sorterg2->doSortLowGain();

        ntotalg2++;
        if (!flagend) {
            hit0->Clear();
            delete hit0;
            datamap_dgtz2.erase(datamap_dgtz2.begin(),++it_datamap_dgtz2);
            break;
        }
      }
    }

    if (!flagend){
        if (data->b==0){
            NIGIRI* datac=new NIGIRI;
            data->Copy(*datac);
            datamap_dgtz0.insert(make_pair(data->ts,datac));
            ntriggerg1[0]++;
        }else if(data->b==1){
            NIGIRI* datac=new NIGIRI;
            data->Copy(*datac);
            datamap_dgtz1.insert(make_pair(data->ts,datac));
            ntriggerg1[1]++;
        }else if(data->b==2){
            NIGIRI* datac=new NIGIRI;
            data->Copy(*datac);
            datamap_dgtz2.insert(make_pair(data->ts,datac));
            ntriggerg2[0]++;
        }else if(data->b==3){
            NIGIRI* datac=new NIGIRI;
            data->Copy(*datac);
            datamap_dgtz3.insert(make_pair(data->ts,datac));
            ntriggerg2[1]++;
        }
    }
}

