
#include <iostream>
#include <pmonitor/pmonitor.h>
#include "myanalysis.h"

#include "libDataStruct.h"
#include "dpp.h"

#include <TH1.h>
#include <TH2.h>

#define V1740_HDR 6
#define V1740_N_CH 64
#define V1740_N_BOARD 4
#define V1740_PACKET_0 100

#define DGTZ_CLK_RES 10
#define V1730_EVENT_TYPE 1
#define V1740_EVENT_TYPE 2
#define TDC_EVENT_TYPE 3


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


int init_done = 0;

using namespace std;

TH1F *h1;
TH2F *h2;


NIGIRI* data;

int pinit()
{
  if (init_done) return 1;
  init_done = 1;
  // initialize data
  data = new NIGIRI;

  h1 = new TH1F ( "h1","test histogram", 500, 0, 500);
  h2 = new TH2F ( "h2","test histogram 2D", 500, 0, 500,500,0,4000);

  return 0;
}

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

          if (data->b==0){
              NIGIRIHit* chdata=data->GetHit(32);
              int itcnt=0;
              for (std::vector<UShort_t>::iterator it = chdata->pulse.begin() ; it != chdata->pulse.end(); ++it){
                  if (itcnt<N_MAX_WF_LENGTH){
                      h1->SetBinContent(itcnt+1,*it);
                      h2->Fill(itcnt,*it);
                  }
                  itcnt++;
              }
          }
          //data->Print();
          delete p1740[bb];
      }//if there is packet data

  }//loop through all posible packet data
  return 0;
}

