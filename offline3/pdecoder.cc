
#include "myanalysis.cc"
#include <stdint.h>
#include "libDataStruct.h"

int CFD_delay=1;//in Tclk unit, 1 Tclk=2ns
double CFD_fraction=0.5;
double LED_threshold=100;
double LED_threshold_LED=100;
double LED_threshold_740=100;
double LED_threshold_LED_740=100;
int gateOffset=10;
int shortGate=20;
int longGate=150;
int nBaseline=16;
int minVarBaseline=100; //criteria for baseline determination
int mode_selection=2;


//! dont care about the following code
uint16_t seperateint(int inint,bool islsb){
  if (islsb) return (uint16_t) (inint&0xFFFF);
  else return (uint16_t) (inint>>16);
}

NIGIRI* data;

int init_done = 0;
int pinit()
{
  if (init_done) return 1;
  init_done = 1;
  gROOT->ProcessLine(".L libDataStruct.so");
  data = new NIGIRI;
  Init();
  return 0;
}

int pclose(){
    CloseMe();
    cout<<"\n\nFinish!"<<endl;
    return 0;
}


int process_event (Event * e)
{  
  //! v1740 packet
  Packet *p1740_0=e->getPacket(V1740_PACKET_0);
  if (p1740_0)
    {
      int* temp;
      int* gg;
      gg=(int*) p1740_0->getIntArray(temp);
      int size=p1740_0->getPadding();
      //! content
      int k=V1740_HDR+V1740_N_MAX_CH;
      data->Clear();
      data->evt_type = V1740_EVENT_TYPE;
      data->b = 0;//for sorter
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
            oj->baselineMean(nBaseline,minVarBaseline);
            double timeData = 0;
            if (mode_selection==0) {
                timeData=oj->ledPos(LED_threshold_LED_740);
            }else if (mode_selection==1) {
                timeData=oj->ledWithCorr(LED_threshold_740);
            }else if (mode_selection==2){
                oj->makeCFD(LED_threshold_740,CFD_delay,CFD_fraction);
                timeData=oj->cfdFast();
            }else if (mode_selection==3){
                oj->makeCFD(LED_threshold_740,CFD_delay,CFD_fraction);
                timeData=oj->cfdSlow();
            }

            //oj->chargeInter(cShort,cLong,timeData,gateOffset,shortGate,longGate);
            chdata->cshort = 0;
            chdata->clong = oj->maxAdcPos(N_MAX_WF_LENGTH)-oj->bL;
            chdata->baseline = oj->bL;
            chdata->finets = timeData;
            delete oj;
        }
        data->AddHit(chdata);
      }
      ProcessEvent(data);
       delete p1740_0;
  }

  //! v1740 packet
  Packet *p1740_1=e->getPacket(V1740_PACKET_1);
  if (p1740_1)
    {
      int* temp;
      int* gg;
      gg=(int*) p1740_1->getIntArray(temp);
      int size=p1740_1->getPadding();
      //! content
      int k=V1740_HDR+V1740_N_MAX_CH;
      data->Clear();
      data->evt_type = V1740_EVENT_TYPE;
      data->b = 1;//for sorter
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
            oj->baselineMean(nBaseline,minVarBaseline);
            double timeData = 0;
            if (mode_selection==0) {
                timeData=oj->ledPos(LED_threshold_LED_740);
            }else if (mode_selection==1) {
                timeData=oj->ledWithCorr(LED_threshold_740);
            }else if (mode_selection==2){
                oj->makeCFD(LED_threshold_740,CFD_delay,CFD_fraction);
                timeData=oj->cfdFast();
            }else if (mode_selection==3){
                oj->makeCFD(LED_threshold_740,CFD_delay,CFD_fraction);
                timeData=oj->cfdSlow();
            }

            //oj->chargeInter(cShort,cLong,timeData,gateOffset,shortGate,longGate);
            chdata->cshort = 0;
            chdata->clong = oj->maxAdcPos(N_MAX_WF_LENGTH)-oj->bL;
            chdata->baseline = oj->bL;
            chdata->finets = timeData;
            delete oj;
        }
        data->AddHit(chdata);
      }
      ProcessEvent(data);
       delete p1740_1;
  }

  //! v1740 packet
  Packet *p1740_2=e->getPacket(V1740_PACKET_2);
  if (p1740_2)
    {
      int* temp;
      int* gg;
      gg=(int*) p1740_2->getIntArray(temp);
      int size=p1740_2->getPadding();
      //! content
      int k=V1740_HDR+V1740_N_MAX_CH;
      data->Clear();
      data->evt_type = V1740_EVENT_TYPE;
      data->b = 2;//for sorter
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
            oj->baselineMean(nBaseline,minVarBaseline);
            double timeData = 0;
            if (mode_selection==0) {
                timeData=oj->ledPos(LED_threshold_LED_740);
            }else if (mode_selection==1) {
                timeData=oj->ledWithCorr(LED_threshold_740);
            }else if (mode_selection==2){
                oj->makeCFD(LED_threshold_740,CFD_delay,CFD_fraction);
                timeData=oj->cfdFast();
            }else if (mode_selection==3){
                oj->makeCFD(LED_threshold_740,CFD_delay,CFD_fraction);
                timeData=oj->cfdSlow();
            }

            //oj->chargeInter(cShort,cLong,timeData,gateOffset,shortGate,longGate);
            chdata->cshort = 0;
            chdata->clong = oj->maxAdcPos(N_MAX_WF_LENGTH)-oj->bL;
            chdata->baseline = oj->bL;
            chdata->finets = timeData;
            delete oj;
        }
        data->AddHit(chdata);
      }
      ProcessEvent(data);
       delete p1740_2;
  }

  //! v1740 packet
  Packet *p1740_3=e->getPacket(V1740_PACKET_3);
  if (p1740_3)
    {
      int* temp;
      int* gg;
      gg=(int*) p1740_3->getIntArray(temp);
      int size=p1740_3->getPadding();
      //! content
      int k=V1740_HDR+V1740_N_MAX_CH;
      data->Clear();
      data->evt_type = V1740_EVENT_TYPE;
      data->b = 3;//for sorter
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
            oj->baselineMean(nBaseline,minVarBaseline);
            double timeData = 0;
            if (mode_selection==0) {
                timeData=oj->ledPos(LED_threshold_LED_740);
            }else if (mode_selection==1) {
                timeData=oj->ledWithCorr(LED_threshold_740);
            }else if (mode_selection==2){
                oj->makeCFD(LED_threshold_740,CFD_delay,CFD_fraction);
                timeData=oj->cfdFast();
            }else if (mode_selection==3){
                oj->makeCFD(LED_threshold_740,CFD_delay,CFD_fraction);
                timeData=oj->cfdSlow();
            }

            //oj->chargeInter(cShort,cLong,timeData,gateOffset,shortGate,longGate);
            chdata->cshort = 0;
            chdata->clong = oj->maxAdcPos(N_MAX_WF_LENGTH)-oj->bL;
            chdata->baseline = oj->bL;
            chdata->finets = timeData;
            delete oj;
        }
        data->AddHit(chdata);
      }
      ProcessEvent(data);
       delete p1740_3;
  }
  return 0;
}



