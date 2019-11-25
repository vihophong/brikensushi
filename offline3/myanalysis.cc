/* ************************ DO NOT REMOVE  ****************************** */
#include <iostream>
#include <pmonitor.h>
#include <TSystem.h>
#include <TROOT.h>
#include "pdecoder.h"
#include <dpp.h>
#include "libDataStruct.h"
#include "DataStruct_dict.h"

#define DGTZ_CLK_RES 8
#define V1730_EVENT_TYPE 1
#define V1740_EVENT_TYPE 2
#define TDC_EVENT_TYPE 3

#define V1730_N_MAX_BOARD 5
#define V1730_N_MAX_CH 16
#define N_MAX_INFO 9

#define TDC_BOARD_N 50
#define TDC_PACKET 10
#define TDC_N_CH 64



#define V1740_HDR 6
#define V1740_N_CH 64
#define V1740_N_BOARD 4
#define V1740_N_BOARD_G1 2
#define V1740_N_BOARD_G2 2
#define V1740_PACKET_0 100
#define V1740_PACKET_1 101
#define V1740_PACKET_2 102
#define V1740_PACKET_3 103
#define V1740_N_MAX_CH 64
#define NSBL 8
// 5000000 - 15 % memory = 1.2 Gb

/* ********************************************************************** */

#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <vector>
#include <map>
#include <TCanvas.h>

#define N_REFESH_WF 1000
#define N_MAX_WF_LENGTH 500

#define COUNTS_TO_DISPLAY 10000
#define COUNTS_TO_CAL 100
#define MAX_N_SAMPLE 300

#define MAX_MAP_LENGTH_G1 50
#define MAX_MAP_LENGTH_G2 50
#define TS_WIN_LOW_G1 100
#define TS_WIN_HIGH_G1 100
#define TS_WIN_LOW_G2 100
#define TS_WIN_HIGH_G2 100


using namespace std;

//! trees
TFile* fout;
//!group 1 data container
TTree* tree;
Int_t evt[V1740_N_BOARD_G1];
ULong64_t dgtz_ts[V1740_N_BOARD_G1];
UShort_t dgtz_nsample[V1740_N_CH*V1740_N_BOARD_G1];
Double_t dgtz_e[V1740_N_CH*V1740_N_BOARD_G1];
Double_t dgtz_bl[V1740_N_CH*V1740_N_BOARD_G1];
Int_t dgtz_ch[V1740_N_CH*V1740_N_BOARD_G1];
UShort_t dgtz_sample[V1740_N_CH*V1740_N_BOARD_G1][MAX_N_SAMPLE];
UShort_t dgtz_waveform[V1740_N_CH*V1740_N_BOARD_G1][MAX_N_SAMPLE];


//!group 2 data container
TTree* tree_g2;
Int_t evt_g2[V1740_N_BOARD_G2];
ULong64_t dgtz_ts_g2[V1740_N_BOARD_G2];
UShort_t dgtz_nsample_g2[V1740_N_CH*V1740_N_BOARD_G2];
Double_t dgtz_e_g2[V1740_N_CH*V1740_N_BOARD_G2];
Double_t dgtz_bl_g2[V1740_N_CH*V1740_N_BOARD_G2];
Int_t dgtz_ch_g2[V1740_N_CH*V1740_N_BOARD_G2];
UShort_t dgtz_sample_g2[V1740_N_CH*V1740_N_BOARD_G2][MAX_N_SAMPLE];
UShort_t dgtz_waveform_g2[V1740_N_CH*V1740_N_BOARD_G2][MAX_N_SAMPLE];


Int_t ntotalg1;
Int_t ntriggerg1[V1740_N_BOARD_G1];
Int_t ntriggerg2[V1740_N_BOARD_G2];

Int_t ntotalg2;

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



void Init(){

}



int phsave(const char* filename)
{
        ntotalg1=0;
        ntotalg2=0;
        fout = new TFile(filename,"RECREATE");
        tree = new TTree("group1","group1");
        tree->Branch("evt",&evt,Form("evt[%d]/I",V1740_N_CH*V1740_N_BOARD_G1));
        tree->Branch("dgtz_e",dgtz_e,Form("dgtz_e[%d]/D",V1740_N_CH*V1740_N_BOARD_G1));
        tree->Branch("dgtz_bl",dgtz_bl,Form("dgtz_bl[%d]/D",V1740_N_CH*V1740_N_BOARD_G1));
        tree->Branch("dgtz_ch",dgtz_ch,Form("dgtz_ch[%d]/I",V1740_N_CH*V1740_N_BOARD_G1));
        tree->Branch("dgtz_nsample",dgtz_nsample,Form("dgtz_nsample[%d]/s",V1740_N_CH*V1740_N_BOARD_G1));
        tree->Branch("dgtz_ts",dgtz_ts,Form("dgtz_ts[%d]/l",V1740_N_CH*V1740_N_BOARD_G1));
        tree->Branch("dgtz_waveform",dgtz_waveform,Form("dgtz_waveform[%d][%d]/s",V1740_N_CH*V1740_N_BOARD_G1,MAX_N_SAMPLE));
        tree->Branch("dgtz_sample",dgtz_sample,Form("dgtz_sample[%d][%d]/s",V1740_N_CH*V1740_N_BOARD_G1,MAX_N_SAMPLE));

        tree_g2 = new TTree("group2","group2");
        tree_g2->Branch("evt",&evt_g2,Form("evt[%d]/I",V1740_N_CH*V1740_N_BOARD_G2));
        tree_g2->Branch("dgtz_e",dgtz_e_g2,Form("dgtz_e[%d]/D",V1740_N_CH*V1740_N_BOARD_G2));
        tree_g2->Branch("dgtz_bl",dgtz_bl_g2,Form("dgtz_bl[%d]/D",V1740_N_CH*V1740_N_BOARD_G2));
        tree_g2->Branch("dgtz_ch",dgtz_ch_g2,Form("dgtz_ch[%d]/I",V1740_N_CH*V1740_N_BOARD_G2));
        tree_g2->Branch("dgtz_nsample",dgtz_nsample_g2,Form("dgtz_nsample[%d]/s",V1740_N_CH*V1740_N_BOARD_G2));
        tree_g2->Branch("dgtz_ts",dgtz_ts_g2,Form("dgtz_ts[%d]/l",V1740_N_CH*V1740_N_BOARD_G2));
        tree_g2->Branch("dgtz_waveform",dgtz_waveform_g2,Form("dgtz_waveform[%d][%d]/s",V1740_N_CH*V1740_N_BOARD_G2,MAX_N_SAMPLE));
        tree_g2->Branch("dgtz_sample",dgtz_sample_g2,Form("dgtz_sample[%d][%d]/s",V1740_N_CH*V1740_N_BOARD_G2,MAX_N_SAMPLE));

	for (int i=0;i<V1740_N_BOARD_G1;i++) ntriggerg1[i]=0;
	for (int i=0;i<V1740_N_BOARD_G2;i++) ntriggerg2[i]=0;
        return 1;
}
void CloseEventG1(){    
    tree->Fill();
    //! reset default value
    for (int i=0;i<V1740_N_BOARD_G1;i++){
        dgtz_ts[i]=0;
        evt[i]=0;
    }
    for (int i=0;i<V1740_N_CH*V1740_N_BOARD_G1;i++){
        dgtz_nsample[i]=-9999;
        dgtz_e[i]=-9999;
        dgtz_bl[i]=-9999;
        dgtz_ch[i]=-9999;
        for (int j=0;j<MAX_N_SAMPLE;j++){
            dgtz_sample[i][j]=-9999;
            dgtz_waveform[i][j]=-9999;
        }
    }
}
void CloseEventG2(){
    tree_g2->Fill();
    //! reset default value
    for (int i=0;i<V1740_N_BOARD_G2;i++){
        dgtz_ts_g2[i]=0;
        evt_g2[i]=0;
    }
    for (int i=0;i<V1740_N_CH*V1740_N_BOARD_G2;i++){
        dgtz_nsample_g2[i]=-9999;
        dgtz_e_g2[i]=-9999;
        dgtz_bl_g2[i]=-9999;
        dgtz_ch_g2[i]=-9999;
        for (int j=0;j<MAX_N_SAMPLE;j++){
            dgtz_sample_g2[i][j]=-9999;
            dgtz_waveform_g2[i][j]=-9999;
        }
    }
}

void CopyG1toG2(){
    for (int i=0;i<V1740_N_BOARD_G2;i++){
        dgtz_ts_g2[i]=dgtz_ts[i];
        evt_g2[i]=evt[i];
    }
    for (int i=0;i<V1740_N_CH*V1740_N_BOARD_G2;i++){
        dgtz_nsample_g2[i]=dgtz_nsample[i];
        dgtz_e_g2[i]=dgtz_e[i];
        dgtz_bl_g2[i]=dgtz_bl[i];
        dgtz_ch_g2[i]=dgtz_ch[i];
        for (int j=0;j<MAX_N_SAMPLE;j++){
            dgtz_sample_g2[i][j]=dgtz_sample[i][j];
            dgtz_waveform_g2[i][j]=dgtz_waveform[i][j];
        }
    }
}

int DoSort(bool flagend,NIGIRI*data=0){
    //! digitizer group 1
    if (datamap_dgtz0.size()>MAX_MAP_LENGTH_G1||flagend){
      for(it_datamap_dgtz0=datamap_dgtz0.begin();it_datamap_dgtz0!=datamap_dgtz0.end();it_datamap_dgtz0++){
        Long64_t ts=(Long64_t)it_datamap_dgtz0->first;
        NIGIRI* hit0=(NIGIRI*)it_datamap_dgtz0->second;
        Long64_t ts1 = ts - TS_WIN_LOW_G1;
        Long64_t ts2 = ts + TS_WIN_HIGH_G1;
        Long64_t corrts = 0;
        it_datamap_dgtz1 = datamap_dgtz1.lower_bound(ts1);
        while(it_datamap_dgtz1!=datamap_dgtz1.end()&&it_datamap_dgtz1->first<ts2){
            corrts = (Long64_t) it_datamap_dgtz1->first;
            NIGIRI* hit1=(NIGIRI*)it_datamap_dgtz1->second;
            //! fill data for dgtz1 here
            dgtz_ts[1]=hit1->ts;
            evt[1]=hit1->evt;
            for (unsigned int i=0;i<hit1->GetMult();i++){
                NIGIRIHit* hittemp=hit1->GetHit(i);
                int ch=hittemp->ch+V1740_N_CH;
                dgtz_nsample[ch]=hittemp->nsample;
                dgtz_e[ch]=hittemp->clong;
                dgtz_bl[ch]=hittemp->baseline;
                dgtz_ch[ch]=ch;
                int intcnt=0;
                for (std::vector<UShort_t>::iterator it = hittemp->pulse.begin() ; it != hittemp->pulse.end(); ++it){
                    if(intcnt<N_MAX_WF_LENGTH){
                        dgtz_waveform[ch][intcnt]=*it;
                        dgtz_sample[ch][intcnt]=intcnt;
                    }
                    intcnt++;
                }
            }
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
        //! fill data for dgtz0 here
        dgtz_ts[0]=hit0->ts;
        evt[0]=hit0->evt;
        for (unsigned int i=0;i<hit0->GetMult();i++){
            NIGIRIHit* hittemp=hit0->GetHit(i);
            int ch=hittemp->ch;
            dgtz_nsample[ch]=hittemp->nsample;
            dgtz_e[ch]=hittemp->clong;
            dgtz_bl[ch]=hittemp->baseline;
            dgtz_ch[ch]=ch;
            int intcnt=0;
            for (std::vector<UShort_t>::iterator it = hittemp->pulse.begin() ; it != hittemp->pulse.end(); ++it){
                if(intcnt<N_MAX_WF_LENGTH){
                    dgtz_waveform[ch][intcnt]=*it;
                    dgtz_sample[ch][intcnt]=intcnt;
                }
                intcnt++;
            }
        }

        //!Generate some fake data
        /*
        if (ntotalg1<1000) {
            CopyG1toG2();
            CloseEventG2();
            ntotalg2++;
        }
        */

        CloseEventG1();




        ntotalg1++;
        if (ntotalg1%1000==0&&ntotalg1>0) cout<<ntotalg1<<" events in group one + "<<ntotalg2<<" events in group two are processed \r"<<flush;
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
        while(it_datamap_dgtz3!=datamap_dgtz3.end()&&it_datamap_dgtz3->first<ts2){
            corrts = (Long64_t) it_datamap_dgtz3->first;
            NIGIRI* hit1=(NIGIRI*)it_datamap_dgtz3->second;
            //! fill data for dgtz3 here
            dgtz_ts_g2[1]=hit1->ts;
            evt_g2[1]=hit1->evt;
            for (unsigned int i=0;i<hit1->GetMult();i++){
                NIGIRIHit* hittemp=hit1->GetHit(i);
                int ch=hittemp->ch+V1740_N_CH;
                dgtz_nsample_g2[ch]=hittemp->nsample;
                dgtz_e_g2[ch]=hittemp->clong;
                dgtz_bl_g2[ch]=hittemp->baseline;
                dgtz_ch_g2[ch]=ch;
                int intcnt=0;
                for (std::vector<UShort_t>::iterator it = hittemp->pulse.begin() ; it != hittemp->pulse.end(); ++it){
                    if(intcnt<N_MAX_WF_LENGTH){
                        dgtz_waveform_g2[ch][intcnt]=*it;
                        dgtz_sample_g2[ch][intcnt]=intcnt;
                    }
                    intcnt++;
                }
            }
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
        //! fill data for dgtz2 here
        dgtz_ts_g2[0]=hit0->ts;
        evt_g2[0]=hit0->evt;
        for (unsigned int i=0;i<hit0->GetMult();i++){
            NIGIRIHit* hittemp=hit0->GetHit(i);
            int ch=hittemp->ch;
            dgtz_nsample_g2[ch]=hittemp->nsample;
            dgtz_e_g2[ch]=hittemp->clong;
            dgtz_bl_g2[ch]=hittemp->baseline;
            dgtz_ch_g2[ch]=ch;
            int intcnt=0;
            for (std::vector<UShort_t>::iterator it = hittemp->pulse.begin() ; it != hittemp->pulse.end(); ++it){
                if(intcnt<N_MAX_WF_LENGTH){
                    dgtz_waveform_g2[ch][intcnt]=*it;
                    dgtz_sample_g2[ch][intcnt]=intcnt;
                }
                intcnt++;
            }
        }
        CloseEventG2();
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

void ProcessEvent(NIGIRI* data){
    if (data->evt_type==V1740_EVENT_TYPE){
        DoSort(false,data);
    }
}

void CloseMe(){
    DoSort(true);
    TTree* otree1= (TTree*)fout->Get("group1");
    otree1->Write();
    TTree* otree2= (TTree*)fout->Get("group2");
    otree2->Write();
    fout->Close();
    cout<<"b1 total number of triggers = "<<ntriggerg1[0]<<endl;
    cout<<"b2 total number of triggers = "<<ntriggerg1[1]<<endl;
    cout<<"total number of triggers in group 1 = "<<ntotalg1<<endl;

    cout<<"\nb3 total number of triggers"<<ntriggerg2[0]<<endl;
    cout<<"b4 total number of triggers"<<ntriggerg2[1]<<endl;
    cout<<"total number of triggers in group 2 = "<<ntotalg2<<endl;    
}




