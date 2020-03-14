#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <sys/time.h>
#include <signal.h>
#include "TMath.h"
#include "TFile.h"
#include "TROOT.h"
#include "TTree.h"
#include "TCutG.h"
#include "TKey.h"
#include "TStopwatch.h"
#include "TClonesArray.h"
#include "CommandLineInterface.hh"
#include "BELEN.hh"
#include "Clover.hh"
#include "TVectorD.h"
#include "BelenReader.hh"

using namespace TMath;
using namespace std;
bool signal_received = false;
void signalhandler(int sig);
double get_time();


int main(int argc, char* argv[]){
    //! Program start time
    double time_start = get_time();
    TStopwatch timer;
    timer.Start();
    //! Add signal handler
    signal(SIGINT,signalhandler);

    cout << "Belen Reader 2" << endl;
    int Verbose = 0;
    int FillFlag = 1;

    char* InputBELEN = NULL;
    char* OutFile = NULL;
    char* MappingFile = NULL;
    int SplitSize = -1;


    CommandLineInterface* interface = new CommandLineInterface();
    interface->Add("-i", "BELEN input list of files", &InputBELEN);
    interface->Add("-o", "output file", &OutFile);
    interface->Add("-v", "verbose level", &Verbose);
    interface->Add("-map", "mapping file (he3 position)", &MappingFile);
    interface->Add("-f", "fill data or not: 1 fill data 0 no fill (default: fill data)", &FillFlag);
    interface->Add("-s", "Maximum size (interger Mb)", &SplitSize);


    interface->CheckFlags(argc, argv);
    //Complain about missing mandatory arguments
    if(InputBELEN == NULL){
      cout << "No BELEN input list of files given " << endl;
      return 1;
    }
    if(MappingFile == NULL){
      cout << "No Mapping file given, try to use default: He3_mapping.txt" << endl;
      MappingFile = new char[500];
      strcpy(MappingFile,"He3_mapping.txt");
      //return 1;
    }
    if(OutFile == NULL){
      cout << "No output ROOT file given " << endl;
      return 2;
    }
    cout<<"output file: "<<OutFile<< endl;
    TFile* ofile = new TFile(OutFile,"recreate");
    ofile->cd();

    //! Book tree and histograms
    TTree* treeneutron=new TTree("neutron","tree neutron");
    TTree* treegamma=new TTree("gammas","tree gammas");
    TTree* treeanc=new TTree("anc","tree anc");

    BELENHit* neutron=new BELENHit;
    CloverHit* gamma=new CloverHit;
    BELENHit* anc=new BELENHit;
    treeneutron->Branch("neutron",&neutron);
    treegamma->Branch("gamma",&gamma);
    treeanc->Branch("anc",&anc);
    treeneutron->BranchRef();
    treegamma->BranchRef();
    treeanc->BranchRef();

    if (SplitSize>0){
        treeneutron->SetMaxTreeSize(SplitSize*1e6);
        treegamma->SetMaxTreeSize(SplitSize*1e6);
        treeanc->SetMaxTreeSize(SplitSize*1e6);
    }

    //! Read list of files
    string inputfiles[1000];
    ifstream inf(InputBELEN);
    Int_t nfiles;
    inf>>nfiles;

    TVectorD runtime(nfiles+1);
    runtime[0] = 0;
    for (Int_t i=0;i<nfiles;i++){
        runtime[i+1] = 0;
        inf>>inputfiles[i];
        cout<<inputfiles[i]<<endl;
    }

    for (Int_t i=0;i<nfiles;i++){
        BelenReader* blrd = new BelenReader;
        blrd->SetGeoMapping(MappingFile);
        if (FillFlag) blrd->BookTree(treeneutron,treegamma,treeanc,neutron,gamma,anc);
        blrd->Init((char*)inputfiles[i].c_str());
        int ctr=0;
        int total = blrd->GetNEvents();
        long long tstart;
        long long tend;

        double time_last = (double) get_time();
        int start=0;
        double local_time_start = get_time();

        //! event loop
        while(blrd->GetNextEvent()){
            ctr=blrd->GetCurrentEvent();
            if(ctr%1000000 == 0){
                int nevtneutron = blrd->GetCurrentNeutronEvent();
                int nevtgamma = blrd->GetCurrentGammaEvent();
                int nevtanc = blrd->GetCurrentAncEvent();

                double time_end = get_time();
                //! get time

                cout << inputfiles[i] << setw(5) << setiosflags(ios::fixed) << setprecision(1) << (100.*ctr)/total<<" % done  " <<
                  (Float_t)(ctr/1000)/(time_end - local_time_start) << "k events/s " <<
                  (Float_t)(nevtneutron/1000)/(time_end - local_time_start) <<"k neutrons/s  "<<
                   (Float_t)(nevtgamma/1000)/(time_end - local_time_start) <<"k gammas/s "<<
                    (Float_t)(nevtanc/1000)/(time_end - local_time_start) <<"k linhtinhs/s "<<
                   (total-ctr)*(time_end - local_time_start)/(Float_t)ctr << "s to go \r "<<flush;

                time_last = time_end;
            }
            if(signal_received){
                blrd->CloseReader();
              break;
            }
        }

        cout<<"All Hits= "<<blrd->GetCurrentEvent()<<endl;
        cout<<"  Neutron Hits= "<<blrd->GetCurrentNeutronEvent()<<endl;
        cout<<"  Gamma Hits= "<<blrd->GetCurrentGammaEvent()<<endl;
        cout<<"  Anc Hits= "<<blrd->GetCurrentAncEvent()<<endl;
        blrd->CloseReader();
        delete blrd;
        ofile->cd();
    }

    if (FillFlag){
        ofile->cd();
        treeneutron->Write();
        treegamma->Write();
        treeanc->Write();
    }


    ofile->Close();

    cout<<"\n**********************SUMMARY**********************\n"<<endl;
    cout<<"Total run length = "<<runtime[0]<< " seconds"<<endl;
    cout<<"Sub runs length"<<endl;
    for (Int_t i=0;i<nfiles;i++){
        cout<<inputfiles[i]<<" - "<<runtime[i+1]<< " seconds"<<endl;
    }

    //! Finish----------------
    double time_end = get_time();
    cout << "\nProgram Run time: " << time_end - time_start << " s." << endl;
    timer.Stop();
    cout << "CPU time: " << timer.CpuTime() << "\tReal time: " << timer.RealTime() << endl;

}
void signalhandler(int sig){
  if (sig == SIGINT){
    signal_received = true;
  }
}

double get_time(){
    struct timeval t;
    gettimeofday(&t, NULL);
    double d = t.tv_sec + (double) t.tv_usec/1000000;
    return d;

    cout << "AIDA event builder" << endl;


}
