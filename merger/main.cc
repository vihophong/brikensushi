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
#include "TVectorD.h"

#include "Merger.hh"

using namespace TMath;
using namespace std;
bool signal_received = false;
void signalhandler(int sig);
double get_time();

int main(int argc, char* argv[]){
    cout<<"BRIKEN Sushi Merger"<<endl;

    char* InputWASABI = NULL;
    char* InputBELEN = NULL;
    char* InputBIGRIPS = NULL;
    char* InputPID = NULL;
    char* OutFile = NULL;

    CommandLineInterface* interface = new CommandLineInterface();
    interface->Add("-w", "Wasabi input file", &InputWASABI);
    interface->Add("-bl", "BELEN input file", &InputBELEN);
    interface->Add("-br", "Bigrips input file", &InputBIGRIPS);
    interface->Add("-pid", "Input PID file", &InputPID);
    interface->Add("-o", "output file", &OutFile);

    interface->CheckFlags(argc, argv);

    //Complain about missing mandatory arguments
    if(InputBELEN == NULL){
      cout << "No BELEN input list of files given " << endl;
      //return 1;
    }
    if(InputWASABI == NULL){
      cout << "No AIDA input list of files given " << endl;
      //return 1;
    }
    if(InputBIGRIPS == NULL){
      cout << "No Bigrips input list of files given " << endl;
      //return 1;
    }
    if(OutFile == NULL){
      cout << "No output ROOT file given " << endl;
      //return 2;
    }

    //! IonBeta corr (test)
    if (!(InputWASABI==NULL)&&InputBIGRIPS==NULL&&InputBELEN==NULL){
        TFile* ofile = new TFile(OutFile,"recreate");
        ofile->cd();
        Merger* merge=new Merger;
        merge->SetWasabiFile(InputWASABI);
        merge->Init();
        merge->ReadWasabi();

        ofile->cd();
        merge->BookIonBetaTree();

        merge->MergeIonBeta();
        ofile->cd();
        merge->GetTreeRI(-1)->Write();
        merge->getProbeHisto1D(0)->Write();
        ofile->Close();
    }

    //! execution implant corr
    if (!(InputWASABI==NULL)&&!(InputBIGRIPS==NULL)&&!(InputBELEN==NULL)){
        TFile* ofile = new TFile(OutFile,"recreate");
        ofile->cd();
        Merger* merge=new Merger;
        merge->SetWasabiFile(InputWASABI);
        merge->SetBigripsFile(InputBIGRIPS);
        merge->SetBrikenFile(InputBELEN);
        merge->Init();
        merge->ReadBigrips();
        merge->ReadBRIKEN();
        merge->ReadWasabi(0);

        merge->DoAddback();

        ofile->cd();
        if (InputPID!=NULL) merge->ReadPID(InputPID);
        merge->BookImplantTree();
        merge->MergeImplant();
        ofile->cd();

        for (Int_t i=0;i<merge->GetNri();i++){
            merge->GetTreeImpRI(i)->Write();
        }

        merge->GetTreeImpRI(-1)->Write();
        merge->getProbeHisto1D(0)->Write();
        ofile->Close();
    }

    return 0;
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
}

