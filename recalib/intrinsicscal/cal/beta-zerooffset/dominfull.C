//! ROOT::Minuit2::Minuit2Minimizer min;//( ROOT::Minuit2::kMigrad );//works only for  ROOT5

#include <TStyle.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <unistd.h>
#include <map>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

#include "TH1.h"
#include "TH2.h"

#include "TCutG.h"
#include "TLatex.h"
#include "TMath.h"
#include "TCanvas.h"

#include "TGraph.h"
#include "Math/GSLMinimizer.h"
#include "Math/Functor.h"
//#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"
#include "Math/GSLSimAnMinimizer.h"
#include "Math/Functor.h"


Double_t ioffset[16][16];
Double_t islope[16][16];
Double_t ioffseterr[16][16];
Double_t islopeerr[16][16];

Double_t iadcoffsetx[16];
Double_t iadcoffsety[16];



double RosenBrock(const double *xx )
{
  Double_t sumchi2=0.;
  for (int i=0;i<16;i++){
      //if (i<127){
          for (int j=0;j<16;j++){
              sumchi2=sumchi2+(islope[i][j]-xx[32+i]/xx[j])*(islope[i][j]-xx[32+i]/xx[j])/islopeerr[i][j]/islopeerr[i][j]+
                      (ioffset[i][j]-((xx[48+i]-xx[16+j])/xx[j]))*(ioffset[i][j]-((xx[48+i]-xx[16+j])/xx[j]))/ioffseterr[i][j]/ioffseterr[i][j];
          }
      //}else{
      //    for (int j=0;j<128;j++){
      //        sumchi2=sumchi2+(islope[i][j]-1./xx[j])*(islope[i][j]-1./xx[j])/islopeerr[i][j]/islopeerr[i][j]+
      //                (ioffset[i][j]-((1.-xx[128+j])/xx[j]))*(ioffset[i][j]-((1.-xx[128+j])/xx[j]))/ioffseterr[i][j]/ioffseterr[i][j];
      //    }
      //}
  }
  return sumchi2;
}



void dominfull(char* infile, char* outfile, int dssd,double dgain=0.7601200527,double doffset=0.000000015)
{

    //ReadPulserCalibTable("cal_R28_29_may2017.txt",dssd);
    //ReadPulserCalibTable("cal_fake.txt",dssd);
    std::ifstream ifs(infile);
    Int_t nlines=0;
    Double_t offset[256];
    Double_t slope[256];
    Double_t offseterr[256];
    Double_t slopeerr[256];


    for (int i=0;i<256;i++){
        offset[i]=0.;
        slope[i]=1.2;
        offseterr[i]=0.1;
        slopeerr[i]=0.01;
        //assign average value
    }

    Int_t ich;
    while (!ifs.eof()){
        ifs>>ich;
        ifs>>slope[ich]>>slopeerr[ich];
        cout<<ich<<"\t"<<slope[ich]<<"\t"<<slopeerr[ich]<<endl;;
        nlines++;
    }

    for (int i=0;i<16;i++){ //dimension: 128x128
        for (int j=0;j<16;j++){
            int idxx=i*16+j;
            ioffset[i][j]=offset[idxx];
            ioffseterr[i][j]=offseterr[idxx];
            islope[i][j]=slope[idxx];
            islopeerr[i][j]=slopeerr[idxx];
        }
    }

    ROOT::Minuit2::Minuit2Minimizer min;//( ROOT::Minuit2::kMigrad );//works only for  ROOT5
    //ROOT::Math::GSLMinimizer min( ROOT::Math::kVectorBFGS );

    //min.SetMaxFunctionCalls(1000000);
    //min.SetMaxIterations(100000);
    min.SetMaxFunctionCalls(1000000000);
    min.SetMaxIterations(100000000);
    min.SetTolerance(0.001);

    ROOT::Math::Functor f(&RosenBrock,64);
    min.SetFunction(f);
      // Set the free variables to be minimized!
    for (int i=0;i<64;i++){
        if (i<16||(i>=32&&i<48)) min.SetVariable(i,Form("x%d",i),dgain, 0.01);//gain
        else min.SetFixedVariable(i,Form("x%d",i),doffset);//offset
    }
    min.SetFixedVariable(39,"x39",dgain);//gain of X strip 7 as ref
    min.SetFixedVariable(55,"x55",doffset);//offset of X strip 7 (from a pulser run) as ref

    cout<<"Be patient! We are busy feeding our birds..."<<endl;
    //min.SetPrintLevel(1);
    double xxx[64];
    for (int i=0;i<64;i++){
        if (i<16||(i>=32&&i<48)) xxx[i]=1.;
        else xxx[i]=0.;
    }
    cout<<"start chisquare ="<<RosenBrock(xxx)<<endl;
    min.Minimize();
    const double *xs = min.X();
    cout<<"minimized chisquare ="<<RosenBrock(xs)<<endl;
    double fslope[4][32];
    double foffset[4][32];
    //! set default value
    for (int i=0;i<4;i++){
        for (int j=0;j<32;j++){
            fslope[i][j]=1.;
            foffset[i][j]=0.;
        }
    }
    for (int i=0;i<64;i++){
        //cout<<i<<"\t"<<xs[i]<<endl;
        if (i>=32&&i<48) {//x gain
            fslope[dssd][i-32]=xs[i];
        }else if (i<16){//y gain
            fslope[dssd][i+16]=xs[i];
        }else if (i>=48){//x offset
            foffset[dssd][i-48]=xs[i];
        }else{//y offset
            foffset[dssd][i-16+16]=xs[i];
        }
    }
    std::ofstream ofs(outfile,std::ofstream::out | std::ofstream::app);
    for (int i=0;i<4;i++){
        if (i!=dssd) continue;
        for (int j=0;j<32;j++){
            //if (j<128) foffset[i][j]=-fslope[i][j]*iadcoffsetx[j];
            //else foffset[i][j]=-fslope[i][j]*iadcoffsety[j-128];
            cout<<"dssd = "<<i<<"\t"<<j<<"\t"<<foffset[i][j]<<"\t"<<fslope[i][j]<<endl;
            ofs<<i<<"\t"<<j<<"\t"<<foffset[i][j]<<"\t"<<fslope[i][j]<<endl;
        }
    }

}

