#include "TH1.h"
#include "TH2.h"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TF1.h"
#include "TSpectrum.h"
#include "TMath.h"
#include <iostream>
using namespace std;
Double_t fcn(Double_t *x, Double_t *par) {
    //! par0 mean
    //! par1 sigma1
    //! par2 sigma2
    //! par3 height
    //! par4 const background
    double left=par[3]*exp(-(x[0]-par[0])*(x[0]-par[0])/2/par[1]/par[1]);
    double right=par[3]*exp(-(x[0]-par[0])*(x[0]-par[0])/2/par[2]/par[2]);
    if (x[0]<par[0]) return left+par[4];
    else return right+par[4];
}

void fitsplitnormal(TH1F* hist,Double_t low, Double_t high, Double_t height, Double_t sigma=2.,Double_t bkglvl=0)
{
    TF1* ffinal=new TF1("fcn_splitnormal",fcn,low,high,5);

    ffinal->SetParameter(0,low/2+high/2);
    ffinal->SetParameter(1,sigma);
    ffinal->SetParLimits(1,0,sigma*10);
    ffinal->SetParameter(2,sigma);
    ffinal->SetParLimits(2,0,sigma*10);
    ffinal->SetParameter(3,height);
    ffinal->SetParLimits(3,0,height*10);
    ffinal->SetParameter(4,bkglvl);
    ffinal->SetParLimits(4,0,bkglvl*10);

    hist->Draw();
    hist->Fit("fcn_splitnormal","LQR");
    cout<<"final fit result"<<endl;
    for (Int_t i=0;i<5;i++) {
        cout<<i<<"\t"<<ffinal->GetParameter(i)<<"\t"<<ffinal->GetParError(i)<<endl;
    }
}





