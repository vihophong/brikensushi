#include "TH1.h"
#include "TH2.h"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TF1.h"
#include "TSpectrum.h"
#include <iostream>
using namespace std;

Double_t fitpeaksingle(TH1F* hist,Double_t low, Double_t high,Double_t bkglvl,Double_t sigma)
{
    TString fdef("pol1(0)+gaus(2)");
    TF1* ffinal=new TF1("ffinal",(char*)fdef.Data(),low,high);
    //TF1* fgaus=new TF1("fgaus","gaus",low,high);

    ffinal->SetParameter(0,bkglvl);
    //ffinal->SetParameter(1,0.);
    ffinal->FixParameter(1,0.);



    //hist->SetAxisRange(low,high);
    TSpectrum* sp=new TSpectrum();
    sp->Search(hist);
    Double_t * xpeaks=sp->GetPositionX();

    Double_t peakcenter=xpeaks[0];
    Double_t peakheight=hist->GetBinContent(hist->GetXaxis()->FindBin(xpeaks[0]))-bkglvl;

    //cout<<xpeaks[0]<<endl;

    ffinal->SetParameter(2,peakheight);// peak heigh
    ffinal->SetParameter(3,peakcenter);// centroid
    ffinal->SetParameter(4,sigma);// sigma

    ffinal->SetRange(peakcenter-0.5*sigma,high);
    hist->Fit("ffinal","LQR");
//    cout<<"final fit result"<<endl;
//    for (Int_t i=0;i<5;i++) {
//        cout<<i<<"\t"<<ffinal->GetParameter(i)<<"\t"<<ffinal->GetParError(i)<<endl;
//    }
    //Double_t peakarea=ffinal->GetParameter(2)*ffinal->GetParameter(4)*sqrt(2*TMath::Pi());
    //Double_t peakareaerr=peakarea*sqrt(pow((ffinal->GetParError(2)/ffinal->GetParameter(2)),2.) + pow((ffinal->GetParError(4)/ffinal->GetParameter(4)),2));
    Double_t peaksigma=ffinal->GetParameter(4);
    Double_t peaksigmaerr=ffinal->GetParError(4);
    Double_t peakcentroid=ffinal->GetParameter(3);
    Double_t peakcentroiderr=ffinal->GetParError(3);

    //fgaus->SetParameter(0,ffinal->GetParameter(2));
    //fgaus->SetParameter(1,ffinal->GetParameter(3));
    //fgaus->SetParameter(2,ffinal->GetParameter(4));
    //Double_t peakarea2=fgaus->Integral(low,high);
    //cout<<"peak area="<<peakarea<<"\t error="<<peakareaerr<<" rel. err="<<peakareaerr/peakarea<<endl;

    cout<<hist->GetName()<<"\t"<<peakcentroid<<"\t"<<peakcentroiderr<<"\t"<<peaksigma<<"\t"<<peaksigmaerr<<endl;
    return peakcentroid+peaksigma*3;
}
void fitpeakgetthreshold(char* infile, char* outfile, Double_t low, Double_t high,Double_t bkglvl=0.,Double_t sigma=17.)
{
    TFile* file1=TFile::Open(infile);
    TH2F* h2=(TH2F*) file1->Get("h2");

    //! output
    TFile* fout = new TFile(outfile,"RECREATE");

    TH1F* hproj[128];
    TH1F* hprojthr[128];
    for (Int_t i=0;i<128;i++){
        hproj[i]=(TH1F*) h2->ProjectionY(Form("ch%d",i),i+1,i+1);
        hprojthr[i]=(TH1F*) hproj[i]->Clone();
        hprojthr[i]->SetName(Form("bellow_thr_%s",hproj[i]->GetName()));
        Double_t thr=fitpeaksingle(hproj[i],low,high,bkglvl,sigma);

        for (Int_t j=0;j<hprojthr[i]->GetNbinsX();j++) {
            if (hprojthr[i]->GetBinCenter(j+1)>thr) hprojthr[i]->SetBinContent(j,1000.);
        }

        hproj[i]->Write();
        hprojthr[i]->Write();
    }

    fout->Close();


}

