#include "TSpectrum.h"
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TFile.h"

#include "TH1.h"
#include "TH2.h"
#include "TH1F.h"

#include "TCut.h"
#include "TCutG.h"
#include "TLatex.h"
#include "TSpectrum.h"
#include "TF1.h"
#include <algorithm>

#include <TROOT.h>


void calibspec(TH1F *h,Double_t rangemin1=150.,Double_t rangeplus1=100,Double_t rangemin2=30,Double_t rangeplus2=50) {
   gROOT->Reset();
   gROOT->Clear();

//   Double_t defaultheight=500;
//   Double_t defaultsigma=5;
//   Double_t max=1000;
//   Double_t min=200;

   Double_t defaultheight=500;
   Double_t defaultsigma=17;
   Double_t max=2000;
   Double_t min=500;


   Double_t source[10000];
   TCanvas* c1;
   c1=new TCanvas("c1","c1",1200,900);
   c1->Divide(2,2);
   c1->cd(1);   
   gStyle->SetOptStat(11111111);
   h->Draw();
   TH1F* d1=(TH1F*)h->Clone();
   d1->Reset();
   h->Draw();
   TSpectrum *s = new TSpectrum();
   Int_t k=0;
   Int_t minbin,maxbin;



   for (Int_t i = 0; i < h->GetNbinsX(); i++) {
       if (h->GetXaxis()->GetBinCenter(i+1)>min&&h->GetXaxis()->GetBinCenter(i+1)<max){
           source[k]=h->GetBinContent(i + 1);
           if (k==0) minbin=i+1;
           maxbin=k+1;
           k++;
       }
   }

   s->Background(source,h->GetNbinsX(),15,TSpectrum::kBackDecreasingWindow,TSpectrum::kBackOrder8,kTRUE,
   TSpectrum::kBackSmoothing5,kTRUE);

   for (Int_t i = 0; i < h->GetNbinsX(); i++) {
       if (i>=minbin&&i<=maxbin)
       d1->SetBinContent(i + 1,source[i-minbin]);
       else
      d1->SetBinContent(i + 1,h->GetBinContent(i+1));

   }

   c1->cd(2);
   d1->SetLineColor(kRed);
   d1->Draw("SAME L");
   TH1F* d2=(TH1F*)h->Clone();
   d2->Add(d1,-1);
   d2->SetLineColor(kGreen);
   d2->Draw("same");

   s->Search(d2);
   Double_t* peaks=s->GetPositionX();

   //! sorting
   std::vector <Double_t> mappeak;
   std::vector <Double_t>::iterator mappeak_it;
   for (Int_t i=0;i<s->GetNPeaks();i++){
       mappeak.push_back(peaks[i]);
       cout<<peaks[i]<<endl;
   }
   std::sort (mappeak.begin(), mappeak.begin()+4);
   k=0;
   for (mappeak_it=mappeak.begin();mappeak_it!=mappeak.end();mappeak_it++){
       peaks[k]=mappeak[k];
       cout<<peaks[k]<<endl;
       k++;
   }



   Double_t peaklib[]={482,555,976,1049};


   TString fdef("gaus(0)+gaus(3)");
   TF1* ffinal=new TF1("ffinal",(char*)fdef.Data(),peaks[0]-rangemin1,peaks[1]+rangeplus1);

   ffinal->SetParameter(0,defaultheight);// peak heigh
   ffinal->SetParameter(1,peaks[0]);// centroid
   ffinal->SetParameter(2,defaultsigma);// sigma
   ffinal->SetParameter(3,defaultheight);// peak heigh
   ffinal->SetParameter(4,peaks[1]);// centroid
   ffinal->SetParameter(5,defaultsigma);// sigma

   d2->Fit("ffinal","LQR","");
   for (Int_t i=0;i<6;i++) {
       cout<<i<<"\t"<<ffinal->GetParameter(i)<<endl;
   }

   TF1* ffinal2=new TF1("ffinal2",(char*)fdef.Data(),peaks[2]-rangemin2,peaks[3]+rangeplus2);

   ffinal2->SetParameter(0,defaultheight);// peak heigh
   ffinal2->SetParameter(1,peaks[2]);// centroid
   ffinal2->SetParameter(2,defaultsigma);// sigma
   ffinal2->SetParameter(3,defaultheight);// peak heigh
   ffinal2->SetParameter(4,peaks[3]);// centroid
   ffinal2->SetParameter(5,defaultsigma);// sigma
   d2->Fit("ffinal2","LQR","");
   for (Int_t i=0;i<6;i++) {
       cout<<i<<"\t"<<ffinal->GetParameter(i)<<endl;
   }
   ffinal->Draw("same");

   //replace peaks with fit parameter;

   Double_t peaks_wfit[4];

   peaks_wfit[0]=ffinal->GetParameter(1);
   peaks_wfit[1]=ffinal->GetParameter(4);

   peaks_wfit[2]=ffinal2->GetParameter(1);
   peaks_wfit[3]=ffinal2->GetParameter(4);
   cout<<"peak wit fit = "<<peaks_wfit[0]<<"\t"<<peaks_wfit[1]<<"\t"<<peaks_wfit[2]<<"\t"<<peaks_wfit[3]<<"\t"<<endl;


   TLatex latex;
   latex.SetTextSize(0.025);
   latex.SetTextAlign(13);  //align at top


   c1->cd(3);
   gStyle->SetOptStat(11111111);
   TGraph *grcal=new TGraph(4,peaks,peaklib);
   grcal->SetNameTitle("grpeakfind","TSpectrum peak finding");
    grcal->Draw("AP*");    
    grcal->Fit("pol1");
    Double_t offset=grcal->GetFunction("pol1")->GetParameter(0);
    Double_t gain=grcal->GetFunction("pol1")->GetParameter(1);
    Double_t x2ndf=grcal->GetFunction("pol1")->GetChisquare()/grcal->GetFunction("pol1")->GetNDF();
    latex.DrawLatex(peaks[1],peaklib[1],Form("#chi^{2}/NDF=%.2f...gain=%.2f...offset=%.2f",x2ndf,gain,offset));

    c1->cd(4);
    gStyle->SetOptStat(11111111);
    TGraph *grcal_wfit=new TGraph(4,peaks_wfit,peaklib);
    grcal_wfit->SetNameTitle("grpeakfit","Peak fitting");
    grcal_wfit->Draw("AP*");
    grcal_wfit->Fit("pol1");

     offset=grcal_wfit->GetFunction("pol1")->GetParameter(0);
     gain=grcal_wfit->GetFunction("pol1")->GetParameter(1);
     x2ndf=grcal_wfit->GetFunction("pol1")->GetChisquare()/grcal_wfit->GetFunction("pol1")->GetNDF();
     latex.DrawLatex(peaks_wfit[1],peaklib[1],Form("#chi^{2}/NDF=%.2f...gain=%.2f...offset=%.2f",x2ndf,gain,offset));
    cout<<"resolution"<<endl;
    cout<<"1\t"<<ffinal->GetParameter(2)*gain<<endl;
    cout<<"2\t"<<ffinal->GetParameter(5)*gain<<endl;
    cout<<"3\t"<<ffinal2->GetParameter(2)*gain<<endl;
    cout<<"4\t"<<ffinal2->GetParameter(5)*gain<<endl;

    std::ofstream ofs("result_calibtable.txt",std::ofstream::out | std::ofstream::app);
    ofs<<h->GetName()<<"\t"<<offset<<"\t"<<gain<<"\t"<<x2ndf<<endl;
}
void fitcalibcurve4(Int_t ch,Double_t peak1=482,Double_t peak2=555,Double_t peak3=976,Double_t peak4=1049) {
    TLatex latex;
    latex.SetTextSize(0.025);
    latex.SetTextAlign(13);  //align at top
    Double_t peaks[]={peak1,peak2,peak3,peak4};
    Double_t peaklib[]={482,555,976,1049};
    gStyle->SetOptStat(11111111);
    TGraph *grcal=new TGraph(4,peaks,peaklib);
    grcal->SetNameTitle("grpeakfind","TSpectrum peak finding");
     grcal->Draw("AP*");
     grcal->Fit("pol1");
     Double_t offset=grcal->GetFunction("pol1")->GetParameter(0);
     Double_t gain=grcal->GetFunction("pol1")->GetParameter(1);
     Double_t x2ndf=grcal->GetFunction("pol1")->GetChisquare()/grcal->GetFunction("pol1")->GetNDF();
     latex.DrawLatex(peaks[1],peaklib[1],Form("#chi^{2}/NDF=%.2f...gain=%.2f...offset=%.2f",x2ndf,gain,offset));
     std::ofstream ofs("result_calibtable.txt",std::ofstream::out | std::ofstream::app);
     ofs<<Form("h%d",ch)<<"\t"<<offset<<"\t"<<gain<<"\t"<<x2ndf<<endl;
}

void fitcalibcurve3(Int_t ch,Double_t peak1=482,Double_t peak2=555,Double_t peak3=976) {
    TLatex latex;
    latex.SetTextSize(0.025);
    latex.SetTextAlign(13);  //align at top
    Double_t peaks[]={peak1,peak2,peak3};
    Double_t peaklib[]={482,976,1049};
    gStyle->SetOptStat(11111111);
    TGraph *grcal=new TGraph(3,peaks,peaklib);
    grcal->SetNameTitle("grpeakfind","TSpectrum peak finding");
     grcal->Draw("AP*");
     grcal->Fit("pol1");
     Double_t offset=grcal->GetFunction("pol1")->GetParameter(0);
     Double_t gain=grcal->GetFunction("pol1")->GetParameter(1);
     Double_t x2ndf=grcal->GetFunction("pol1")->GetChisquare()/grcal->GetFunction("pol1")->GetNDF();
     latex.DrawLatex(peaks[1],peaklib[1],Form("#chi^{2}/NDF=%.2f...gain=%.2f...offset=%.2f",x2ndf,gain,offset));
     std::ofstream ofs("result_calibtable.txt",std::ofstream::out | std::ofstream::app);
     ofs<<Form("h%d",ch)<<"\t"<<offset<<"\t"<<gain<<"\t"<<x2ndf<<endl;
}


void calibspecPrint(TH1F *h,Double_t rangemin1=150.,Double_t rangeplus1=100,Double_t rangemin2=30,Double_t rangeplus2=50) {
   gROOT->Reset();
   gROOT->Clear();

   //   Double_t defaultheight=500;
   //   Double_t defaultsigma=5;
   //   Double_t max=1000;
   //   Double_t min=200;

      Double_t defaultheight=500;
      Double_t defaultsigma=17;
      Double_t max=2000;
      Double_t min=500;


   Double_t source[10000];
   TCanvas* c1;
   c1=new TCanvas("c1","c1",1200,900);
   c1->Divide(2,2);
   c1->cd(1);
   gStyle->SetOptStat(11111111);
   h->Draw();
   TH1F* d1=(TH1F*)h->Clone();
   d1->Reset();
   h->Draw();
   TSpectrum *s = new TSpectrum();
   Int_t k=0;
   Int_t minbin,maxbin;



   for (Int_t i = 0; i < h->GetNbinsX(); i++) {
       if (h->GetXaxis()->GetBinCenter(i+1)>min&&h->GetXaxis()->GetBinCenter(i+1)<max){
           source[k]=h->GetBinContent(i + 1);
           if (k==0) minbin=i+1;
           maxbin=k+1;
           k++;
       }
   }

   s->Background(source,h->GetNbinsX(),15,TSpectrum::kBackDecreasingWindow,TSpectrum::kBackOrder8,kTRUE,
   TSpectrum::kBackSmoothing5,kTRUE);

   for (Int_t i = 0; i < h->GetNbinsX(); i++) {
       if (i>=minbin&&i<=maxbin)
       d1->SetBinContent(i + 1,source[i-minbin]);
       else
      d1->SetBinContent(i + 1,h->GetBinContent(i+1));

   }

   c1->cd(2);
   d1->SetLineColor(kRed);
   d1->Draw("SAME L");
   TH1F* d2=(TH1F*)h->Clone();
   d2->Add(d1,-1);
   d2->SetLineColor(kGreen);
   d2->Draw("same");

   s->Search(d2);
   Double_t* peaks=s->GetPositionX();

   //! sorting
   std::vector <Double_t> mappeak;
   std::vector <Double_t>::iterator mappeak_it;
   for (Int_t i=0;i<s->GetNPeaks();i++){
       mappeak.push_back(peaks[i]);
       cout<<peaks[i]<<endl;
   }
}


