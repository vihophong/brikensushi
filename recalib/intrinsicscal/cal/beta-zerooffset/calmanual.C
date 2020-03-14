#include <TH2.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TROOT.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <TString.h>
#include <string>
#include <TSpectrum.h>
#include <TText.h>
#include <TF1.h>
#include "TExec.h"
#include <vector>
#include "TSystem.h"
#include "TGraph.h"
#include "TAxis.h"

#include "gsl/gsl_fit.h"

void calmanual(Int_t dssd=0)
{
    TFile *f=TFile::Open(Form("dssd%d_beta.root",dssd));
    TH2F* h2g[256];
    TCanvas* c1=new TCanvas("hisname","hisname",1200,900);
    for (Int_t i=0;i<256;i++){
        h2g[i]=(TH2F*) f->Get(Form("h%d",i));
        gStyle->SetOptStat(0);
        h2g[i]->Draw("colz");

        c1->Update();
        c1->SetName(Form("his%d",i));
        TExec ex("ex1",".x exec1.C");
        ex.Draw();
        c1->WaitPrimitive();
        gSystem->Sleep(500);
    }
    c1->WaitPrimitive();
}
void checkfit(Int_t dssd=0)
{
    std::ifstream ifslowpoints(Form("dssd%d_beta_lowpoints.txt",dssd));
    Double_t xlow[256];
    Double_t ylow[256];
    std::string tempstr;
    for (Int_t i=0;i<256;i++){
        ifslowpoints>>tempstr>>xlow[i]>>ylow[i];
    }
    ifslowpoints.close();
    std::ifstream ifshighpoints(Form("dssd%d_beta_highpoints.txt",dssd));
    Double_t xhigh[256];
    Double_t yhigh[256];
    for (Int_t i=0;i<256;i++){
        ifshighpoints>>tempstr>>xhigh[i]>>yhigh[i];
    }
    ifshighpoints.close();

    TFile *f=TFile::Open(Form("dssd%d_beta.root",dssd));
    TH2F* h2g[256];
    TF1* f2g[256];
    TGraph* g2g[256];

    TCanvas* c1=new TCanvas("hisname","hisname",800,600);

    std::ofstream ofs(Form("dssd%d_beta_coefs.txt",dssd));
    for (Int_t i=0;i<256;i++){
        h2g[i]=(TH2F*) f->Get(Form("h%d",i));
        h2g[i]->Draw("colz");
        //cout<<xlow[i]<<"\t"<<ylow[i]<<"\t"<<xhigh[i]<<"\t"<<yhigh[i]<<endl;
        Double_t x[]={xlow[i],xhigh[i]};
        Double_t y[]={ylow[i],yhigh[i]};
        g2g[i]=new TGraph(2,x,y);
        Double_t a0,a1;
        Int_t ifail;
        g2g[i]->LeastSquareLinearFit(2,a0,a1,ifail);

        f2g[i]=new TF1(Form("f%d",i),"pol1",h2g[i]->GetXaxis()->GetXmin(),h2g[i]->GetXaxis()->GetXmax());
        f2g[i]->FixParameter(0,a0);
        f2g[i]->FixParameter(1,a1);
        g2g[i]->SetMarkerStyle(20);
        g2g[i]->SetMarkerSize(1.1);
        g2g[i]->SetMarkerColor(1);
        f2g[i]->Draw("SAME");
        g2g[i]->Draw("P SAME");
        c1->Update();

        cout<<"ch="<<i<<"\ta0="<<a0<<"\ta1="<<a1<<endl;
        ofs<<i<<"\t"<<a0<<"\t"<<a1<<"\t"<<a0/100<<"\t"<<a1/100<<endl;//pretent error 1%
        gSystem->Sleep(1000);
    }
}



void refine_fit(Int_t dssd=0,Double_t width_limit=150)
{

    std::ifstream ifslowpoints(Form("dssd%d_beta_lowpoints.txt",dssd));
    Double_t xlow[256];
    Double_t ylow[256];
    std::string tempstr;
    for (Int_t i=0;i<256;i++){
        ifslowpoints>>tempstr>>xlow[i]>>ylow[i];
    }
    ifslowpoints.close();
    std::ifstream ifshighpoints(Form("dssd%d_beta_highpoints.txt",dssd));
    Double_t xhigh[256];
    Double_t yhigh[256];
    for (Int_t i=0;i<256;i++){
        ifshighpoints>>tempstr>>xhigh[i]>>yhigh[i];
    }
    ifshighpoints.close();

    TFile *f=TFile::Open(Form("dssd%d_beta.root",dssd));
    TH2F* h2g[256];
    TGraph* graphin[256];
    TF1* f2g[256];
    TF1* f2g_refine[256];

    TGraph* g2g[256];


    TGraph* refineGraph[256];

    TCanvas* c1=new TCanvas("hisname","hisname",800,600);

    std::ofstream ofs(Form("dssd%d_beta_coefs.txt",dssd));
    for (Int_t i=0;i<256;i++){
        h2g[i]=(TH2F*) f->Get(Form("h%d",i));
        graphin[i]=(TGraph*) f->Get(Form("g%d",i));
        h2g[i]->Draw("colz");
        //cout<<xlow[i]<<"\t"<<ylow[i]<<"\t"<<xhigh[i]<<"\t"<<yhigh[i]<<endl;
        Double_t x[]={xlow[i],xhigh[i]};
        Double_t y[]={ylow[i],yhigh[i]};
        g2g[i]=new TGraph(2,x,y);
        Double_t a0,a1;
        Int_t ifail;
        g2g[i]->LeastSquareLinearFit(2,a0,a1,ifail);
        f2g[i]=new TF1(Form("f%d",i),"pol1",h2g[i]->GetXaxis()->GetXmin(),h2g[i]->GetXaxis()->GetXmax());
        f2g[i]->FixParameter(0,a0);
        f2g[i]->FixParameter(1,a1);
        g2g[i]->SetMarkerStyle(20);
        g2g[i]->SetMarkerSize(1.1);
        g2g[i]->SetMarkerColor(1);
        f2g[i]->Draw("SAME");
        g2g[i]->Draw("P SAME");

        //! refine

        Double_t* xx=graphin[i]->GetX();
        Double_t* yy=graphin[i]->GetY();

        Double_t xnew[graphin[i]->GetN()];
        Double_t ynew[graphin[i]->GetN()];
        Int_t k=0;

        for (Int_t j=0;j<graphin[i]->GetN();j++){
            if ((yy[j]>a0-width_limit+a1*xx[j])&&(yy[j]<a0+width_limit+a1*xx[j])){
                xnew[k]=xx[j];
                ynew[k]=yy[j];
                k++;
            }
        }
        refineGraph[i]=new TGraph(k,xnew,ynew);
        refineGraph[i]->SetName(Form("gnew%d",i));
        refineGraph[i]->Fit("pol1","Q");
        refineGraph[i]->GetFunction("pol1")->SetLineColor(3);
	refineGraph[i]->SetMarkerStyle(20);
	refineGraph[i]->SetMarkerSize(0.2);
        refineGraph[i]->Draw("P SAME");

        /*
        refineGraph[i]->LeastSquareLinearFit(k,a0,a1,ifail);
        f2g_refine[i]=new TF1(Form("fre%d",i),"pol1",h2g[i]->GetXaxis()->GetXmin(),h2g[i]->GetXaxis()->GetXmax());
        f2g_refine[i]->FixParameter(0,a0);
        f2g_refine[i]->FixParameter(1,a1);
        f2g_refine[i]->SetLineColor(3);
        f2g_refine[i]->Draw("SAME");
        */

        cout<<"ch="<<i<<"\ta0="<<refineGraph[i]->GetFunction("pol1")->GetParameter(0)<<"\ta1="<<refineGraph[i]->GetFunction("pol1")->GetParameter(1)
           <<"\ta0_err="<<refineGraph[i]->GetFunction("pol1")->GetParError(0)<<"\ta1_err="<<refineGraph[i]->GetFunction("pol1")->GetParError(1)<<endl;

        c1->Update();
        ofs<<i<<"\t"<<refineGraph[i]->GetFunction("pol1")->GetParameter(0)<<"\t"<<refineGraph[i]->GetFunction("pol1")->GetParameter(1)<<"\t"<<refineGraph[i]->GetFunction("pol1")->GetParError(0)<<"\t"<<refineGraph[i]->GetFunction("pol1")->GetParError(1)<<endl;
        gSystem->Sleep(100);
    }
}

void refine_fit_mul(Int_t dssd=0,Double_t width_limit=150)
{

    std::ifstream ifslowpoints(Form("dssd%d_beta_lowpoints.txt",dssd));
    Double_t xlow[256];
    Double_t ylow[256];
    std::string tempstr;
    for (Int_t i=0;i<256;i++){
        ifslowpoints>>tempstr>>xlow[i]>>ylow[i];
    }
    ifslowpoints.close();
    std::ifstream ifshighpoints(Form("dssd%d_beta_highpoints.txt",dssd));
    Double_t xhigh[256];
    Double_t yhigh[256];
    for (Int_t i=0;i<256;i++){
        ifshighpoints>>tempstr>>xhigh[i]>>yhigh[i];
    }
    ifshighpoints.close();

    TFile *f=TFile::Open(Form("dssd%d_beta.root",dssd));
    TH2F* h2g[256];
    TGraph* graphin[256];
    TF1* f2g[256];
    TF1* f2g_refine[256];

    TGraph* g2g[256];


    TGraph* refineGraph[256];

    TCanvas* c1=new TCanvas("hisname","hisname",800,600);

    std::ofstream ofs(Form("dssd%d_beta_coefs.txt",dssd));
    for (Int_t i=0;i<256;i++){
        h2g[i]=(TH2F*) f->Get(Form("h%d",i));
        graphin[i]=(TGraph*) f->Get(Form("g%d",i));
        h2g[i]->Draw("colz");
        //cout<<xlow[i]<<"\t"<<ylow[i]<<"\t"<<xhigh[i]<<"\t"<<yhigh[i]<<endl;
        Double_t x[]={xlow[i],xhigh[i]};
        Double_t y[]={ylow[i],yhigh[i]};
        g2g[i]=new TGraph(2,x,y);
        Double_t a0,a1;
        Int_t ifail;
        g2g[i]->LeastSquareLinearFit(2,a0,a1,ifail);
        f2g[i]=new TF1(Form("f%d",i),"pol1",h2g[i]->GetXaxis()->GetXmin(),h2g[i]->GetXaxis()->GetXmax());
        f2g[i]->FixParameter(0,a0);
        f2g[i]->FixParameter(1,a1);
        g2g[i]->SetMarkerStyle(20);
        g2g[i]->SetMarkerSize(1.1);
        g2g[i]->SetMarkerColor(1);
        f2g[i]->Draw("SAME");
        g2g[i]->Draw("P SAME");

        //! refine

        Double_t* xx=graphin[i]->GetX();
        Double_t* yy=graphin[i]->GetY();

        Double_t xnew[graphin[i]->GetN()];
        Double_t ynew[graphin[i]->GetN()];
        Int_t k=0;

        for (Int_t j=0;j<graphin[i]->GetN();j++){
            if ((yy[j]>a0-width_limit+a1*xx[j])&&(yy[j]<a0+width_limit+a1*xx[j])){
                xnew[k]=xx[j];
                ynew[k]=yy[j];
                k++;
            }
        }

        refineGraph[i]=new TGraph(k,xnew,ynew);
        refineGraph[i]->SetName(Form("gnew%d",i));
        refineGraph[i]->SetMarkerStyle(20);
        refineGraph[i]->SetMarkerSize(0.2);
        refineGraph[i]->Draw("P SAME");

        //double c1,cov11,sumsq;
        //gsl_fit_mul(xnew,1,ynew,1,k,&c1,&cov11,&sumsq);

        f2g_refine[i]=new TF1(Form("fre%d",i),"[0]*x",h2g[i]->GetXaxis()->GetXmin(),h2g[i]->GetXaxis()->GetXmax());
        f2g_refine[i]->SetParameter(0,a1);
        f2g_refine[i]->SetLineColor(3);

        refineGraph[i]->Fit(f2g_refine[i],"Q");
        f2g_refine[i]->Draw("SAME");

        c1->Update();
        ofs<<i<<"\t"<<f2g_refine[i]->GetParameter(0)<<"\t"<<f2g_refine[i]->GetParError(0)<<endl;
        gSystem->Sleep(100);
    }
}
