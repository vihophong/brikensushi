#include "TChain.h"
#include "TLatex.h"

#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH2.h"
#include "TRandom.h"
#include "TString.h"
#include "TLine.h"
#include "TSpectrum.h"
#include "TF1.h"
#include <fstream>


void plotpid(char* listfile="listall.txt",char* pidfile="pidfile.txt"){
    std::ifstream ifspid(pidfile);
    Int_t nri=0;
    Int_t nbinszet,nbinsaoq;
    Double_t zetrange[2];
    Double_t aoqrange[2];
    ifspid>>nbinsaoq>>aoqrange[0]>>aoqrange[1]>>nbinszet>>zetrange[0]>>zetrange[1];
    ifspid>>nri;

    Int_t enablepid[nri];
    Int_t enablepid2[nri];

    TString nameri[nri];
    TString latexnametri[nri];
    Double_t parmsri[nri][7];
    TString tempriname,tempria;
    TCutG* cutg[nri];
    TLatex* pidtag[nri];

    Double_t halflife[nri];

    for (Int_t i=0;i<nri;i++){
        ifspid>>enablepid[i]>>enablepid2[i]>>tempria>>tempriname>>halflife[i];
        for(Int_t j=0;j<7;j++) ifspid>>parmsri[i][j];
        nameri[i]=tempriname+tempria;
        latexnametri[i]=TString("^{")+tempria+TString("}"+tempriname);
        cout<<nameri[i]<<"\t"<<latexnametri[i]<<"\t"<<halflife[i];
        for(Int_t j=0;j<7;j++) cout<<"\t"<<parmsri[i][j];
        cout<<endl;

        pidtag[i]=new TLatex(parmsri[i][0],parmsri[i][1]+0.2,latexnametri[i]);
        pidtag[i]->SetTextSize(0.025);
        pidtag[i]->SetTextColor(2);
    }
    std::ifstream ifs(listfile);
    string filelist[1000];
    Int_t nfiles=0;
    while (!ifs.eof()){
        ifs>>filelist[nfiles];
        cout<<filelist[nfiles]<<endl;
        nfiles++;
    }

    Int_t ncutpts=20;// number of cut points
    for (Int_t i=0;i<nri;i++){
        cutg[i]=new TCutG(nameri[i],ncutpts);
        cutg[i]->SetTitle(nameri[i]);
        cutg[i]->SetVarX("decayaoq");
        cutg[i]->SetVarY("decay.zet");

        Double_t theta=0;
        Double_t x1=parmsri[i][0];
        Double_t y1=parmsri[i][1];
        Double_t r1=parmsri[i][2];
        Double_t r2=parmsri[i][3];

        Double_t phi1 = TMath::Min(0,360);
        Double_t phi2 = TMath::Max(0,360);
        Double_t kPI = 3.14159265358979323846;
        Double_t angle,dx,dy;

        Double_t x[ncutpts], y[ncutpts];
        Double_t dphi = (phi2-phi1)*kPI/(180*ncutpts);
        Double_t ct   = TMath::Cos(kPI*theta/180);
        Double_t st   = TMath::Sin(kPI*theta/180);

        for (Int_t j=0;j<ncutpts;j++) {
           angle = phi1*kPI/180 + Double_t(j)*dphi;
           dx    = r1*TMath::Cos(angle);
           dy    = r2*TMath::Sin(angle);
           x[j]  = x1 + dx*ct - dy*st;
           y[j]  = y1 + dx*st + dy*ct;
           //cout<<x[j]<<"<"<<y[j]<<endl;
           cutg[i]->SetPoint(j,x[j],y[j]);
        }
        cutg[i]->SetPoint(ncutpts,x[0],y[0]);
    }

    TChain* ch = new TChain("treeimp");
    nfiles=nfiles-1;
    cout<<"There are "<<nfiles<<" files in total!"<<endl;

    for (Int_t i=0;i<nfiles;i++){
        char tempchar2[1000];
        sprintf(tempchar2,"%s/treeimp",filelist[i].c_str());
        ch->Add(tempchar2);
    }
    TCanvas* c1=new TCanvas("pid","pid",900,700);
    gStyle->SetOptStat(1111111);
    ch->Draw(Form("zet:aoq>>hpid(%d,%f,%f,%d,%f,%f)",nbinsaoq,aoqrange[0],aoqrange[1],nbinszet,zetrange[0],zetrange[1]),"zet>0&&aoq>0&&ion.z>=0","colz");
    TH1F* hpid=(TH1F*)gDirectory->Get("hpid");
    hpid->SaveAs("pid.root");

    for (Int_t i=0;i<nri;i++){
        cutg[i]->Draw("same");
        pidtag[i]->Draw("same");
    }
}


void projectpid(char* infile,Double_t Zcenter,Double_t width,char* parmsfile, char* outfile,Int_t opt=0)
{
    TFile* file1=TFile::Open(infile);
    TH2F* hpid=(TH2F*)file1->Get("hpid");
    hpid->Draw("colz");
    Int_t startbin=hpid->GetYaxis()->FindBin(Zcenter-width/2);
    Int_t stopbin=hpid->GetYaxis()->FindBin(Zcenter+width/2);
    cout<<startbin<<"-"<<stopbin<<endl;
    TH1F* hproj=(TH1F*) hpid->ProjectionX("prj",startbin,stopbin);
    TCanvas* c1=new TCanvas("c1","c1",900,700);
    c1->SetLogz();
    TLine* l1=new TLine(hpid->GetXaxis()->GetXmin(), Zcenter-width/2,hpid->GetXaxis()->GetXmax(),Zcenter-width/2);
    TLine* l2=new TLine(hpid->GetXaxis()->GetXmin(), Zcenter+width/2,hpid->GetXaxis()->GetXmax(),Zcenter+width/2);
    l1->SetLineColor(2);
    l2->SetLineColor(2);
    hpid->Draw("colz");
    l1->Draw("same");
    l2->Draw("same");

    TCanvas* c2=new TCanvas("c2","c2",900,700);
    hproj->Draw();
    TSpectrum *s=new TSpectrum();
    s->Search(hproj,10,"",0.001);

    //! fitting
    //! single component fit
//    TF1 *fa = new TF1("fa","[0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2]))",2.738,2.75);
//    fa->SetParameter(0,7.31852e+02);
//    fa->SetParameter(1,2.74145e+00);
//    fa->SetParameter(2,1.34367e-03);

    //! double component fit
    TF1* ffit[100];
    TF1* fcontamination[100];
    Int_t is_one_peak[100];
    TString namefit[100];
    Double_t par0[100];
    Double_t par1[100];
    Double_t par2[100];
    Double_t par3[100];
    Double_t par4[100];

    Double_t fitrangelow[100];
    Double_t fitrangeup[100];



    std::ifstream ifs;
    ifs.open(parmsfile);



    Int_t nfits=0;


    while (!ifs.eof()){
        ifs>>namefit[nfits]>>is_one_peak[nfits]>>par0[nfits]>>par1[nfits]>>par2[nfits]>>par3[nfits]>>par4[nfits]>>fitrangelow[nfits]>>fitrangeup[nfits];
        cout<<par4[nfits]<<endl;
        nfits++;
    }
    nfits=nfits-1;
    cout<<"There are "<<nfits<<" lines in input file"<<endl;



/*
    nfits=3;
    is_one_peak[0]=0;
    namefit[0]="Sn135";
    par0[0]=10;
    par1[0]=2.734;
    par2[0]=7.31852e+02;
    par3[0]=2.74145e+00;
    par4[0]=1.36919e-03;
    fitrangelow[0]=2.732;
    fitrangeup[0]=2.75;

    is_one_peak[1]=0;
    namefit[1]="Sn136";
    par0[1]=10;
    par1[1]=2.714;
    par2[1]=2500;
    par3[1]=2.721;
    par4[1]=1.36919e-03;
    fitrangelow[1]=2.711;
    fitrangeup[1]=2.727;

    is_one_peak[2]=1;
    namefit[2]="Sn137";
    par0[2]=180;
    par1[2]=2.701;
    par2[2]=0;
    par3[2]=0;
    par4[2]=1.36919e-03;
    fitrangelow[2]=2.696;
    fitrangeup[2]=2.706;
    */


    std::ofstream str;
    if (opt==0){
        str.open(outfile);
    }else{
        str.open(outfile,std::ofstream::out | std::ofstream::app);
    }

    for (Int_t i=0;i<nfits;i++){
        ffit[i]=new TF1(Form("fit_%s",namefit[i].Data()),"[0]*exp(-0.5*((x-[1])/[4])*((x-[1])/[4]))+[2]*exp(-0.5*((x-[3])/[4])*((x-[3])/[4]))",fitrangelow[i],fitrangeup[i]);
        fcontamination[i]=new TF1(Form("fitcont_%s",namefit[i].Data()),"gaus",fitrangelow[i],fitrangeup[i]);

        ffit[i]->SetParameter(0,par0[i]);
        ffit[i]->SetParameter(1,par1[i]);

        if (is_one_peak[i]!=0){
            ffit[i]->FixParameter(2,par2[i]);
            ffit[i]->FixParameter(3,par3[i]);
        }else{
            ffit[i]->SetParameter(2,par2[i]);
            ffit[i]->SetParameter(3,par3[i]);
        }

        ffit[i]->SetParameter(4,par4[i]);
        hproj->Fit(ffit[i],"LER+");

        Double_t nsigma=3.;
        fcontamination[i]->SetParameter(0,ffit[i]->GetParameter(2));
        fcontamination[i]->SetParameter(1,ffit[i]->GetParameter(3));
        fcontamination[i]->SetParameter(2,ffit[i]->GetParameter(4));
        fcontamination[i]->SetLineColor(3);
        fcontamination[i]->Draw("same");
        //! scan for optimum charge state contaminant

        //if (i==nfits-1) nsigma=2.1; //for Ag131
        //if (i==nfits-2) nsigma=2.1; //for Cd133
        //if (i==nfits-1) nsigma=1.05; //for Cd134
        //if (i==nfits-2) nsigma=1.4; //for In136
        //if (i==nfits-3) nsigma=2.2; //for Sn138
        //if (i==nfits-2) nsigma=0.55; //for Sn139


        Double_t rightgate=ffit[i]->GetParameter(1)+nsigma*ffit[i]->GetParameter(4);
        Double_t leftgate=ffit[i]->GetParameter(1);
        Double_t ncontamination=fcontamination[i]->Integral(leftgate,rightgate);
        Double_t ngated=ffit[i]->Integral(leftgate,rightgate);


        //! optimization of left cut
        for (Int_t k=0;k<1001;k++){
            Double_t lsigma=k*nsigma/1000.;
            Double_t ileftgate=ffit[i]->GetParameter(1)-lsigma*ffit[i]->GetParameter(4);
            Double_t incontamination=fcontamination[i]->Integral(ileftgate,rightgate);
            Double_t ingated=ffit[i]->Integral(ileftgate,rightgate);
            if (incontamination/ingated>0.01) {
                //break;
            }
            leftgate=ileftgate;
            ncontamination=incontamination;
            ngated=ingated;
        }



        TLine* l1=new TLine(leftgate,0,leftgate,100000);
        TLine* l2=new TLine(rightgate,0,rightgate,100000);
        l1->Draw();
        l2->Draw();

        str<<namefit[i]<<"\t"<<is_one_peak[i]<<"\t"<<ffit[i]->GetParameter(0)<<"\t"<<ffit[i]->GetParameter(1)<<"\t"<<ffit[i]->GetParameter(2)<<"\t"<<ffit[i]->GetParameter(3)<<"\t"<<ffit[i]->GetParameter(4)<<"\t"<<ffit[i]->GetParError(4)<<"\t"<<fitrangelow[i]<<"\t"<<fitrangeup[i]<<"\t"<<leftgate<<"\t"<<rightgate<<"\t"<<ncontamination/ngated*100.<<endl;
    }


    /*
    TF1 *fa = new TF1("fa","[0]*exp(-0.5*((x-[1])/[4])*((x-[1])/[4]))+[2]*exp(-0.5*((x-[3])/[4])*((x-[3])/[4]))",2.732,2.75);
    fa->SetParameter(0,10);
    fa->SetParameter(1,2.734);
    fa->SetParameter(2,7.31852e+02);
    fa->SetParameter(3,2.74145e+00);
    fa->SetParameter(4,1.34367e-03);
    hproj->Fit(fa,"LER+");
    */

}

void fitHZ(TH1D* hzin){

    hzin->Draw();

    //high in
    Double_t zpeakh[]={88,2290,11551,20057,18715,9480,4024,9341,50564,47973,19731,7995,2626,1319,175};

    TF1* fa=new TF1("fa","gaus(0)+gaus(3)+gaus(6)+gaus(9)+gaus(12)+gaus(15)+gaus(18)+gaus(21)+gaus(24)+gaus(27)+gaus(30)+gaus(33)+gaus(36)+gaus(39)+gaus(42)",31.4,46.6);
    for (Int_t i=0;i<15;i++){
        fa->SetParameter(i*3,zpeakh[i]);
        fa->SetParameter(i*3+1,i+32);
        fa->SetParameter(i*3+2,0.19);
    }

    fa->SetNpx(2000);
    fa->Draw("same");
    hzin->Fit(fa,"LER+");
    for (Int_t i=0;i<15;i++){
        cout<<fa->GetParameter(i*3+1)<<"\t"<<fa->GetParameter(i*3+2)<<"\t"<<fa->GetParError(i*3+2)<<endl;
    }

}
