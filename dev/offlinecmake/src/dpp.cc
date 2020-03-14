#include <TF1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TString.h>
#include <TGaxis.h>
#include <TRandom.h>
#include <stdio.h>
#include <stdlib.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TSpline.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string.h>
#include <math.h>
#include <time.h>
#include <memory.h>
#include <unistd.h>
#include <sys/time.h>

#include "dpp.h"
using namespace std;

dpp::dpp(int fnsize,unsigned short * fpulse)
{
  PZC=1;
  NZC=0;
  armed=-1;
  imax=0;
  imin=0;
  overFlow=false;
  nsize=fnsize;
  pulse=fpulse;
  cfd_ret=new double[nsize];
  icfd=new double[nsize];
  pulseD=new double[nsize];
  bL=fpulse[0];

  fineLED=0;
  pulsePolarity=0; // default negative
}
dpp::~dpp()
{
  delete[]cfd_ret;
  delete[]icfd;
  delete[]pulseD;
}
//----Invert pulse
void dpp::invertPositivePulse()
{
  for (int i=0;i<nsize;i++){
    pulse[i]=16384-pulse[i];
  }
  bL=16384-bL;
}
//----Determine the baseline
void dpp::baselineMean(int fnBaseline,int mVB)
{
  double sum=0;
  double mean;
  int temp=pulse[0];
  for (int i=0;i<nsize;i++){
    if ((abs(pulse[i]-temp)<mVB)&&(i<fnBaseline-1)){
      sum+=pulse[i];
      temp=pulse[i];
    }else {
      if (i!=0) mean=sum/(double)i;//avoid error divide by 0
      else mean=pulse[i];
      break;
    }
  }
  bL=mean;
}
// plot the pulse
int  dpp::plotPulse()
{
  ofstream abc("temp.dat");
  for (int i=0;i<nsize;i++){
    abc<<i<<" "<<pulse[i]<<"\n";
  }
  abc.close();
  FILE* gnuplotPipe=popen("gnuplot -persist","w");
  fprintf(gnuplotPipe,"set xrange [0:500]\n");
  fprintf(gnuplotPipe,"set yrange [0:16000]\n");
  fprintf(gnuplotPipe,"plot 'temp.dat' u 1:2 w p pt 1 ps 0.6 lc 0\n");
  return 0;
}
int  dpp::plotCFD()
{
  ofstream abc("temp.dat");
  for (int i=0;i<nsize;i++){
    abc<<i<<" "<<cfd_ret[i]<<"\n";
  }
  abc.close();
  FILE* gnuplotPipe=popen("gnuplot -persist","w");
  fprintf(gnuplotPipe,"set xrange [0:500]\n");
  fprintf(gnuplotPipe,"set yrange [-5000:5000]\n");
  fprintf(gnuplotPipe,"plot 'temp.dat' u 1:2 w p pt 1 ps 0.6 lc 0\n");
  return 0;
}

void dpp::makeCFD(int LED_threshold,int delay,double fraction)
{
  if (pulsePolarity==1) bL=16384-bL;

  int min=0;
  int max=0;
  int armed_thres=bL-LED_threshold;
  int att,invertDelay; //attenuated signal,invert delayed signal
  //ofstream abc("tempcfd.dat");
  bool stopZC=false;


  if (pulsePolarity!=1) ledMin=pulse[0]; else ledMin=16384-pulse[0];

  for (int i=0;i<nsize;i++){
    if (pulsePolarity==1) pulse[i]=16384-pulse[i];

    //Search for minimum simultaneosly
    if (pulse[i]<ledMin) ledMin=pulse[i];

    //Search for LED simultaneosly
      double h,l;
      if ((bL-pulse[i]>LED_threshold)&&(i>0)){
        h=(double)(bL-pulse[i]-LED_threshold);
        l=(double)(LED_threshold-bL+pulse[i-1]);
        fineLED=(double)(i-1)+l/(h+l);
      }


    pulseD[i]=(double)pulse[i];
    //overFlowCheck
    if (pulse[i]==0) overFlow=true;

    att=round((pulse[i]-bL)*fraction+bL);//attenualted
    if (i<delay){
      invertDelay=bL;
    }else{
      invertDelay=pulse[i-delay];//delayed
      invertDelay=bL-(invertDelay-bL);//inverted
    }
    cfd_ret[i]=(double)(invertDelay+att-2*bL);
    if (pulse[i]<armed_thres) armed=i;
    if (stopZC==false&&armed!=-1&&cfd_ret[i]>0){
      PZC=i;
      NZC=i-1;
      stopZC=true;
    }
    if (cfd_ret[i]<min){
      min=cfd_ret[i];
      imin=i;
    }
    if (cfd_ret[i]>max) {
      max=cfd_ret[i];
      imax=i;
    }
    icfd[i]=i;
  }
}

//----------- Constant fraction discrimination
double  dpp::cfdFast()
{
  double fineCfd=0;
  if (PZC==1&&NZC==0){
    fineCfd=0;
    }else{
    fineCfd=(double)(NZC)-(double)cfd_ret[NZC]/((double)cfd_ret[PZC]-(double)cfd_ret[NZC]);
  }
 return fineCfd;
}
double  dpp::interpolationFunc(double* x,double* par){
  return ip->Eval(x[0]);
}
double dpp::cfdSlow()
{
  if (armed==-1) return 0;
  double fine_cfd;
  ip=new TSpline3("Cubic Spline",icfd,cfd_ret,nsize,"b2e2",0,0);
  sp=new TF1("f",this,&dpp::interpolationFunc,0,nsize,1,"dpp","interpolationFunc");
  fine_cfd=sp->GetX(0.,NZC,PZC,1e-10,1000);
  delete ip;
  delete sp;
  return fine_cfd;
}
double* dpp::getInp(int pNum)
{
  double fine_cfd[pNum];
  ip=new TSpline3("Cubic Spline",icfd,cfd_ret,nsize,"b2e2",0,0);
  double x=0;
  for (int i=0;i<pNum;i++){
    x+=(double)nsize/(double)pNum;
    fine_cfd[i]=ip->Eval(x);
  }
  delete ip;
  return fine_cfd;
}
//----- Leading edge discrimination with pulse size correction (run after make
double dpp::ledWithCorr(int minthres)
{
    double fine_led=0;
    double h,l;
    double corrThres=(double)(bL-ledMin)*0.3;
    for (int i=0;i<nsize;i++){
      if (((double)(bL-pulse[i])>corrThres)&&(bL-pulse[i]>minthres)&&(i>0)){
        h=(double)(bL-pulse[i])-corrThres;
        l=(double)(corrThres-bL+pulse[i-1]);
        fine_led=(double)(i-1)+l/(h+l);
        break;
      }
    }
    return fine_led;
}
//----- Leading edge discrimination
double  dpp::led(int threshold)
{
  double fine_led=0;
  double h,l;
  if (pulsePolarity==1) bL=16384-bL;
//cout<<endl;
  for (int i=0;i<nsize;i++){
    //cout<<i<<" "<<pulse[i]<<endl;
    if (pulsePolarity==1) pulse[i]=16384-pulse[i];
    if ((bL-pulse[i]>threshold)&&(i>0)){
      h=(double)(bL-pulse[i]-threshold);
      l=(double)(threshold-bL+pulse[i-1]);
      fine_led=(double)(i-1)+l/(h+l);
      break;
    }
  }
//cout<<endl;
  return fine_led;
}

//----- Leading edge discrimination for positive pulse
double  dpp::ledPos(int threshold)
{
  double fine_led=0;
  double h,l;
  if (pulsePolarity==1) bL=16384-bL;
//cout<<endl;
  for (int i=0;i<nsize;i++){
    //cout<<i<<" "<<pulse[i]<<endl;
    if (pulsePolarity==1) pulse[i]=16384-pulse[i];
    if ((pulse[i]-bL>threshold)&&(i>0)){
      h=(double)(pulse[i]-bL-threshold);
      l=(double)(threshold-pulse[i-1]+bL);
      fine_led=(double)(i-1)+l/(h+l);
      break;
    }
  }
//cout<<endl;
  return fine_led;
}


int dpp::maxfind(int nsamp) //work only for positive
{
    int max=0;
    for (int i=0;i<nsize/nsamp;i++){
        if (pulse[i]*nsamp>max) max=pulse[i];
    }
    return max;
}
int dpp::minAdcNeg(int nsamp) //work only for negative pulse
{
    int min=pulse[0];
    if (nsamp>nsize) nsamp=nsize;
    for (int i=0;i<nsamp;i++){
        if (pulse[i]<min) min=pulse[i];
    }
    return min;
}
int dpp::maxAdcPos(int nsamp) //work only for positive pulse
{
    int max=0;
    if (nsamp>nsize) nsamp=nsize;
    for (int i=0;i<nsamp;i++){
        if (pulse[i]>max) max=pulse[i];
    }
    return max;
}
int  dpp::emptyFlag(int threshold)
{
  int ledF=1;
  if (pulsePolarity==1) bL=16384-bL;
  for (int i=0;i<nsize;i++){
    if (pulsePolarity==1) pulse[i]=16384-pulse[i];
    if ((bL-pulse[i]>threshold)&&(i>0)){
      ledF=0;
      break;
    }
  }
  return ledF;
}

int  dpp::emptyFlagPositive(int threshold) //no need to change polarity
{
  int ledF=1;
  for (int i=0;i<nsize;i++){
    if ((pulse[i]-bL>threshold)&&(i>0)){
      ledF=0;
      break;
    }
  }
  return ledF;
}

//-----Charge intergration within shortGate and longGate
void  dpp::chargeInter(int& cShort,int& cLong,int fineTS,int fgateOffset,int fshortGate,int flongGate){
  cShort=0;
  cLong=0;
  for (int i=0;i<nsize;i++){
    if ((i>=fineTS-fgateOffset)&&(i<=fineTS-fgateOffset+fshortGate)&&fineTS>0){
      cShort=bL-pulse[i]+cShort;
    }
    if ((i>=fineTS-fgateOffset)&&(i<=fineTS-fgateOffset+flongGate)&&fineTS>0){
      cLong=bL-pulse[i]+cLong;
    }
  }
}
//Individual operation
void dpp::storedCfdSpline(TString filename)
{
  ip=new TSpline3("Cubic Spline",icfd,cfd_ret,nsize,"b2e2",0,0);
  sp=new TF1("f",this,&dpp::interpolationFunc,0,nsize,1,"dpp","interpolationFunc");
  TH1* hist=sp->GetHistogram();
  TFile* f1=new TFile(filename,"RECREATE");
  hist->Write();
  f1->Close();
  delete sp;
  delete ip;
}
void dpp::storePulseSpline(TString filename)
{
  ip=new TSpline3("Cubic Spline",icfd,pulseD,nsize,"b2e2",0,0);
  sp=new TF1("f",this,&dpp::interpolationFunc,0,nsize,1,"dpp","interpolationFunc");
  TH1* hist=sp->GetHistogram();
  TFile* f1=new TFile(filename,"RECREATE");
  hist->Write();
  f1->Close();
  delete sp;
  delete ip;
}
