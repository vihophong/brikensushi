#include "stdio.h"
#include "string.h"
void dofitchannel(Int_t dssd,Double_t ecutminx,Double_t ecutminy,Double_t ecutmax,Int_t npoint,char* outname,Int_t startch=0,Int_t stopch=256) {
  gROOT->ProcessLine("gSystem->Load(\"libMathMore\");");
  gROOT->ProcessLine(".L fitchannel.C");
  gROOT->ProcessLine("fitchannel ee;");
  gROOT->ProcessLine(Form("ee.Loop(%d,%f,%f,%f,%d,\"%s\",%d,%d);",dssd,ecutminx,ecutminy,ecutmax,npoint,outname,startch,stopch));
}
