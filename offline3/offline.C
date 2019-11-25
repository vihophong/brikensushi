#include "stdio.h"
#include "string.h"
void offline(char* input,char* output) {
  gROOT->ProcessLine(".L libmyanalysis.so");
  gROOT->ProcessLine(Form("pfileopen(\"%s\");",input));
  gROOT->ProcessLine(Form("phsave(\"%s\");",output));
  gROOT->ProcessLine("prun();");
}
