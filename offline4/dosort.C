#include "stdio.h"
#include "string.h"
void dosort(char* input,char* output) {
  gROOT->ProcessLine(".L aidaclass.h+");
  gROOT->ProcessLine(".L wasabiclass.h+");
  gROOT->ProcessLine(".L sorter.C+");
  gROOT->ProcessLine(Form("sorter(\"%s\",\"%s\");",input,output));
}
