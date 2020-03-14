#include "stdio.h"
#include "string.h"
void dosort_lightparticle(char* input,char* output) {
  gROOT->ProcessLine(".L aidaclass.h+");
  gROOT->ProcessLine(".L sorter_lightparticle.C+");
  gROOT->ProcessLine(Form("sorter(\"%s\",\"%s\");",input,output));
}
