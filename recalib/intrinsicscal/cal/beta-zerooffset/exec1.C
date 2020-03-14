#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TFile.h"
#include "Riostream.h"
#include <TROOT.h>
#include "TExec.h"
#include "TMarker.h"
void exec1()
{
  //example of macro called when a pad is redrawn
  //one must create a TExec object in the following way
  // TExec ex("ex",".x exec1.C");
  // ex.Draw();
  // this macro prints the bin number and the bin content when one clicks
  //on the histogram contour of any histogram in a pad
  //Author: Rene Brun
  ofstream str("result_endpoints.txt",ios::app);
  int event = gPad->GetEvent();

  //if (event != 11) return;
  int px = gPad->GetEventX();
  int py = gPad->GetEventY();
  TObject *select = gPad->GetSelected();
  if (!select) return;
  if (select->InheritsFrom("TH2")) {
    TH2 *h = (TH2*)select;
    Float_t xx = gPad->AbsPixeltoX(px);
    Float_t x  = gPad->PadtoX(xx);
    Float_t yy = gPad->AbsPixeltoY(py);
    Float_t y  = gPad->PadtoX(yy);
    TMarker marker;
    marker.SetMarkerStyle(20);
    marker.SetMarkerColor(2);
    marker.SetMarkerSize(1.5);
    marker.DrawMarker(x,y);
    printf("event=%d, hist:%s, x=%f, content=%f\n",event,h->GetName(),x,y);
    str<<h->GetName()<<"\t"<<x<<"\t"<<y<<endl;
  }
  str.close();
}
