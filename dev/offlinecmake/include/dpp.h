#include <TF1.h>
#include <TString.h>
#include <TSpline.h>
#ifndef DPP_h
#define DPP_h 1
class dpp
{
public:
  dpp(int,unsigned short*);
  virtual ~dpp();
  void baselineMean(int fnBaseline,int mVB);
  void makeCFD(int LED_threshold,int delay,double fraction);
  int plotPulse();
  int plotCFD();
  double cfdFast();
  double cfdSlow();
  double interpolationFunc(double* x,double* par);
  double led(int threshold);
  double ledPos(int threshold);
  double ledWithCorr(int minthres);
  double* getCFD(){return cfd_ret;}
  double* getInp(int pNum);
  void chargeInter(int& cShort,int& cLong,int fineTS,int fgateOffset,int flongGate,int fshortGate);
  void storedCfdSpline(TString filename);
  int maxfind(int nsample);
  int minAdcNeg(int nsample);
  int maxAdcPos(int nsample);
  void storePulseSpline(TString filenam);
  int emptyFlag(int threshold);
  int emptyFlagPositive(int threshold);
  int nsize;
  unsigned short* pulse;
  int PZC;
  int NZC;
  double bL;
  bool overFlow;

  void invertPositivePulse();
  int pulsePolarity;
  int ledMin;
  double fineLED;
private:
  double* icfd;
  double* pulseD;
  double* cfd_ret;
  int armed;
  int imax;
  int imin;
  TSpline3 * ip;
  TF1* sp;
};
#endif

