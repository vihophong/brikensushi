#include "libDataStruct.h"

#include "TH1.h"
#include "TH2.h"

#ifndef SORTER_h
#define SORTER_h 1
class sorter
{
public:
  sorter();
  virtual ~sorter();

  void setCurrentDgtzHits(NIGIRI *hitdgtz1, NIGIRI *hitdgtz2);

  void doSortLowGain();
  void doSortHighGain();

  void printCurrent();

  void initHistos();
  void writeHistos();

private:
  NIGIRI* fhitDgtz1;
  NIGIRI* fhitDgtz2;

  NIGIRI* fcurrentHitDgtz1;
  NIGIRI* fcurrentHitDgtz2;

  TH1F* hh1;
};
#endif
