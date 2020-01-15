
#include <iostream>
#include <pmonitor/pmonitor.h>
#include "myanalysis.h"

#include <TH1.h>
#include <TH2.h>

int init_done = 0;

using namespace std;

TH1F *h1;
TH2F *h2;


int pinit()
{

  if (init_done) return 1;
  init_done = 1;

  h1 = new TH1F ( "h1","test histogram", 400, -50, 50);
  h2 = new TH2F ( "h2","test histogram 2D", 100, -50.5, 49.5, 100, -500, 500);

  return 0;

}

int process_event (Event * e)
{

  Packet *p = e->getPacket(1003);
  if (p)
    {

      //  h1->Fill ( p->iValue(3)/1000. );
      //  h2->Fill ( p->iValue(0), p->iValue(1) );

      delete p;

    }
  return 0;
}

