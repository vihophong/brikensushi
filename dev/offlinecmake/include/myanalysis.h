#ifndef __MYANALYSIS_H__
#define __MYANALYSIS_H__

#include <pmonitor/pmonitor.h>
#include <Event/Event.h>
#include <Event/EventTypes.h>
#include "libDataStruct.h"

int process_event (Event *e); //++CINT 

int mergedata(bool flagend,NIGIRI*data=0);

#endif /* __MYANALYSIS_H__ */
