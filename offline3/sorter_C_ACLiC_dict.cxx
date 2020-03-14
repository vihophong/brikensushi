// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dIdata01dIusersdatadIphongdIprojectsdIbrikenfall2017dIbrikensushidIoffline3dIsorter_C_ACLiC_dict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "/data01/usersdata/phong/projects/brikenfall2017/brikensushi/offline3/./sorter.C"

// Header files passed via #pragma extra_include

namespace {
  void TriggerDictionaryInitialization_sorter_C_ACLiC_dict_Impl() {
    static const char* headers[] = {
"./sorter.C",
0
    };
    static const char* includePaths[] = {
"/opt/cernroot/root_v6.08.00/include",
"/opt/cernroot/root_v6.08.00/etc",
"/opt/cernroot/root_v6.08.00/etc/cling",
"/opt/cernroot/root_v6.08.00/include",
"/usr/lib/gcc/x86_64-redhat-linux/4.8.5/../../../../include/c++/4.8.5",
"/usr/lib/gcc/x86_64-redhat-linux/4.8.5/../../../../include/c++/4.8.5/x86_64-redhat-linux",
"/usr/lib/gcc/x86_64-redhat-linux/4.8.5/../../../../include/c++/4.8.5/backward",
"/data01/usersdata/phong/projects/brikenfall2017/brikensushi/offline3/",
"/opt/cernroot/root_v6.08.00/include",
"/data01/usersdata/phong/projects/brikenfall2017/brikensushi/offline3/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "sorter_C_ACLiC_dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "sorter_C_ACLiC_dict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif
#ifndef __ACLIC__
  #define __ACLIC__ 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "./sorter.C"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"sorter", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("sorter_C_ACLiC_dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_sorter_C_ACLiC_dict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_sorter_C_ACLiC_dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_sorter_C_ACLiC_dict() {
  TriggerDictionaryInitialization_sorter_C_ACLiC_dict_Impl();
}
