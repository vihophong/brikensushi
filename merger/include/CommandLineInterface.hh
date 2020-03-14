////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////
////////////                       GrROOT
////////////
////////////          Purpose:
////////////                   To assist in the analysis of data from
////////////                 the gretina/S800 experimental setup.
////////////                          
////////////          Current Maintainers:
////////////                 Kathrin Wimmer  (wimmer@phys.s.u-tokyo.ac.jp)
////////////                 Eric Lunderberg (lunderberg@nscl.msu.edu)
////////////
////////////          Distribution:
////////////                   Please do not redistribute this software directly.
////////////                   If someone new wants a copy of this software,
////////////                 email one of the maintainers for the download link.
////////////                   This allows us to keep track of who has the software,
////////////                 and properly distribute updates and bug fixes.
////////////                 
////////////          Suggestions:
////////////                   We view the development of the software as a collaborative
////////////                 effort, and as such, welcome and appreciate any suggestions
////////////                 for bug fixes and improvements.
////////////
////////////          Disclaimer:
////////////                 This software is provided as-is, with no warranty.
////////////                 No current or future support is guaranteed for this software.
////////////
////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef __COMMAND_LINE_INTERFACE
#define __COMMAND_LINE_INTERFACE

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


using namespace std;

class CommandLineInterface
{
public:
  CommandLineInterface();
  ~CommandLineInterface(){};

  //main functions to check all flags from command line
  bool CheckFlags(uint,char*[],const bool& Debug = false);

  //functions to add flags
  void Add(const char*);
  void Add(const char*, const char*, bool*);
  void Add(const char*, const char*, char**);
  void Add(const char*, const char*, string*);
  void Add(const char*, const char*, int*);
  void Add(const char*, const char*, size_t*);
  void Add(const char*, const char*, long long*);
  void Add(const char*, const char*, double*, double factor = 1.);
  void Add(const char*, const char*, vector<char*>*);
  void Add(const char*, const char*, vector<string>*);
  void Add(const char*, const char*, vector<int>*);
  void Add(const char*, const char*, vector<long long>*);
  void Add(const char*, const char*, vector<double>*, double factor = 1.);

  friend ostream& operator <<(ostream &,const CommandLineInterface &);

private:
  uint fMaximumFlagLength;
  vector<string> fFlags;
  vector<void*>  fValues;
  uint fMaximumTypeLength;
  vector<string> fTypes;
  uint fMaximumCommentLength;
  vector<string> fComments;
  vector<double> fFactors;
};

#endif
