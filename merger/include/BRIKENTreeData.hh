#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <string.h>

#ifndef BrikenTreeData_hh
#define BrikenTreeData_hh
class BrikenTreeData {
public:
    BrikenTreeData(){}
    ~BrikenTreeData(){}
//    std::vector<double> CorrE;
//    std::vector<double> CorrT;
//    std::vector<uint64_t> CorrTS;
//    std::vector<uint16_t> CorrId;
//    std::vector<uint16_t> CorrType;
    Double_t E;
    ULong_t T;
    UShort_t Id;
    UShort_t type;
    UShort_t Index1;
    UShort_t Index2;
    UShort_t InfoFlag;
    std::string Name;
    std::vector<unsigned short> Samples;

    void clear(){}

    /// \cond CLASSIMP
    ClassDef(BrikenTreeData,1);
    /// \endcond
    ///
};

#endif
