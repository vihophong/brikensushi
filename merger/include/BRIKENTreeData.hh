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

class AidaTreeData {
public:
    AidaTreeData( ){};
    ~AidaTreeData(){};


    ULong64_t       T;
    ULong64_t       Tfast;
    Double_t        E;
    Double_t        EX;
    Double_t        EY;
    Double_t        x;
    Double_t        y;
    Double_t        z;
    Int_t         nx;
    Int_t         ny;
    Int_t         nz;
    UChar_t         ID;
    void clear(){
        T=0;
        Tfast=0;
        E=0;
        EX=0;
        EY=0;
        x=0;
        y=0;
        z=0;
        ID=0;

    }
    /// \cond CLASSIMP
    ClassDef(AidaTreeData,1);
    /// \endcond
    ///

};

typedef AidaTreeData YSOData;

#endif
