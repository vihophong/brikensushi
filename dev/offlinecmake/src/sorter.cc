#include "sorter.h"
#include <iostream>

sorter::sorter()
{

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

sorter::~sorter()
{

}
void sorter::initHistos()
{
    hh1=new TH1F("hh1","hh1",200,0,100);
}
void sorter::writeHistos()
{
    hh1->Write();
}

void sorter::setCurrentDgtzHits(NIGIRI *hitdgtz1, NIGIRI *hitdgtz2)
{
    fcurrentHitDgtz1=hitdgtz1;
    fcurrentHitDgtz2=hitdgtz2;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void sorter::printCurrent()
{
    fcurrentHitDgtz1->Print();
    if (fcurrentHitDgtz2==0){
        std::cout<<"No correlation with dgtz 2"<<std::endl;
    }else{
        fcurrentHitDgtz2->Print();
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void sorter::doSortLowGain()
{

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void sorter::doSortHighGain()
{

}
