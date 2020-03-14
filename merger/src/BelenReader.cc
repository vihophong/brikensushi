#include "BelenReader.hh"

void genRndCircle(Double_t &x,Double_t &y,Double_t a,Double_t b,Double_t xpos,Double_t ypos,Double_t R){
    if (b<a){
        Double_t temp=a;
        a=b;
        b=temp;
    }
    x=xpos+b*R*TMath::Cos(2*TMath::Pi()*a/b);
    y=ypos+b*R*TMath::Sin(2*TMath::Pi()*a/b);
}

BelenReader::BelenReader():rr()
{

    for (Int_t i=0;i<MaxID;i++){
        fHe3Ecal[i][0]=1.;
        fHe3Ecal[i][1]=1.;
        fHe3Id2posX[i]=0;
        fHe3Id2posY[i]=0;
        fHe3Id2posZ[i]=0;
        fHe3Id2diameter[i]=0;
        fHe3Id2ring[i]=0;
        fHe3Id2length[i]=0;
    }

    for (Int_t i=0;i<MaxIndex1;i++){
        for (Int_t j=0;j<MaxIndex2;j++){
            fCrystalId2posX[i][j]=0;
            fCrystalId2posY[i][j]=0;
            fCrystalId2posZ[i][j]=0;
        }
    }

    ftreedataNeuron = NULL;
    ftreedataGamma = NULL;
    ftreedataAnc = NULL;

    ftreedataYSO = NULL;

    fflag_filldata = false;
}

BelenReader::~BelenReader()
{
    if (!fflag_filldata){
        delete flocalNeutron;
        delete flocalGamma;
        delete flocalAnc;
    }
    delete finfile;
}


void BelenReader::Init(char* belenfile){
    if (!fflag_filldata){
        flocalNeutron = new BELENHit;
        flocalGamma = new CloverHit;
        flocalAnc = new BELENHit;
    }
    fBLAncEntry = 0;
    fBLGamEntry = 0;
    fBLAncEntry = 0;

    finfile = new TFile(belenfile);
    ftree = (TTree*) finfile->Get("BRIKENTree");
    finfile->cd();
    fnentries = ftree->GetEntries();
    cout<<"There are "<<fnentries<<" entries in Belen: "<< belenfile<<endl;

    //! brach tree

    ftree->SetBranchAddress("Neutrons.",&ftreedataNeuron);
    ftree->SetBranchAddress("Gamma.",&ftreedataGamma);
    ftree->SetBranchAddress("Ancillary.",&ftreedataAnc);


    ftreeYSO = (TTree*) finfile->Get("YSOTree");
    finfile->cd();
    fnentriesYSO = ftreeYSO->GetEntries();
    cout<<"There are "<<fnentriesYSO<<" entries in YSO"<<endl;

    ftreeYSO->SetBranchAddress("Yso.",&ftreedataYSO);

    fcurentry = 0;
    fcurentryYSO=0;

    fBLNeuEntry = 0;
    fBLGamEntry = 0;
    fBLAncEntry = 0;
    fBLYSOEntryIon = 0;
    fBLYSOEntryBeta = 0;
    GetGeoMapping();
    //if (!GetNextEvent()) exit(1);
}

void BelenReader::ClearAncHits(){
    /*s
    for (unsigned int idx=0;idx<flocalAncUpstreamPL.size();idx++){
        delete flocalAncUpstreamPL[idx];
    }
    for (unsigned int idx=0;idx<flocalAncAIDAPL.size();idx++){
        delete flocalAncAIDAPL[idx];
    }
    for (unsigned int idx=0;idx<flocalAncdE.size();idx++){
        delete flocalAncdE[idx];
    }
    for (unsigned int idx=0;idx<flocalAncF11PL.size();idx++){
        delete flocalAncF11PL[idx];
    }
    flocalAncUpstreamPL.clear();
    flocalAncAIDAPL.clear();
    flocalAncdE.clear();
    flocalAncF11PL.clear();
    */
    flocalAnc->Clear();
}

void BelenReader::GetGeoMapping(){
    std::ifstream inpf(fmappingfile);
    if (inpf.fail()){
        cout<<"No BELEN Mapping file is given"<<endl;
        return;
    }
    cout<<"Start reading BELEN Mapping file: "<<fmappingfile<<endl;

    Int_t id,index1,index2;
    UShort_t ring;
    Double_t x,y,z;
    Double_t d,length;
    Double_t ftoadc,ecal1;
    Int_t mm=0;

    while (inpf.good()){
        inpf>>id>>index1>>index2>>d>>x>>y>>z>>ring>>length>>ftoadc>>ecal1;
        if (id<=500){//for he3
            fHe3Ecal[id][0]=ftoadc;
            fHe3Ecal[id][1]=ecal1;
            fHe3Id2posX[id]=x;
            fHe3Id2posY[id]=y;
            fHe3Id2posZ[id]=z;
            fHe3Id2diameter[id]=d;
            fHe3Id2ring[id]=ring;
            fHe3Id2length[id]=length;
        }else if(id>500){ //for clover
            fCrystalId2posX[index1][index2]=x;
            fCrystalId2posY[index1][index2]=y;
            fCrystalId2posZ[index1][index2]=z;
        }
        mm++;
    }
    cout<<"Read "<<mm<<" line"<<endl;
    inpf.close();
}

void BelenReader::BookTree(TTree* treeNeutron, TTree *treeGamma, TTree *treeAnc, BELENHit* neutron,CloverHit* gamma, BELENHit* anc)
{
    //! initilize output
    fmtrNeutron = treeNeutron;
    fmtrGamma = treeGamma;
    fmtrAnc = treeAnc;
    flocalNeutron = neutron;
    flocalGamma = gamma;
    flocalAnc = anc;
    fflag_filldata=true;
}

void BelenReader::BookYSOTree(TTree* treeYSOion,TTree* treeYSObeta)
{
    fmtrYSOion = treeYSOion;
    fmtrYSOion->Branch("ysoion",&ftreedataYSO);
    fmtrYSOion->BranchRef();
    fmtrYSObeta = treeYSObeta;
    fmtrYSObeta->Branch("ysobeta",&ftreedataYSO);
    fmtrYSObeta->BranchRef();
}


bool BelenReader::GetNextEvent(){
    ftree->GetEntry(fcurentry);

    fE = ftreedataNeuron->E + ftreedataGamma->E + ftreedataAnc->E;
    fT = ftreedataNeuron->T + ftreedataGamma->T + ftreedataAnc->T;
    fId = ftreedataNeuron->Id + ftreedataGamma->Id + ftreedataAnc->Id;
    ftype = ftreedataNeuron->type + ftreedataGamma->type + ftreedataAnc->type;
    fIndex1 = ftreedataNeuron->Index1 + ftreedataGamma->Index1 + ftreedataAnc->Index1;
    fIndex2 = ftreedataNeuron->Index2 + ftreedataGamma->Index2 + ftreedataAnc->Index2;

    fInfoFlag = ftreedataNeuron->InfoFlag + ftreedataGamma->InfoFlag + ftreedataAnc->InfoFlag;

    fName = ftreedataNeuron->Name + ftreedataGamma->Name + ftreedataAnc->Name;

    fcurentry++;

    //! fill data if needed
    if (ftreedataNeuron->type==1){
        flocalNeutron->SetTimestamp(ftreedataNeuron->T);
        flocalNeutron->SetDaqID(ftreedataNeuron->Id);
        Int_t id = atoi(fName.substr(2,3).c_str());
        //Int_t id =ftreedataNeuron->Id+1;
        flocalNeutron->SetID(id);
        flocalNeutron->SetADC(ftreedataNeuron->E/fHe3Ecal[id][0]);
        flocalNeutron->SetEnergy(ftreedataNeuron->E/fHe3Ecal[id][0]*fHe3Ecal[id][1]);
        flocalNeutron->SetRing(fIndex1);
        flocalNeutron->SetType(fIndex2);
        PerturbateHe3(id);
        flocalNeutron->SetRndPos(fposX,fposY,fposZ);
        if (fflag_filldata) fmtrNeutron->Fill();
        fBLNeuEntry++;
    }else if (ftype==2){
        flocalGamma->SetEnergy(ftreedataGamma->E);
        flocalGamma->SetTimestamp(ftreedataGamma->T);
        flocalGamma->SetDaqID(ftreedataGamma->Id);
        if (fIndex1==Index1Clover1) flocalGamma->SetClover(1);
        else if (fIndex1==Index1Clover2) flocalGamma->SetClover(2);
        flocalGamma->SetCloverLeaf(fIndex2);
        flocalGamma->SetID((flocalGamma->GetClover()-1)*4 + flocalGamma->GetCloverLeaf());
        PerturbateClover(fIndex1,fIndex2);
        flocalGamma->SetPos(fposX,fposY,fposZ);
        if (fflag_filldata) fmtrGamma->Fill();
        fBLGamEntry++;
    }
    else if (ftype==3){
        //clear hits first!
        //ClearAncHits();

//        BELENHit* hit = new BELENHit();
//        hit->SetEnergy(fE);
//        hit->SetTimestamp(fT);
//        hit->SetDaqID(fId);
//        hit->SetID(fIndex2);
//        hit->SetRing(0);
//        hit->SetPos(0,0,0);

        flocalAnc->SetEnergy(ftreedataAnc->E);
        flocalAnc->SetTimestamp(ftreedataAnc->T);
        flocalAnc->SetDaqID(ftreedataAnc->Id);
        flocalAnc->SetID(fIndex2);

        Bool_t fflaganc=false;
        if (fIndex1==Index1UPlastic){// not present in current experiment
            flocalAnc->SetPos(0,0,1);
            flocalAnc->SetRing(2);
            flocalAnc->SetType(ScintillatorType);
            fflaganc=true;
            //hit->SetType(ScintillatorType);
            //flocalAncUpstreamPL.push_back(hit);
        }else if(fIndex1==Index1F11){
            flocalAnc->SetRing(1);
            flocalAnc->SetType(ScintillatorType);
            flocalAnc->SetPos(0,0,1);
            fflaganc=true;
            //hit->SetType(ScintillatorType);
            //flocalAncF11PL.push_back(hit);
        }else if (fIndex1==Index1AIDAPL){
            flocalAnc->SetID(5); // after lowgain+highgain F11 plastic
            flocalAnc->SetRing(4);
            flocalAnc->SetType(ScintillatorType);
            flocalAnc->SetPos(0,0,4);
            fflaganc=true;
            //hit->SetType(ScintillatorType);
            //flocalAncAIDAPL.push_back(hit);
        }else if (fIndex1==Index1dE){// not present in current experiment
            flocalAnc->SetRing(3);
            flocalAnc->SetType(SilliconType);
            flocalAnc->SetPos(0,0,3);
            fflaganc=true;
            //hit->SetType(SilliconType);
            //flocalAncdE.push_back(hit);
        }else if (fIndex1==Index1F11LG){
            flocalAnc->SetID(fIndex2+2);//start after f11 high gain
            flocalAnc->SetRing(1);
            flocalAnc->SetType(ScintillatorType);
            flocalAnc->SetPos(0,0,1);
            fflaganc=true;
        }
        if (fId==IdDTPulser){
            flocalAnc->SetID(6); // after lowgain+highgain F11 plastic
            flocalAnc->SetRing(5);
            flocalAnc->SetType(PulserType);
            flocalAnc->SetPos(0,0,5);
            fflaganc=true;
        }else if (fId==IdSyncPulser){// not present in current experiment
            flocalAnc->SetRing(6);
            flocalAnc->SetType(PulserType);
            flocalAnc->SetPos(0,0,6);
            fflaganc=true;
        }
        if (fflag_filldata&&fflaganc) {
            fmtrAnc->Fill();
            fBLAncEntry++;
        }
    }
    if (fcurentry>fnentries) return false;
    return true;
}

bool BelenReader::GetNextNeutronEvent(){
    if (!GetNextEvent()) return false;
    while (ftype!=1){
        if (!GetNextEvent()) return false;
    }
    fBLNeuEntry++;
    return true;
}

bool BelenReader::GetNextGammaEvent(){
    if (!GetNextEvent()) return false;
    while (ftype!=2){
        if (!GetNextEvent()) return false;
    }
    fBLGamEntry++;
    return true;
}

bool BelenReader::GetNextAncEvent(){
    if (!GetNextEvent()) return false;
    while (ftype!=3){
        if (!GetNextEvent()) return false;
    }
    fBLAncEntry++;
    return true;
}

bool BelenReader::GetNextYSOEvent(){
    ftreeYSO->GetEntry(fcurentryYSO);
    //! reco YSO here!
    if (ftreedataYSO->ID==4){
        fmtrYSOion->Fill();
        fBLYSOEntryIon++;
    }else{
        fmtrYSObeta->Fill();
        fBLYSOEntryBeta++;
    }

    fcurentryYSO++;
    if (fcurentryYSO>fnentriesYSO) return false;
    return true;
}

void BelenReader::PerturbateHe3(UShort_t He3Id){
    //!there is nothing here for the moment!
    fposX = fHe3Id2posX[He3Id];
    fposY = fHe3Id2posY[He3Id];
    fposZ = fHe3Id2posZ[He3Id];
    flocalNeutron->SetPos(fposX,fposY,fposZ);
    //! pertubating
    Double_t a,b,x,y,r;
    r=fHe3Id2diameter[He3Id]/2;
    a=rr.Rndm();
    b=rr.Rndm();
    genRndCircle(x,y,a,b,fposX,fposY,r);
    fposX = x;
    fposY = y;
    fposZ = rr.Rndm()*fHe3Id2length[He3Id]+fposZ-fHe3Id2length[He3Id]/2;
}

void BelenReader::PerturbateClover(UShort_t Index1, UShort_t Index2){
    //!there is nothing here for the moment!
    fposX = fCrystalId2posX[Index1][Index2];
    fposY = fCrystalId2posZ[Index1][Index2];
    fposZ = fCrystalId2posZ[Index1][Index2];
    //! pertubating
}
