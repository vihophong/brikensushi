#include "Merger.hh"
Merger::Merger():implantdata()
{
    finputWasabi=NULL;
    finputBigrips=NULL;
    finputBriken=NULL;
    fTW_IonBetalow=20000000000;
    fTW_IonBetaup=10000000000;

    fTW_IonPidlow = 8500;
    fTW_IonPidup = -8000;

    fTW_IonF11LRlow=8000;
    fTW_IonF11LRup=-7000;


    fTW_IondElow = 8000;
    fTW_IondEup = -7000;

    fTW_IonDownVetolow = 8000;
    fTW_IonDownVetoup = -7000;

    fTW_PIDGammalow = 400000;
    fTW_PIDGammaup = 400000;

    fTW_IonYSOionlow = 8200;
    fTW_IonYSOionup = -8000;


    fsquareSpatialWindowX=3*3;
    fsquareSpatialWindowY=3*3;

    h1dProbe[0]=new TH1F("h1","h1",2000,-fTW_IonBetaup,fTW_IonBetalow);
    fdataOutWasabiBeta=new wasabiHit;
    fdataOutWasabiIon=new wasabiHit;

    //! pulser stuffs
    for (Int_t i=0;i<141;i++){
        fhpulser[i]=new TH1F(Form("hpulser%d",i+1),Form("hpulser%d",i+1),4000,1500,3500);
        fhpulserall[i]=new TH1F(Form("fhpulserall%d",i+1),Form("fhpulserall%d",i+1),4000,1500,3500);
    }
    fh2deadtime=new TH2F("h2deadtime","h2deadtime with veto using calculated start-stop time",141,0,141,2000,0,20);
    fh1deadtime=new TH1F("h1deadtime","h1deadtime with veto using calculated start-stop time",2000,0,20);

    fh2deadtime2=new TH2F("h2deadtime2","h2deadtime without veto using calculated start-stop time",141,0,141,2000,0,20);
    fh1deadtime2=new TH1F("h1deadtime2","h1deadtime without veto using calculated start-stop time",2000,0,20);

    fh2deadtime3=new TH2F("h2deadtime3","h2deadtime with veto using dtpulser",141,0,141,2000,0,20);
    fh1deadtime3=new TH1F("h1deadtime3","h1deadtime with veto using dtpulser",2000,0,20);

    fh2deadtime4=new TH2F("h2deadtime4","h2deadtime without veto using dtpulser",141,0,141,2000,0,20);
    fh1deadtime4=new TH1F("h1deadtime4","h1deadtime without veto using dtpulser",2000,0,20);

    fh1dtpulser=new TH1F("h1dtpulser","h1dtpulser",2000,0,2000);

    ftreeallRI=0;
    ftreeimplantAll=0;

    isslewcorr=false;


    //! gamma re-calibration parameters (not yet updated)
    //! does not have an effect now!
    Double_t sep[8] = {4.912531E+05,4.998953E+05,5.248885E+05,7.622378E+05,1.108068E+06,1.147523E+06,1.164587E+06,9.553438E+05}; //separation points
    Double_t low_offset[8] = {0.310056,0.602804	,0.594689,0.00575245,-0.409776,-0.739931,-0.231043,-0.410475};//low offset
    Double_t low_gain[8] = {0.00158189,0.00154124,0.00153239,0.000668465,0.000707824,0.000687903,0.000672509,0.000700449};//low gain
    Double_t low_se[8]= {2.5638600000E-12,5.1029500000E-12,8.0936100000E-12,-2.5440800000E-12,-3.2340000000E-12,-4.4067300000E-12,-3.9413800000E-12,-2.6372600000E-12};// low second order
    Double_t high_offset[8] = {1.32609,4.65127,-17.4998,4.21305,-1.68140e+00,-1.13557,-2.79797,6.44019};//high offset
    Double_t high_gain[8] = {0.00158113,0.00153239,0.00159751,0.000659638,7.05388e-04,0.000682542,0.000670503,0.000688148};//high gain
    Double_t high_se[8] = {-9.9227100000E-14,6.6060100000E-12,-5.0294000000E-11,1.7949000000E-12,1.00018e-16,5.6552400000E-13,-3.2623400000E-13,2.7326500000E-12};//high second order

    //! updated for briken2 experiment
    Double_t gain_old[8] =
    {0.001558636,0.0006599344,0.0015311753,0.0015807021,0.0007046135,0.000683305,0.0006685615,0.0006968099};
    Double_t offset_old[8] =
    {-0.6545742444,0.4330534066,2.001278013,1.0994163913,-0.5577366054,-0.6254548603,-0.1966707051,-0.2212697416};

    for (int i=0;i<8;i++){
        fsep[i]=sep[i];
        flow_offset[i]=low_offset[i];
        flow_gain[i]=low_gain[i];
        flow_se[i]=low_se[i];
        fhigh_offset[i]=high_offset[i];
        fhigh_gain[i]=high_gain[i];
        fhigh_se[i]=high_se[i];
        fcgainold[i]=gain_old[i];
        fcoffsetold[i]=offset_old[i];
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Merger::~Merger()
{

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Merger::ReadPID(char* pidfile, Int_t ncutpts){
    std::ifstream ifspid(pidfile);
    ifspid>>nbinsaoq>>aoqrange[0]>>aoqrange[1]>>nbinszet>>zetrange[0]>>zetrange[1];
    ifspid>>nri;
    TString tempriname,tempria;
    for (Int_t i=0;i<nri;i++){
        ifspid>>enablepid[i]>>enablepid2[i]>>tempria>>tempriname>>halflife[i];
        for(Int_t j=0;j<7;j++) ifspid>>parmsri[i][j];
        nameri[i]=tempriname+tempria;
        latexnametri[i]=TString("^{")+tempria+TString("}"+tempriname);
        cout<<nameri[i]<<"\t"<<latexnametri[i]<<"\t"<<halflife[i];
        for(Int_t j=0;j<7;j++) cout<<"\t"<<parmsri[i][j];
        cout<<endl;
    }
     //! Setup PID cut
    for (Int_t i=0;i<nri;i++){
        pidtag[i]=new TLatex(parmsri[i][0],parmsri[i][1]+0.2,latexnametri[i]);
        pidtag[i]->SetTextSize(0.025);
        pidtag[i]->SetTextColor(2);

        cutg[i]=new TCutG(nameri[i],ncutpts);
        cutg[i]->SetTitle(nameri[i]);
        cutg[i]->SetVarX("decay.aoq");
        cutg[i]->SetVarY("decay.zet");

        Double_t theta=0;
        Double_t x1=parmsri[i][0];
        Double_t y1=parmsri[i][1];
        Double_t r1=parmsri[i][2];
        Double_t r2=parmsri[i][3];

        Double_t phi1 = TMath::Min(0,360);
        Double_t phi2 = TMath::Max(0,360);
        Double_t kPI = 3.14159265358979323846;
        Double_t angle,dx,dy;

        Double_t x[ncutpts], y[ncutpts];
        Double_t dphi = (phi2-phi1)*kPI/(180*ncutpts);
        Double_t ct   = TMath::Cos(kPI*theta/180);
        Double_t st   = TMath::Sin(kPI*theta/180);

        for (Int_t j=0;j<ncutpts;j++) {
           angle = phi1*kPI/180 + Double_t(j)*dphi;
           dx    = r1*TMath::Cos(angle);
           dy    = r2*TMath::Sin(angle);
           x[j]  = x1 + dx*ct - dy*st;
           y[j]  = y1 + dx*st + dy*ct;
           //cout<<x[j]<<"<"<<y[j]<<endl;
           cutg[i]->SetPoint(j,x[j],y[j]);
        }
        cutg[i]->SetPoint(ncutpts,x[0],y[0]);
        cutg[i]->SetName((char*)nameri[i].Data());
    }

    for (Int_t i=0;i<nri;i++){
        ftreeRI[i]=new TTree(Form("tree%s",(char*)nameri[i].Data()),Form("tree%s",(char*)nameri[i].Data()));
        ftreeimplantRI[i]=new TTree(Form("treeimp%s",(char*)nameri[i].Data()),Form("treeimp%s",(char*)nameri[i].Data()));
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Merger::Init()
{
    fWasabiBeta = NULL;
    fWasabiIon = NULL;

    fbigrips = NULL;
    fclover = NULL;
    fneutron = NULL;
    fanc = NULL;
    fysoion = NULL;
    fysobeta = NULL;


    fnentriesBigrips = 0;
    fnentriesGamma = 0;
    fnentriesNeutron = 0;
    fnentriesAnc = 0;

    fnentriesYSOIon = 0;
    fnentriesYSOBeta = 0;

    if (finputWasabi!=NULL){
        //! init wasabi
        fWasabiBeta=new Wasabi(finputWasabi,true);
        fnentriesWasabiBeta=fWasabiBeta->GetEntries();
        fWasabiIon=new Wasabi(finputWasabi,false);
        fnentriesWasabiIon=fWasabiIon->GetEntries();

        cout<<"Reading "<<fnentriesWasabiIon<<" ions, "
           <<fnentriesWasabiBeta<<" betas in Wasabi file"<<endl;
        cout<<"Printing first few timestamp:"<<endl;
        for (Long64_t i=0;i<10;i++){
            fWasabiBeta->GetEntry(i);
            fWasabiIon->GetEntry(i);
            cout<<"ionts "<<fWasabiIon->timestamp<<"- betats "<<fWasabiBeta->timestamp<<endl;
        }
    }
    if (finputBriken!=NULL){
        //! init briken
        fBrikenFile = new TFile(finputBriken);
        fBrikenFile->GetObject("neutron",ftrNeutron);
        fBrikenFile->GetObject("gammas",ftrGamma);
        fBrikenFile->GetObject("anc",ftrAnc);
        fBrikenFile->GetObject("ysoion",ftrYSOIon);
        fBrikenFile->GetObject("ysobeta",ftrYSOBeta);
        ftrNeutron->SetBranchAddress("neutron",&fneutron);
        ftrGamma->SetBranchAddress("gamma",&fclover);
        ftrAnc->SetBranchAddress("anc",&fanc);
        ftrYSOIon->SetBranchAddress("ysoion",&fysoion);
        ftrYSOBeta->SetBranchAddress("ysobeta",&fysobeta);

        fnentriesNeutron = ftrNeutron->GetEntries();
        fnentriesGamma = ftrGamma->GetEntries();
        fnentriesAnc = ftrAnc->GetEntries();
        fnentriesYSOIon = ftrYSOIon->GetEntries();
        fnentriesYSOBeta = ftrYSOBeta->GetEntries();

        cout<<"Reading "<<fnentriesNeutron<<" neutrons, "
           <<fnentriesGamma<<" gammas, "<<fnentriesAnc<<" anc hits, "<<fnentriesYSOIon<<" yso ion hits, "<<fnentriesYSOBeta<<" and yso beta hits in AIDA tree"<<endl;
        cout<<"Printing first few timestamp:"<<endl;
        for (Long64_t i=0;i<10;i++){
            ftrNeutron->GetEvent(i);
            ftrGamma->GetEvent(i);
            ftrAnc->GetEvent(i);
            ftrYSOIon->GetEvent(i);
            ftrYSOBeta->GetEvent(i);
            cout<<"neutronts "<<fneutron->GetTimestamp()<<"- gammats "<<fclover->GetTimestamp()<<"- ancts"<<fanc->GetTimestamp()<<"- ysoionts"<<fysoion->T<<"- ysobetats"<<fysobeta->T<<endl;
        }
    }

    if (finputBigrips!=NULL){
        //! init bigrips
        fBigripsFile = new TFile(finputBigrips);
        cout<<finputBigrips<<endl;
        fBigripsFile->GetObject("tree",ftrBigrips);
        ftrBigrips->SetBranchAddress("bigrips",&fbigrips);
        fnentriesBigrips = ftrBigrips->GetEntries();

        cout<<"Reading "<<fnentriesBigrips<<" bigrips items"<<endl;
        cout<<"Printing first few timestamp:"<<endl;
        for (Long64_t i=0;i<10;i++){
            ftrBigrips->GetEvent(i);
            cout<<"bigripsts "<<fbigrips->ts<<endl;
        }
        cout<<"\n\t\t\t****************************************\n"<<endl;
    }


    ftsbeginpulser=0;
    ftsendpulser=0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Merger::ReadWasabi(Int_t opt)
{
    if (opt==0||opt==2){
        for (Long64_t ientry=0;ientry<fnentriesWasabiIon;ientry++){
            fWasabiIon->GetEntry(ientry);
            for (unsigned int i=0;i<fWasabiIon->nhits;i++){
                if (fWasabiIon->z[i]==fWasabiIon->zmax){// last layer fired
                    wasabiHit* ionHit=new wasabiHit;
                    ionHit->ts=fWasabiIon->timestamp;
    #ifdef POS_CORR_CENTER_OF_GRAVITY
                    ionHit->x=fWasabiIon->xcog[i];
                    ionHit->y=fWasabiIon->ycog[i];
    #else
                    ionHit->x=fWasabiIon->xmax[i];
                    ionHit->y=fWasabiIon->ymax[i];
    #endif
                    ionHit->z=fWasabiIon->z[i];
                    ionHit->nx=fWasabiIon->nx[i];
                    ionHit->ny=fWasabiIon->ny[i];
                    ionHit->nz=fWasabiIon->nz;
                    ionHit->zmax=fWasabiIon->zmax;
                    ionHit->ex=fWasabiIon->ex[i];
                    ionHit->ey=fWasabiIon->ey[i];
                    ionHit->blxmax=fWasabiIon->blxmax[i];
                    ionHit->blymax=fWasabiIon->blymax[i];
                    fwasabiIonMap.insert(make_pair(fWasabiIon->timestamp,ionHit));
                }
            }
        }
        cout<<"Finished reading Wasabi ion data with "<<fwasabiIonMap.size()<<" events filled in the time map"<<endl;
    }
    if (opt==1||opt==2){
        for (Long64_t ientry=0;ientry<fnentriesWasabiBeta;ientry++){
            fWasabiBeta->GetEntry(ientry);
            for (unsigned int i=0;i<fWasabiBeta->nhits;i++){
                wasabiHit* betaHit=new wasabiHit;
                betaHit->ts=fWasabiBeta->timestamp;
    #ifdef POS_CORR_CENTER_OF_GRAVITY
                    betaHit->x=fWasabiBeta->xcog[i];
                    betaHit->y=fWasabiBeta->ycog[i];
    #else
                    betaHit->x=fWasabiBeta->xmax[i];
                    betaHit->y=fWasabiBeta->ymax[i];
    #endif
                    betaHit->z=fWasabiBeta->z[i];
                    betaHit->nx=fWasabiBeta->nx[i];
                    betaHit->ny=fWasabiBeta->ny[i];
                    betaHit->nz=fWasabiBeta->nz;
                    betaHit->zmax=fWasabiBeta->zmax;
                    betaHit->ex=fWasabiBeta->ex[i];
                    betaHit->ey=fWasabiBeta->ey[i];
                    betaHit->blxmax=fWasabiBeta->blxmax[i];
                    betaHit->blymax=fWasabiBeta->blymax[i];
                    fwasabiBetaMap.insert(make_pair(fWasabiBeta->timestamp,betaHit));
            }
        }
        cout<<"Finished reading Wasabi beta data with "<<fwasabiBetaMap.size()<<" events filled in the time map"<<endl;
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Merger::ReadBigrips()
{
    for (unsigned int jentry = 0;jentry < (unsigned int) fnentriesBigrips;jentry++){
        ftrBigrips->GetEvent(jentry);
        fbigripsMap.insert(make_pair(fbigrips->ts,jentry));
    }
    cout<<"Finished reading bigrips ts table with "<<fbigripsMap.size()<<" rows"<<endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Merger::ReadBRIKEN(unsigned int startN, unsigned int stopN,unsigned int startG, unsigned int stopG,unsigned int startA, unsigned int stopA)
{
    //! read neutron
    unsigned int sstartN,sstopN;
    if (startN==0) sstartN = 0; else sstartN = startN;
    if (stopN==0) sstopN = (unsigned int) fnentriesNeutron; else sstopN = stopN;

    //! stuff for calculating cycle of time stamp reset
    unsigned long long prevts=0;
    unsigned short currcycle=0;
    for (unsigned int jentry = sstartN;jentry < sstopN;jentry++){
        ftrNeutron->GetEvent(jentry);
        if (fneutron->GetTimestamp()<prevts) currcycle++;// if time jump detected -> increase cycle by 1
        prevts=fneutron->GetTimestamp();
        BELENHit* hit=new BELENHit;
        fneutron->Copy(*hit);
        hit->SetHitsAdded(currcycle);
        fhe3Map.insert(make_pair(fneutron->GetTimestamp(),hit));
        //if (fneutron->GetEnergy()<fmaxneue&&fneutron->GetEnergy()>fminneue) fhe3Map.insert(make_pair(fneutron->GetTimestamp(),hit));
    }

    cout<<"Finished reading neutron  ts table with "<<fhe3Map.size()<<" rows"<<endl;
    //! read gamma
    unsigned int sstartG,sstopG;
    if (startG==0) sstartG = 0; else sstartG = startG;
    if (stopG==0) sstopG = (unsigned int) fnentriesGamma; else sstopG = stopG;

    for (unsigned int jentry = sstartG;jentry < sstopG;jentry++){
        ftrGamma->GetEvent(jentry);
        fcloverMap.insert(make_pair(fclover->GetTimestamp(),jentry));
    }
    cout<<"Finished reading gamma  ts table with "<<fcloverMap.size()<<" rows"<<endl;

    //! read anc
    unsigned int sstartA,sstopA;
    if (startA==0) sstartA = 0; else sstartA = startA;
    if (stopA==0) sstopA = (unsigned int) fnentriesAnc; else sstopA = stopA;

    prevts=0;
    currcycle=0;
    Long64_t prevtsdtpulser=0;
    Int_t ncountsDTpulser = 0;

    Long64_t firsttspulser1=0;

    for (unsigned int jentry = sstartA;jentry < sstopA;jentry++){
        ftrAnc->GetEvent(jentry);


        if (fanc->GetTimestamp()<prevts) currcycle++;// if time jump detected -> increase cycle by 1
        prevts=fanc->GetTimestamp();
        if (firsttspulser1==0) firsttspulser1=prevts;
        //fancMap.insert(make_pair(fanc->GetTimestamp(),jentry));
        if (fanc->GetMyPrecious()==1&&fanc->GetID()==1) {
            fF11RMap.insert(make_pair(fanc->GetTimestamp(),jentry));
        }
        if (fanc->GetMyPrecious()==1&&fanc->GetID()==2) fF11LMap.insert(make_pair(fanc->GetTimestamp(),jentry));

        if (fanc->GetMyPrecious()==1&&(fanc->GetID()==1||fanc->GetID()==2)) fF11LRMap.insert(make_pair(fanc->GetTimestamp(),jentry));

        if (fanc->GetMyPrecious()==2&&fanc->GetID()==1) fVetoTopMap.insert(make_pair(fanc->GetTimestamp(),jentry));
        if (fanc->GetMyPrecious()==2&&fanc->GetID()==2) fVetoBotMap.insert(make_pair(fanc->GetTimestamp(),jentry));

        if (fanc->GetMyPrecious()==3&&fanc->GetID()==1) fdETopMap.insert(make_pair(fanc->GetTimestamp(),jentry));
        if (fanc->GetMyPrecious()==3&&fanc->GetID()==2) fdEBotMap.insert(make_pair(fanc->GetTimestamp(),jentry));

        if (fanc->GetMyPrecious()==4) fVetoDownMap.insert(make_pair(fanc->GetTimestamp(),jentry));

        if (fanc->GetMyPrecious()==5) {
            if (ftsbeginpulser==0&&ncountsDTpulser==20) ftsbeginpulser=fanc->GetTimestamp();
            ftsendpulser=fanc->GetTimestamp();

            BELENHit* hit=new BELENHit;
            fanc->Copy(*hit);
            hit->SetID(141);
            hit->SetHitsAdded(currcycle);
            fdtpulserMap.insert(make_pair(fanc->GetTimestamp(),hit));

            if (ftsbeginpulser>0 && prevtsdtpulser!=0) {
                fh1dtpulser->Fill(fanc->GetEnergy());//for dead time calculation
                //fh1->Fill(ftsendpulser-prevtsdtpulser);
            }
            prevtsdtpulser=ftsendpulser;
            ncountsDTpulser++;
        }
    }
    if (ncountsDTpulser==0) {
        ftsbeginpulser=firsttspulser1+1.80000e+11;
        ftsendpulser=prevts;
    }


    //! read yso ion
    for (unsigned int jentry = 0;jentry < fnentriesYSOIon;jentry++){
        ftrYSOIon->GetEvent(jentry);
        YSOData* hit=new YSOData;
        copyYSOHit(fysoion,hit);
        fYSOIonMap.insert(make_pair(hit->T,hit));
    }

    //! read yso beta
    for (unsigned int jentry = 0;jentry < fnentriesYSOBeta;jentry++){
        ftrYSOBeta->GetEvent(jentry);
        YSOData* hit=new YSOData;
        copyYSOHit(fysobeta,hit);
        fYSOBetaMap.insert(make_pair(hit->T,hit));
    }

    //cout<<"Finished reading anc  ts table with "<<fancMap.size()<<" rows"<<endl;
    cout<<"Finished reading F11R  ts table with "<<fF11RMap.size()<<" rows"<<endl;
    cout<<"Finished reading F11L  ts table with "<<fF11LMap.size()<<" rows"<<endl;
    cout<<"Finished reading Veto Top ts table with "<<fVetoTopMap.size()<<" rows"<<endl;
    cout<<"Finished reading Veto Bottom ts table with "<<fVetoBotMap.size()<<" rows"<<endl;
    cout<<"Finished reading dE Top ts table with "<<fdETopMap.size()<<" rows"<<endl;
    cout<<"Finished reading dE Bottom ts table with "<<fdEBotMap.size()<<" rows"<<endl;
    cout<<"Finished reading Veto Down ts table with "<<fVetoDownMap.size()<<" rows"<<endl;
    cout<<"Finished reading YSO ion ts table with "<<fYSOIonMap.size()<<" rows"<<endl;
    cout<<"Finished reading YSO beta ts table with "<<fYSOBetaMap.size()<<" rows"<<endl;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Merger::BookIonBetaTree()
{
    ftreeallRI=new TTree("tree","tree");
    ftreeallRI->Branch("ion",fdataOutWasabiIon);
    ftreeallRI->Branch("beta",fdataOutWasabiBeta);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Merger::BookImplantTree()
{
    for (Int_t i=0;i<nri;i++){
        ftreeimplantRI[i]->Branch("ion",fdataOutWasabiIon);
        ftreeimplantRI[i]->Branch("implant",&implantdata,"yso_t/D:yso_e/D:yso_x/D:yso_y/D:zet/D:aoq/D:beta/D:F11L_T/D:F11L_E/D:F11R_T/D:F11R_E/D:F7_T/D:veto_T/D:veto_E/D");
        ftreeimplantRI[i]->Branch("gc_hit",&implantdata.gc_hit,"gc_hit/I");
        ftreeimplantRI[i]->Branch("gc_E",implantdata.gc_E,"gc_E[gc_hit]/D");
        ftreeimplantRI[i]->Branch("gc_T",implantdata.gc_T,"gc_T[gc_hit]/D");
        ftreeimplantRI[i]->Branch("gc_Tslew",implantdata.gc_Tslew,"gc_Tslew[gc_hit]/D");
        ftreeimplantRI[i]->Branch("gc_ch",implantdata.gc_ch,"gc_ch[gc_hit]/I");

        ftreeimplantRI[i]->Branch("gc1_hit",&implantdata.gc1_hit,"gc1_hit/I");
        ftreeimplantRI[i]->Branch("gc1_E",implantdata.gc1_E,"gc1_E[gc1_hit]/D");
        ftreeimplantRI[i]->Branch("gc1_T",implantdata.gc1_T,"gc1_T[gc1_hit]/D");
        ftreeimplantRI[i]->Branch("gc1_Tslew",implantdata.gc1_Tslew,"gc1_Tslew[gc1_hit]/D");
        ftreeimplantRI[i]->Branch("gc1_ch",implantdata.gc1_ch,"gc1_ch[gc1_hit]/I");

        ftreeimplantRI[i]->Branch("gc2_hit",&implantdata.gc2_hit,"gc2_hit/I");
        ftreeimplantRI[i]->Branch("gc2_E",implantdata.gc2_E,"gc2_E[gc2_hit]/D");
        ftreeimplantRI[i]->Branch("gc2_T",implantdata.gc2_T,"gc2_T[gc2_hit]/D");
        ftreeimplantRI[i]->Branch("gc2_Tslew",implantdata.gc2_Tslew,"gc2_Tslew[gc2_hit]/D");
        ftreeimplantRI[i]->Branch("gc2_ch",implantdata.gc2_ch,"gc2_ch[gc2_hit]/I");


        ftreeimplantRI[i]->Branch("ab1_hit",&implantdata.ab1_hit,"ab1_hit/I");
        ftreeimplantRI[i]->Branch("ab1_E",implantdata.ab1_E,"ab1_E[ab1_hit]/D");
        ftreeimplantRI[i]->Branch("ab1_T",implantdata.ab1_T,"ab1_T[ab1_hit]/D");
        ftreeimplantRI[i]->Branch("ab1_Tslew",implantdata.ab1_Tslew,"ab1_Tslew[ab1_hit]/D");
        ftreeimplantRI[i]->Branch("ab1_ch",implantdata.ab1_ch,"ab1_ch[ab1_hit]/I");
        ftreeimplantRI[i]->Branch("ab1_mult",implantdata.ab1_mult,"ab1_mult[ab1_hit]/S");

        ftreeimplantRI[i]->Branch("ab2_hit",&implantdata.ab2_hit,"ab2_hit/I");
        ftreeimplantRI[i]->Branch("ab2_E",implantdata.ab2_E,"ab2_E[ab2_hit]/D");
        ftreeimplantRI[i]->Branch("ab2_T",implantdata.ab2_T,"ab2_T[ab2_hit]/D");
        ftreeimplantRI[i]->Branch("ab2_Tslew",implantdata.ab2_Tslew,"ab2_Tslew[ab2_hit]/D");
        ftreeimplantRI[i]->Branch("ab2_ch",implantdata.ab2_ch,"ab2_ch[ab2_hit]/I");
        ftreeimplantRI[i]->Branch("ab2_mult",implantdata.ab2_mult,"ab2_mult[ab2_hit]/S");

        ftreeimplantRI[i]->Branch("neu_hit",&implantdata.neu_hit,"neu_hit/I");
        ftreeimplantRI[i]->Branch("neu_E",implantdata.neu_E,"neu_E[neu_hit]/D");
        ftreeimplantRI[i]->Branch("neu_T",implantdata.neu_T,"neu_T[neu_hit]/D");
        ftreeimplantRI[i]->Branch("neu_ch",implantdata.neu_ch,"neu_ch[neu_hit]/I");

    }
    ftreeimplantAll=new TTree("treeimp","treeimp");
    ftreeimplantAll->Branch("ion",fdataOutWasabiIon);
    ftreeimplantAll->Branch("implant",&implantdata,"yso_t/D:yso_e/D:yso_x/D:yso_y/D:zet/D:aoq/D:beta/D:F11L_T/D:F11L_E/D:F11R_T/D:F11R_E/D:F7_T/D:veto_T/D:veto_E/D");
    ftreeimplantAll->Branch("gc_hit",&implantdata.gc_hit,"gc_hit/I");
    ftreeimplantAll->Branch("gc_E",implantdata.gc_E,"gc_E[gc_hit]/D");
    ftreeimplantAll->Branch("gc_T",implantdata.gc_T,"gc_T[gc_hit]/D");
    ftreeimplantAll->Branch("gc_Tslew",implantdata.gc_Tslew,"gc_Tslew[gc_hit]/D");
    ftreeimplantAll->Branch("gc_ch",implantdata.gc_ch,"gc_ch[gc_hit]/I");

    ftreeimplantAll->Branch("gc1_hit",&implantdata.gc1_hit,"gc1_hit/I");
    ftreeimplantAll->Branch("gc1_E",implantdata.gc1_E,"gc1_E[gc1_hit]/D");
    ftreeimplantAll->Branch("gc1_T",implantdata.gc1_T,"gc1_T[gc1_hit]/D");
    ftreeimplantAll->Branch("gc1_Tslew",implantdata.gc1_Tslew,"gc1_Tslew[gc1_hit]/D");
    ftreeimplantAll->Branch("gc1_ch",implantdata.gc1_ch,"gc1_ch[gc1_hit]/I");

    ftreeimplantAll->Branch("gc2_hit",&implantdata.gc2_hit,"gc2_hit/I");
    ftreeimplantAll->Branch("gc2_E",implantdata.gc2_E,"gc2_E[gc2_hit]/D");
    ftreeimplantAll->Branch("gc2_T",implantdata.gc2_T,"gc2_T[gc2_hit]/D");
    ftreeimplantAll->Branch("gc2_Tslew",implantdata.gc2_Tslew,"gc2_Tslew[gc2_hit]/D");
    ftreeimplantAll->Branch("gc2_ch",implantdata.gc2_ch,"gc2_ch[gc2_hit]/I");


    ftreeimplantAll->Branch("ab1_hit",&implantdata.ab1_hit,"ab1_hit/I");
    ftreeimplantAll->Branch("ab1_E",implantdata.ab1_E,"ab1_E[ab1_hit]/D");
    ftreeimplantAll->Branch("ab1_T",implantdata.ab1_T,"ab1_T[ab1_hit]/D");
    ftreeimplantAll->Branch("ab1_Tslew",implantdata.ab1_Tslew,"ab1_Tslew[ab1_hit]/D");
    ftreeimplantAll->Branch("ab1_ch",implantdata.ab1_ch,"ab1_ch[ab1_hit]/I");
    ftreeimplantAll->Branch("ab1_mult",implantdata.ab1_mult,"ab1_mult[ab1_hit]/S");

    ftreeimplantAll->Branch("ab2_hit",&implantdata.ab2_hit,"ab2_hit/I");
    ftreeimplantAll->Branch("ab2_E",implantdata.ab2_E,"ab2_E[ab2_hit]/D");
    ftreeimplantAll->Branch("ab2_T",implantdata.ab2_T,"ab2_T[ab2_hit]/D");
    ftreeimplantAll->Branch("ab2_Tslew",implantdata.ab2_Tslew,"ab2_Tslew[ab2_hit]/D");
    ftreeimplantAll->Branch("ab2_ch",implantdata.ab2_ch,"ab2_ch[ab2_hit]/I");
    ftreeimplantAll->Branch("ab2_mult",implantdata.ab2_mult,"ab2_mult[ab2_hit]/S");

    ftreeimplantAll->Branch("neu_hit",&implantdata.neu_hit,"neu_hit/I");
    ftreeimplantAll->Branch("neu_E",implantdata.neu_E,"neu_E[neu_hit]/D");
    ftreeimplantAll->Branch("neu_T",implantdata.neu_T,"neu_T[neu_hit]/D");
    ftreeimplantAll->Branch("neu_ch",implantdata.neu_ch,"neu_ch[neu_hit]/I");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void::Merger::ResetImplantData()
{
    implantdata.yso_t=-9999;implantdata.yso_e=-9999;implantdata.yso_x=-9999;implantdata.yso_y=-9999;
    implantdata.zet=-9999.;implantdata.aoq=-9999.;implantdata.beta=-9999;

    implantdata.F11L_E=-9999.;implantdata.F11L_T=-9999.;implantdata.F11R_E=-9999.;implantdata.F11R_T=-9999.;
    implantdata.F7_T=-9999.;implantdata.veto_E=-9999.;implantdata.veto_T=-9999.;
    implantdata.gc_hit=0;
    implantdata.gc1_hit=0;
    implantdata.gc2_hit=0;
    implantdata.ab1_hit=0;
    implantdata.ab2_hit=0;

    for (int i=0;i<kMaxGamma;i++){
        implantdata.gc_E[i]=-9999.;
        implantdata.gc_T[i]=-9999.;
        implantdata.gc_Tslew[i]=-9999.;
        implantdata.gc_ch[i]=-9999;

        implantdata.gc1_E[i]=-9999.;
        implantdata.gc1_T[i]=-9999.;
        implantdata.gc1_Tslew[i]=-9999.;
        implantdata.gc1_ch[i]=-9999;
        implantdata.gc2_E[i]=-9999.;
        implantdata.gc2_T[i]=-9999.;
        implantdata.gc2_Tslew[i]=-9999.;
        implantdata.gc2_ch[i]=-9999;

        implantdata.ab1_E[i]=-9999.;
        implantdata.ab1_ch[i]=-9999.;
        implantdata.ab1_T[i]=-9999;
        implantdata.ab1_Tslew[i]=-9999.;
        implantdata.ab2_E[i]=-9999.;
        implantdata.ab2_ch[i]=-9999.;
        implantdata.ab2_T[i]=-9999;
        implantdata.ab2_Tslew[i]=-9999.;

        implantdata.ab1_mult[i]=0;
        implantdata.ab2_mult[i]=0;
    }

    implantdata.neu_hit=0;
    for (int i=0;i<kMaxNeutron;i++){
        implantdata.neu_E[i]=-9999.;
        implantdata.neu_T[i]=-9999.;
        implantdata.neu_ch[i]=-9999;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Merger::DoAddback()
{
    Bool_t ab_started=false;
    Long64_t ab_beg=0;
    Long64_t ab_end=0;
    Long64_t ab_window=1000;
    Long64_t ts_prev=0;

    Int_t k=0;
    Int_t ch_beg=0;
    Double_t esum=0;
    Double_t e_chbeg=0;

    gammaab* abdata=new gammaab();
    gammaab* abdata2=new gammaab();

    for (fcloverMap_it=fcloverMap.begin();fcloverMap_it!=fcloverMap.end();fcloverMap_it++){
        if (k%100000==0) cout<<"Addback clover D4 "<<k<<"/"<<fcloverMap.size()<<"\r"<<flush;
        Long64_t ts = (Long64_t) fcloverMap_it->first;
        unsigned int entry = fcloverMap_it->second;
        ftrGamma->GetEvent(entry);
        if (fclover->GetClover()==1){
            Double_t ecal=fclover->GetEnergy();
            /*
            //! gamma recalibration goes here
            //convert back to adc
            ecal=(ecal-fcoffsetold[fclover->GetID()-1])/fcgainold[fclover->GetID()-1];
            //apply new calibration
            if (ecal<fsep[fclover->GetID()-1]){//low energy calibration
                ecal=flow_offset[fclover->GetID()-1]+flow_gain[fclover->GetID()-1]*ecal+flow_se[fclover->GetID()-1]*ecal*ecal;
            }else{//high energy calibration
                ecal=fhigh_offset[fclover->GetID()-1]+fhigh_gain[fclover->GetID()-1]*ecal+fhigh_se[fclover->GetID()-1]*ecal*ecal;
            }
            */

            if (!ab_started) {
                ab_beg=ts;
                ab_started=true;
                e_chbeg=ecal;
            }

            if (ch_beg==0) ch_beg=fclover->GetCloverLeaf();

            if (ts-ab_beg>ab_window&&ab_started){// end of event, next event start
                ab_end=ts_prev;
                //--------event operation here
                abdata->ab_ch=ch_beg;
                abdata->ab_E=esum;
                abdata->ab_T=ab_beg;
                //perform slew correction here
                abdata->ab_Tslew=0;
                if (isslewcorr){
                    abdata->ab_Tslew=(d[ch_beg-1]+(a[ch_beg-1]-d[ch_beg-1])/(1+pow(e_chbeg/c[ch_beg-1],b[ch_beg-1])));
                }
                gammaab* abobj=new gammaab();
                CopyAddbackData(abdata,abobj);
                faddbackclover1Map.insert(make_pair(abdata->ab_T,abobj));

                e_chbeg=ecal;
                ch_beg=fclover->GetCloverLeaf();
                ab_beg=ts;// start new event
                esum=0;
                memset(abdata->ab_mult,0,sizeof(abdata->ab_mult));
            }
            //" colleting hits to an event here
            esum+=ecal;
            if (fclover->GetCloverLeaf()==1) abdata->ab_mult[0]++;
            else if (fclover->GetCloverLeaf()==2) abdata->ab_mult[1]++;
            else if (fclover->GetCloverLeaf()==3) abdata->ab_mult[2]++;
            else if (fclover->GetCloverLeaf()==4) abdata->ab_mult[3]++;

            ts_prev=ts;
        }
        k++;
    }
    if (abdata->ab_mult[0]+abdata->ab_mult[1]+abdata->ab_mult[2]+abdata->ab_mult[3]>0){
        abdata->ab_ch=ch_beg;
        abdata->ab_E=esum;
        abdata->ab_T=ab_beg;
        //perform slew correction here
        abdata->ab_Tslew=0;
        if (isslewcorr){
            abdata->ab_Tslew=(d[ch_beg-1]+(a[ch_beg-1]-d[ch_beg-1])/(1+pow(e_chbeg/c[ch_beg-1],b[ch_beg-1])));
        }
        gammaab* abobj=new gammaab();
        CopyAddbackData(abdata,abobj);
        faddbackclover1Map.insert(make_pair(abdata->ab_T,abobj));
    }

    //! add back for 2nd clover
    ab_started=false;
    ab_beg=0;
    ab_end=0;
    ab_window=1000;
    ts_prev=0;

    k=0;
    ch_beg=0;
    esum=0;
    e_chbeg=0;

    for (fcloverMap_it=fcloverMap.begin();fcloverMap_it!=fcloverMap.end();fcloverMap_it++){
        if (k%100000==0) cout<<"Addback clover G7 "<<k<<"/"<<fcloverMap.size()<<"\r"<<flush;
        Long64_t ts = (Long64_t) fcloverMap_it->first;
        unsigned int entry = fcloverMap_it->second;
        ftrGamma->GetEvent(entry);
        if (fclover->GetClover()==2){
            Double_t ecal=fclover->GetEnergy();
            /*
            //! gamma recalibration goes here
            //convert back to adc
            ecal=(ecal-fcoffsetold[fclover->GetID()-1])/fcgainold[fclover->GetID()-1];
            //apply new calibration
            if (ecal<fsep[fclover->GetID()-1]){//low energy calibration
                ecal=flow_offset[fclover->GetID()-1]+flow_gain[fclover->GetID()-1]*ecal+flow_se[fclover->GetID()-1]*ecal*ecal;
            }else{//high energy calibration
                ecal=fhigh_offset[fclover->GetID()-1]+fhigh_gain[fclover->GetID()-1]*ecal+fhigh_se[fclover->GetID()-1]*ecal*ecal;
            }
            */

            if (!ab_started) {
                ab_beg=ts;
                ab_started=true;
                e_chbeg=ecal;
            }

            if (ch_beg==0) ch_beg=fclover->GetCloverLeaf();

            if (ts-ab_beg>ab_window&&ab_started){// end of event, next event start
                ab_end=ts_prev;
                //--------event operation here

                abdata2->ab_ch=ch_beg;
                abdata2->ab_E=esum;
                abdata2->ab_T=ab_beg;
                //perform slew correction here
                abdata2->ab_Tslew=0;
                if (isslewcorr){
                    abdata2->ab_Tslew=(d[ch_beg+3]+(a[ch_beg+3]-d[ch_beg+3])/(1+pow(e_chbeg/c[ch_beg+3],b[ch_beg+3])));
                }
                gammaab* abobj=new gammaab();
                CopyAddbackData(abdata2,abobj);
                faddbackclover2Map.insert(make_pair(abdata2->ab_T,abobj));

                e_chbeg=ecal;
                ch_beg=fclover->GetCloverLeaf();
                ab_beg=ts;// start new event
                esum=0;
                memset(abdata2->ab_mult,0,sizeof(abdata2->ab_mult));
            }
            //" colleting hits to an event here
            esum+=ecal;
            if (fclover->GetCloverLeaf()==1) abdata2->ab_mult[0]++;
            else if (fclover->GetCloverLeaf()==2) abdata2->ab_mult[1]++;
            else if (fclover->GetCloverLeaf()==3) abdata2->ab_mult[2]++;
            else if (fclover->GetCloverLeaf()==4) abdata2->ab_mult[3]++;
            ts_prev=ts;
        }
        k++;
    }
    if (abdata2->ab_mult[0]+abdata2->ab_mult[1]+abdata2->ab_mult[2]+abdata2->ab_mult[3]>0){
        abdata2->ab_ch=ch_beg;
        abdata2->ab_E=esum;
        abdata2->ab_T=ab_beg;
        //perform slew correction here
        abdata2->ab_Tslew=0;
        if (isslewcorr){
            abdata2->ab_Tslew=(d[ch_beg+3]+(a[ch_beg+3]-d[ch_beg+3])/(1+pow(e_chbeg/c[ch_beg+3],b[ch_beg+3])));
        }
        gammaab* abobj=new gammaab();
        CopyAddbackData(abdata2,abobj);
        faddbackclover2Map.insert(make_pair(abdata2->ab_T,abobj));
    }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Merger::ReadSlewCorr()
{
    std::ifstream ifs("slewcorrparms.txt");
    for (Int_t i=0;i<8;i++){
        Int_t temp;
        ifs>>temp>>a[i]>>b[i]>>c[i]>>d[i];
    }
    isslewcorr=true;
    ifs.close();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Merger::MergeIonBeta() //just for testing
{
    //!******************************Implantation****************
    //! Time correlation with implantation (Build Decay Curve)
    Long64_t ktotal=fwasabiBetaMap.size();
    Long64_t k=0;
    Long64_t ntotalionbetacorr=0;
    for (fwasabiBetaMap_it=fwasabiBetaMap.begin();fwasabiBetaMap_it!=fwasabiBetaMap.end();fwasabiBetaMap_it++){
        if (k%10000==0) cout<<k<<"/"<<ktotal<<"\tncorr="<<ntotalionbetacorr<<"\r"<<flush;
        Long64_t ts,ts1,ts2,corrts,correntry,check_time,ncorr;

        ts=fwasabiBetaMap_it->first;
        wasabiHit* hit=(wasabiHit*) fwasabiBetaMap_it->second;

        //! correlate with implant
        ts1 = ts - fTW_IonBetalow;
        ts2 = ts + fTW_IonBetaup;
        corrts = 0;
        correntry = 0;
        check_time = 0;
        ncorr = 0;

        fwasabiIonMap_it = fwasabiIonMap.lower_bound(ts1);
        while(fwasabiIonMap_it!=fwasabiIonMap.end()&&fwasabiIonMap_it->first<ts2){
            corrts=fwasabiIonMap_it->first;
            wasabiHit* corrhit=(wasabiHit*) fwasabiIonMap_it->second;
            if ((hit->z==corrhit->z)&&(pow((hit->x-corrhit->x),2)<fsquareSpatialWindowX)&&(pow((hit->y-corrhit->y),2)<fsquareSpatialWindowY)){
                //! fill data here!
                h1dProbe[0]->Fill(ts-corrts);
                hit->Copy(*fdataOutWasabiBeta);
                corrhit->Copy(*fdataOutWasabiIon);
                if (ftreeallRI!=0) ftreeallRI->Fill();
                check_time=corrts;
                ncorr++;
                ntotalionbetacorr++;
            }
            fwasabiIonMap_it++;
        }
        k++;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Merger::MergeImplant()
{
    Long64_t ktotal=fwasabiIonMap.size();
    Long64_t k=0;
    Long64_t ncorrwbigrips=0;
    Long64_t ncorrwclover=0;
    Long64_t ncorrwyso=0;

    ResetImplantData();
    for (fwasabiIonMap_it=fwasabiIonMap.begin();fwasabiIonMap_it!=fwasabiIonMap.end();fwasabiIonMap_it++){
        if (k%10000==0) cout<<k<<"/"<<ktotal<<"\t ncorr w bigrips "<<ncorrwbigrips<<"\t ncorr w clover "<<ncorrwclover<<"\t ncorr w yso "<<ncorrwyso<<"\r"<<flush;
        Long64_t ts=fwasabiIonMap_it->first;

        wasabiHit* ionhit=(wasabiHit*) fwasabiIonMap_it->second;
        ImplantCorrelationVector* corrvector=new ImplantCorrelationVector;
        corrvector->correntrybrips=-9999;
        corrvector->correntryf11r=-9999;
        corrvector->correntryf11l=-9999;
        corrvector->correntrydEbot=-9999;
        corrvector->correntrydEtop=-9999;
        corrvector->correntryvetodown=-9999;

        corrvector->yso_t=-9999;
        corrvector->yso_e=-9999;
        corrvector->yso_x=-9999;
        corrvector->yso_y=-9999;

        corrvector->gammagc1_vector.clear();
        corrvector->gammagc2_vector.clear();
        corrvector->gammaab1_vector.clear();
        corrvector->gammaab2_vector.clear();

        //! Correlate imp with bigrips
        Long64_t ts1 = ts - fTW_IonPidlow;
        Long64_t ts2 = ts + fTW_IonPidup;
        Long64_t corrts = 0;
        Int_t ncorr=0;
        unsigned int correntry = 0;
        Long64_t check_time = 0;
        fbigripsMap_it = fbigripsMap.lower_bound(ts1);
        while(fbigripsMap_it!=fbigripsMap.end()&&fbigripsMap_it->first<ts2){
            corrts =  fbigripsMap_it->first;
            if (corrts!=check_time){
                check_time=corrts;
                corrvector->correntrybrips = (int) fbigripsMap_it->second;
                //! fill data here
                ncorr++;
                break;
            }
            fbigripsMap_it++;
        }
        if (ncorr>0) ncorrwbigrips++;

        //! Correlate with yso ion
        ts1 = ts - fTW_IonYSOionlow;
        ts2 = ts + fTW_IonYSOionup;
        corrts = 0;
        ncorr=0;
        correntry = 0;
        check_time = 0;
        fYSOIonMap_it = fYSOIonMap.lower_bound(ts1);

        while(fYSOIonMap_it!=fYSOIonMap.end()&&fYSOIonMap_it->first<ts2){
            corrts =  fYSOIonMap_it->first;
            if (corrts!=check_time){
                check_time=corrts;
                YSOData* hit = (YSOData*) fYSOIonMap_it->second;
                corrvector->yso_t=(Double_t)((Long64_t)hit->T-ts);
                corrvector->yso_e=hit->E;
                corrvector->yso_x=hit->x;
                corrvector->yso_y=hit->y;
                ncorr++;
                break;
            }
            fYSOIonMap_it++;
        }
        if (ncorr>0) ncorrwyso++;

        //! Correlate with gamma
        ts1 = ts - fTW_PIDGammalow;
        ts2 = ts + fTW_PIDGammaup;
        corrts = 0;
        ncorr=0;
        correntry = 0;
        check_time = 0;
        fcloverMap_it = fcloverMap.lower_bound(ts1);

        while(fcloverMap_it!=fcloverMap.end()&&fcloverMap_it->first<ts2){
            corrts =  fcloverMap_it->first;
            correntry = fcloverMap_it->second;
            if (corrts!=check_time){
                check_time=corrts;
                ftrGamma->GetEvent(correntry);
                gammahit* ghit=new gammahit;
                Double_t ecal=fclover->GetEnergy();
                /*
                //! gamma recalibration goes here
                //convert back to adc
                ecal=(ecal-fcoffsetold[fclover->GetID()-1])/fcgainold[fclover->GetID()-1];
                //apply new calibration
                if (ecal<fsep[fclover->GetID()-1]){//low energy calibration
                    ecal=flow_offset[fclover->GetID()-1]+flow_gain[fclover->GetID()-1]*ecal+flow_se[fclover->GetID()-1]*ecal*ecal;
                }else{//high energy calibration
                    ecal=fhigh_offset[fclover->GetID()-1]+fhigh_gain[fclover->GetID()-1]*ecal+fhigh_se[fclover->GetID()-1]*ecal*ecal;
                }
                */

                ghit->gc_ch=fclover->GetCloverLeaf();
                ghit->gc_E=ecal;
                ghit->gc_T=(Double_t)(corrts-ts);
                //! time slew correction goes here
                ghit->gc_Tslew=0;
                if (isslewcorr){
                    ghit->gc_Tslew=(d[ghit->gc_ch-1]+(a[ghit->gc_ch-1]-d[ghit->gc_ch-1])/(1+pow(ghit->gc_E/c[ghit->gc_ch-1],b[ghit->gc_ch-1])));
                }

                if (fclover->GetClover()==1)
                    corrvector->gammagc1_vector.push_back(ghit);
                else
                    corrvector->gammagc2_vector.push_back(ghit);
                ncorr++;
            }
            fcloverMap_it++;
        }
        if (ncorr>0) ncorrwclover++;

        //! Correlate with addback gamma data clover D4
        ts1 = ts - fTW_PIDGammalow;
        ts2 = ts + fTW_PIDGammaup;
        corrts = 0;
        ncorr=0;
        correntry = 0;
        check_time = 0;

        faddbackclover1Map_it = faddbackclover1Map.lower_bound(ts1);

        while(faddbackclover1Map_it!=faddbackclover1Map.end()&&faddbackclover1Map_it->first<ts2){
            corrts =  faddbackclover1Map_it->first;
            gammaab* ab1 = faddbackclover1Map_it->second;
            if (corrts!=check_time){
                check_time=corrts;

                ab1->ab_T=(Double_t)(corrts-ts);
                ab1->ab_Tslew=0;
                //! time slew correction goes here

                corrvector->gammaab1_vector.emplace_back(ab1);
                ncorr++;
            }
            faddbackclover1Map_it++;
        }

        //! Correlate with addback gamma data clover G7
        ts1 = ts - fTW_PIDGammalow;
        ts2 = ts + fTW_PIDGammaup;
        corrts = 0;
        ncorr=0;
        correntry = 0;
        check_time = 0;

        faddbackclover2Map_it = faddbackclover2Map.lower_bound(ts1);
        while(faddbackclover2Map_it!=faddbackclover2Map.end()&&faddbackclover2Map_it->first<ts2){
            corrts =  faddbackclover2Map_it->first;
            gammaab* ab2 = faddbackclover2Map_it->second;
            if (corrts!=check_time){
                check_time=corrts;

                ab2->ab_T=(Double_t)(corrts-ts);
                ab2->ab_Tslew=0;
                //! time slew correction goes here

                corrvector->gammaab2_vector.emplace_back(ab2);
                ncorr++;
            }
            faddbackclover2Map_it++;
        }

        //! Correlate imp with f11r
        ts1 = ts - fTW_IonF11LRlow;
        ts2 = ts + fTW_IonF11LRup;
        corrts = 0;
        ncorr=0;
        check_time = 0;
        fF11MapR_it = fF11RMap.lower_bound(ts1);
        while(fF11MapR_it!=fF11RMap.end()&&fF11MapR_it->first<ts2){
            corrts =  fF11MapR_it->first;
            if (corrts!=check_time){
                check_time=corrts;
                corrvector->correntryf11r = (int) fF11MapR_it->second;
                ncorr++;
                break;
            }
            fF11MapR_it++;
        }


        //! Correlate imp with f11l
        ts1 = ts - fTW_IonF11LRlow;
        ts2 = ts + fTW_IonF11LRup;
        corrts = 0;
        ncorr=0;
        check_time = 0;
        fF11MapL_it = fF11LMap.lower_bound(ts1);
        while(fF11MapL_it!=fF11LMap.end()&&fF11MapL_it->first<ts2){
            corrts =  fF11MapL_it->first;
            if (corrts!=check_time){
                check_time=corrts;
                corrvector->correntryf11l = (int) fF11MapL_it->second;
                ncorr++;
                break;
            }
            fF11MapL_it++;
        }

        //! Correlate imp with vetodown
        ts1 = ts - fTW_IonDownVetolow;
        ts2 = ts + fTW_IonDownVetoup;
        corrts = 0;
        ncorr=0;
        check_time = 0;
        fVetoDownMap_it = fVetoDownMap.lower_bound(ts1);
        while(fVetoDownMap_it!=fVetoDownMap.end()&&fVetoDownMap_it->first<ts2){
            corrts =  fVetoDownMap_it->first;
            if (corrts!=check_time){
                check_time=corrts;
                corrvector->correntryvetodown = (int) fVetoDownMap_it->second;
                ncorr++;
                break;
            }
            fVetoDownMap_it++;
        }

        fwasabiImplantMap.insert(make_pair(ts,make_pair(corrvector,ionhit)));
        k++;
    }
    //! Write out implantation tree
    if (ftreeimplantAll!=0){
        for (fwasabiImplantMap_it=fwasabiImplantMap.begin();fwasabiImplantMap_it!=fwasabiImplantMap.end();fwasabiImplantMap_it++){
            wasabiHit* hit=fwasabiImplantMap_it->second.second;
            ImplantCorrelationVector* corrvectorimp=fwasabiImplantMap_it->second.first;
            //! fill wasabidata here
            hit->Copy(*fdataOutWasabiIon);

            implantdata.yso_e=corrvectorimp->yso_e;
            implantdata.yso_t=corrvectorimp->yso_t;
            implantdata.yso_x=corrvectorimp->yso_x;
            implantdata.yso_y=corrvectorimp->yso_y;

            //! fill bigrips data here
            if (corrvectorimp->correntrybrips>=0) {
                ftrBigrips->GetEvent(corrvectorimp->correntrybrips);
                implantdata.zet=fbigrips->zet;
                implantdata.aoq=fbigrips->aoq;
                implantdata.beta=fbigrips->beta;
                implantdata.F7_T=(Double_t)((Long64_t)fbigrips->ts-hit->ts);
            }else{
                implantdata.zet=-9999;
                implantdata.aoq=-9999;
                implantdata.beta=-9999;
                implantdata.F7_T=-9999;
            }

            //! fill anc data here
            if (corrvectorimp->correntryf11r>=0){
                ftrAnc->GetEvent(corrvectorimp->correntryf11r);
                implantdata.F11R_E=fanc->GetEnergy();
                implantdata.F11R_T=(Double_t)((Long64_t)fanc->GetTimestamp()-hit->ts);
            }else{
                implantdata.F11R_E=-9999;
                implantdata.F11R_T=-9999;
            }

            if (corrvectorimp->correntryf11l>=0){
                ftrAnc->GetEvent(corrvectorimp->correntryf11l);
                implantdata.F11L_E=fanc->GetEnergy();
                implantdata.F11L_T=(Double_t)((Long64_t)fanc->GetTimestamp()-hit->ts);
            }else{
                implantdata.F11L_E=-9999;
                implantdata.F11L_T=-9999;
            }

            if (corrvectorimp->correntryvetodown>=0){
                ftrAnc->GetEvent(corrvectorimp->correntryvetodown);
                implantdata.veto_E=fanc->GetEnergy();
                implantdata.veto_T=(Double_t)((Long64_t)fanc->GetTimestamp()-hit->ts);
            }else{
                implantdata.veto_E=-9999;
                implantdata.veto_T=-9999;
            }

            //! fill gamma data here
            std::vector<gammahit*>::iterator gammagc1_vector_it;
            Int_t nhit=0;
            for (gammagc1_vector_it=corrvectorimp->gammagc1_vector.begin();gammagc1_vector_it!=corrvectorimp->gammagc1_vector.end();gammagc1_vector_it++){
                gammahit* gchit=*gammagc1_vector_it;
                implantdata.gc1_ch[nhit]=gchit->gc_ch;
                implantdata.gc1_E[nhit]=gchit->gc_E;
                implantdata.gc1_T[nhit]=gchit->gc_T;
                implantdata.gc1_Tslew[nhit]=gchit->gc_Tslew;
                nhit++;
            }
            implantdata.gc1_hit=nhit;

            std::vector<gammahit*>::iterator gammagc2_vector_it;
            nhit=0;
            for (gammagc2_vector_it=corrvectorimp->gammagc2_vector.begin();gammagc2_vector_it!=corrvectorimp->gammagc2_vector.end();gammagc2_vector_it++){
                gammahit* gchit=*gammagc2_vector_it;
                implantdata.gc2_ch[nhit]=gchit->gc_ch;
                implantdata.gc2_E[nhit]=gchit->gc_E;
                implantdata.gc2_T[nhit]=gchit->gc_T;
                implantdata.gc2_Tslew[nhit]=gchit->gc_Tslew;
                nhit++;
            }
            implantdata.gc2_hit=nhit;

            //! fill gamma addback data here
            std::vector<gammaab*>::iterator gammaab1_vector_it;
            nhit=0;
            for (gammaab1_vector_it=corrvectorimp->gammaab1_vector.begin();gammaab1_vector_it!=corrvectorimp->gammaab1_vector.end();gammaab1_vector_it++){
                gammaab* abhit=*gammaab1_vector_it;
                implantdata.ab1_ch[nhit]=abhit->ab_ch;
                implantdata.ab1_E[nhit]=abhit->ab_E;
                implantdata.ab1_T[nhit]=abhit->ab_T;
                implantdata.ab1_Tslew[nhit]=abhit->ab_Tslew;
                implantdata.ab1_mult[nhit]=abhit->ab_mult[0]+abhit->ab_mult[1]+abhit->ab_mult[2]+abhit->ab_mult[3];
                nhit++;
            }
            implantdata.ab1_hit=nhit;

            std::vector<gammaab*>::iterator gammaab2_vector_it;
            nhit=0;
            for (gammaab2_vector_it=corrvectorimp->gammaab2_vector.begin();gammaab2_vector_it!=corrvectorimp->gammaab2_vector.end();gammaab2_vector_it++){
                gammaab* abhit=*gammaab2_vector_it;
                implantdata.ab2_ch[nhit]=abhit->ab_ch;
                implantdata.ab2_E[nhit]=abhit->ab_E;
                implantdata.ab2_T[nhit]=abhit->ab_T;
                implantdata.ab2_Tslew[nhit]=abhit->ab_Tslew;
                implantdata.ab2_mult[nhit]=abhit->ab_mult[0]+abhit->ab_mult[1]+abhit->ab_mult[2]+abhit->ab_mult[3];
                nhit++;
            }
            implantdata.ab2_hit=nhit;
            //! fill neutron data here
            ftreeimplantAll->Fill();
            //! separate tree and fill
            for (Int_t j=0;j<nri;j++){
                if (!enablepid2[j]) continue;
                if (cutg[j]->IsInside(implantdata.aoq,implantdata.zet)){
                    ftreeimplantRI[j]->Fill();
                }
            }

        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Merger::MakeF11NeutronVeto()
{

}
