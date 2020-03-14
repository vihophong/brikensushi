//only command line only

//X side
root -l outrootfiles/bi207_calib_00042to00045.root
.L calibspec.C
TH1F* h[100];Int_t i=0;
i=0;aida->Draw(Form("ex>>h%d(200,0,2000)",i),Form("z==3&&x==%d",i),"goff");h[i]=(TH1F*)gDirectory->Get(Form("h%d",i));calibspec(h[i],30,30,40,40);
i=1;fitcalibcurve4(1,317.5,372.5,667.5,722.5);
i=2;aida->Draw(Form("ex>>h%d(200,0,2000)",i),Form("z==3&&x==%d",i),"goff");h[i]=(TH1F*)gDirectory->Get(Form("h%d",i));calibspec(h[i],30,30,40,40);
i=3;aida->Draw(Form("ex>>h%d(200,0,2000)",i),Form("z==3&&x==%d",i),"goff");h[i]=(TH1F*)gDirectory->Get(Form("h%d",i));calibspec(h[i],30,30,40,40);
i=4;i=4;aida->Draw(Form("ex>>h%d(200,0,2000)",i),Form("z==3&&x==%d",i),"goff");h[i]=(TH1F*)gDirectory->Get(Form("h%d",i));calibspec(h[i],30,50,40,40);//some difference in begin-end of fit
i=5;aida->Draw(Form("ex>>h%d(200,0,2000)",i),Form("z==3&&x==%d",i),"goff");h[i]=(TH1F*)gDirectory->Get(Form("h%d",i));calibspec(h[i],30,30,40,40);
i=6;aida->Draw(Form("ex>>h%d(200,0,2000)",i),Form("z==3&&x==%d",i),"goff");h[i]=(TH1F*)gDirectory->Get(Form("h%d",i));calibspec(h[i],30,30,40,40);
i=7;fitcalibcurve3(7,604.361,1250.84,1345.92);
i=8;aida->Draw(Form("ex>>h%d(200,0,2000)",i),Form("z==3&&x==%d",i),"goff");h[i]=(TH1F*)gDirectory->Get(Form("h%d",i));calibspec(h[i],30,30,40,40);
i=9;aida->Draw(Form("ex>>h%d(200,0,2000)",i),Form("z==3&&x==%d",i),"goff");h[i]=(TH1F*)gDirectory->Get(Form("h%d",i));calibspec(h[i],30,30,40,40);
i=10;aida->Draw(Form("ex>>h%d(200,0,2000)",i),Form("z==3&&x==%d",i),"goff");h[i]=(TH1F*)gDirectory->Get(Form("h%d",i));calibspec(h[i],30,30,40,40);
i=11;aida->Draw(Form("ex>>h%d(200,0,2000)",i),Form("z==3&&x==%d",i),"goff");h[i]=(TH1F*)gDirectory->Get(Form("h%d",i));calibspec(h[i],30,30,40,40);
i=12;aida->Draw(Form("ex>>h%d(200,0,2000)",i),Form("z==3&&x==%d",i),"goff");h[i]=(TH1F*)gDirectory->Get(Form("h%d",i));calibspec(h[i],30,30,40,40);
i=13;aida->Draw(Form("ex>>h%d(200,0,2000)",i),Form("z==3&&x==%d",i),"goff");h[i]=(TH1F*)gDirectory->Get(Form("h%d",i));calibspec(h[i],30,30,40,40);
i=14;aida->Draw(Form("ex>>h%d(200,0,2000)",i),Form("z==3&&x==%d",i),"goff");h[i]=(TH1F*)gDirectory->Get(Form("h%d",i));calibspec(h[i],30,30,40,40);
i=15;fitcalibcurve3(15,563.312,1179.11,1275.19);

//Y side
root -l outrootfiles/bi207_calib_00042to00045.root
.L calibspec.C
TH1F* h[100];Int_t i;
i=0;fitcalibcurve3(0,605,1285,1405);
i=1;fitcalibcurve3(1,595,1275,1385);
i=2;aida->Draw(Form("ey>>h%d(200,0,2000)",i),Form("z==3&&y==%d",i),"goff");h[i]=(TH1F*)gDirectory->Get(Form("h%d",i));calibspec(h[i],30,30,40,40);
i=3;aida->Draw(Form("ey>>h%d(200,0,2000)",i),Form("z==3&&y==%d",i),"goff");h[i]=(TH1F*)gDirectory->Get(Form("h%d",i));calibspec(h[i],30,30,40,40);
i=4;aida->Draw(Form("ey>>h%d(200,0,2000)",i),Form("z==3&&y==%d",i),"goff");h[i]=(TH1F*)gDirectory->Get(Form("h%d",i));calibspec(h[i],30,30,40,40);
i=5;i=5;aida->Draw(Form("ey>>h%d(200,0,2000)",i),Form("z==3&&y==%d",i),"goff");h[i]=(TH1F*)gDirectory->Get(Form("h%d",i));calibspec(h[i],50,50,50,50);
i=6;aida->Draw(Form("ey>>h%d(200,0,2000)",i),Form("z==3&&y==%d",i),"goff");h[i]=(TH1F*)gDirectory->Get(Form("h%d",i));calibspec(h[i],30,30,40,40);
i=7;aida->Draw(Form("ey>>h%d(200,0,2000)",i),Form("z==3&&y==%d",i),"goff");h[i]=(TH1F*)gDirectory->Get(Form("h%d",i));calibspec(h[i],30,30,40,40);
i=8;aida->Draw(Form("ey>>h%d(200,0,2000)",i),Form("z==3&&y==%d",i),"goff");h[i]=(TH1F*)gDirectory->Get(Form("h%d",i));calibspec(h[i],30,30,40,40);
i=9;i=9;aida->Draw(Form("ey>>h%d(200,0,2000)",i),Form("z==3&&y==%d",i),"goff");h[i]=(TH1F*)gDirectory->Get(Form("h%d",i));calibspec(h[i],40,30,50,40);
i=10;fitcalibcurve3(10,612.831,1275.32,1378.58);
i=11;aida->Draw(Form("ey>>h%d(200,0,2000)",i),Form("z==3&&y==%d",i),"goff");h[i]=(TH1F*)gDirectory->Get(Form("h%d",i));calibspec(h[i],30,30,40,40);
i=12;aida->Draw(Form("ey>>h%d(200,0,2000)",i),Form("z==3&&y==%d",i),"goff");h[i]=(TH1F*)gDirectory->Get(Form("h%d",i));calibspec(h[i],30,30,40,40);
i=13;fitcalibcurve3(13,614.576,1278.97,1379.77);
i=14;aida->Draw(Form("ey>>h%d(200,0,2000)",i),Form("z==3&&y==%d",i),"goff");h[i]=(TH1F*)gDirectory->Get(Form("h%d",i));calibspec(h[i],30,30,40,40);
i=15;aida->Draw(Form("ey>>h%d(200,0,2000)",i),Form("z==3&&y==%d",i),"goff");h[i]=(TH1F*)gDirectory->Get(Form("h%d",i));calibspec(h[i],30,30,40,40);






