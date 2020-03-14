#ifndef BELEN_H
#define BELEN_H
#include <iostream>
#include <vector>
#include <cstdlib>
#include <math.h>

#include "TObject.h"
#include "TVector3.h"
#include "TMath.h"
#include "BELENdefs.hh"

using namespace std;


class BELENHit : public TObject{
public:
    BELENHit(){
        Clear();
    }
    BELENHit(Double_t posx, Double_t posy, Double_t posz, short daqid, short id, short ring,short type, unsigned long long ts, int adc, int en, unsigned short hitsadded)
    {
        fdaqid = daqid;
        fid = id;
        fring = ring;
        ftype = type;
        fts = ts;
        fadc = adc;
        fen = en;
        fpos.SetXYZ(posx,posy,posz);
        fhitsadded = hitsadded;
    }
    virtual void Clear(){
        fid = -1;
        fring = -1;
        ftype = -1;
        fdaqid = -1;
        fpos.SetXYZ(-9999,-9999,-9999);
        frndpos.SetXYZ(-9999,-9999,-9999);
        fts = 0;
        fadc = -1;
        fen = -1;
        fhitsadded = 0;

        //for veto
        fdvetotime = -999999.;
        ff11time = -999999.;
        fvetotime = -999999.;
    }
    virtual void Copy(BELENHit& obj){
        obj.SetID(fid);
        obj.SetRing(fring);
        obj.SetType(ftype);
        obj.SetDaqID(fdaqid);
        obj.SetTimestamp(fts);
        obj.SetPos(fpos.X(),fpos.Y(),fpos.Z());
        obj.SetRndPos(frndpos.X(),frndpos.Y(),frndpos.Z());
        obj.SetADC(fadc);
        obj.SetEnergy(fen);
        obj.SetHitsAdded(fhitsadded);
        //for veto
        obj.SetF11Time(ff11time);
        obj.SetDownstreamVetoTime(fdvetotime);
        obj.SetFinalVetoTime(fvetotime);
    }

    //! Set the energy
    void SetEnergy(double energy){fen = energy;}

    //! Set the raw ADC value
    void SetADC(int adc){fadc = adc;}

    //! Set the counter ID
    void SetID(short id){fid = id;}
    //! Set the counter daq ID
    void SetDaqID(short daqid){fdaqid = daqid;}
    //! Set the counter ring
    void SetRing(short ring){fring = ring;}
    //! Set the counter type
    void SetType(short type){ftype = type;}

    //! Set the timestamp
    void SetTimestamp(unsigned long long int ts){fts = ts;}

    //! Set the He3 position
    void SetPos(Double_t x, Double_t y, Double_t z){fpos.SetXYZ(x,y,z);}

    //! Set the He3 position
    void SetRndPos(Double_t x, Double_t y, Double_t z){frndpos.SetXYZ(x,y,z);}

    //! Set current hits
    void SetHitsAdded(unsigned short hitsadded){fhitsadded = hitsadded;}


    //! for veto
    void SetF11Time(double f11time){ff11time = f11time;}
    void SetDownstreamVetoTime(double dvetotime){fdvetotime = dvetotime;}
    void SetFinalVetoTime(double vetotime){fvetotime = vetotime;}


    //! Get the ID
    short GetID(){return fid;}
    //! Get the ID
    short GetDaqID(){return fdaqid;}
    //! Get the ring (my precious!)
    short GetMyPrecious(){return fring;}
    //! Get the type
    short GetType(){return ftype;}

    //! Get the energy
    double GetEnergy(){return fen;}
    //! Get the timestamp
    unsigned long long int GetTimestamp(){return fts;}
    //! Get the raw ADC value
    int GetADC(){return fadc;}
    //! Get 3He position
    TVector3 GetPosition(){return fpos;}
    TVector3 GetRndPosition(){return frndpos;}

    //! Get current hits
    unsigned short GetHitsAdded(){return fhitsadded;}


    //! for veto
    double GetF11Time(){return ff11time;}
    double GetDownstreamVetoTime(){return fdvetotime;}
    double GetFinalVetoTime(){return fvetotime;}

    //! Printing information
    void Print(Option_t *option = "") const {
      cout << "ID " << fid;
      cout << "daq ID " << fdaqid;
      cout << "ring " << fring;
      cout << "type" <<ftype;
      cout << "\tX pos " << fpos.X();
      cout << "\tY pos " << fpos.Y();
      cout << "\tZ pos " << fpos.Z();
      cout << "\tX pos " << frndpos.X();
      cout << "\tY pos " << frndpos.Y();
      cout << "\tZ pos " << frndpos.Z();
      cout << "\tadc " << fadc;
      cout << "\tenergy " << fen;
      cout << "\ttimestamp " << fts;
      cout << "\thits added " << fhitsadded << endl;
      return;
    }

protected:

    //! current hits
    unsigned short fhitsadded;
    //! Position of 3He counter
    TVector3 fpos;
    //! Position of 3He counter with random generator
    TVector3 frndpos;
    //! daq ID number of 3He counter
    short fdaqid;
    //! physical ID number of 3He counter
    short fid;
    //! ring number of the 3He counter
    short fring;
    //! type of tube: 0: riken, 1: upc 1 inch 2: ornl 1 inch 3: ornl 2 inch
    short ftype;

    //! ADC value
    int fadc;
    //! Energy calibrated value
    double fen;
    //! timestamp value
    unsigned long long fts;

    //! special added
    double ff11time;//in us
    double fdvetotime;//in us
    double fvetotime; //in us


    /// \cond CLASSIMP
    ClassDef(BELENHit,1);
    /// \endcond
    ///
};

class BELEN : public TObject
{
public:
    //! default constructor
    BELEN(){
        Clear();
    }
    //! clear BELEN info
    virtual void Clear(){
        fmult = 0;
        //! Dealocating memory
        for (size_t idx=0;idx<fhits.size();idx++){
            delete fhits[idx];
        }
        fhits.clear();
    }

    //! Set time stamp
    void SetTimestamp(unsigned long long ts){fbelents = ts;}
    //! Set Multiplicity;
    void SetMult(unsigned short mult) {fmult = mult;}

    //! Add more hits
    void AddHits(vector<BELENHit*> hits){
      fmult += hits.size();
      for(vector<BELENHit*>::iterator hit=hits.begin(); hit!=hits.end(); hit++){
        //set hit add here!
        //
        fhits.push_back(*hit);
      }
    }

    //! Add a hit
    void AddHit(BELENHit* hit){
      //!newly added
      hit->SetHitsAdded(fmult);
      fhits.push_back(hit);
      fmult++;
    }

    //! Set all hits
    void SetHits(vector<BELENHit*> hits){
      fmult = hits.size();
      fhits = hits;
    }

    //! Returns timestamp
    unsigned long long GetTimestamp(){return fbelents;}
    //! Returns the multiplicity of the event
    unsigned short GetMult(){return fmult;}


    //! Returns the whole vector of hits
    vector<BELENHit*> GetHits(){return fhits;}
    //! Returns the hit number n
    BELENHit* GetHit(unsigned short n){return fhits.at(n);}


    void Print(Option_t *option = "") const {
        cout <<"timestamp " << fbelents << endl;
        cout << "multiplicity " << fmult << endl;
    }
protected:
    //! total multiplicity
    unsigned short fmult;
    //! the ealiest time stamp found
    unsigned long long  fbelents;
    //! vector with the hits
    vector<BELENHit*> fhits;

    /// \cond CLASSIMP
    ClassDef(BELEN,1);
    /// \endcond
};

#endif // BELLEN_H
