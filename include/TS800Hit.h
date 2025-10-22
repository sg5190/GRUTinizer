#ifndef TSEIGHTHUNDRADHIT_H
#define TSEIGHTHUNDRADHIT_H

#include <TObject.h>
#include <TClass.h>
#include <iostream>

#include <TDetectorHit.h>
#include <TMath.h>
#include "TRandom.h"
#include "GValue.h"


#define MAXCRDC 513

class TS800Channel : public TDetectorHit {
  public:
    TS800Channel()                     { Clear();    }
    TS800Channel(short value)          { Set(value); }
    TS800Channel(unsigned short value) { Set(value); }
    ~TS800Channel()                    {}

    void Set(short value)          { fValue = (unsigned short)value; }
    void Set(unsigned short value) { fValue = value; }

    virtual short GetId()      const { return (((fValue)&0xf000)>>12);  }
    short GetValue()           const { return fValue&0x0fff;  }

    virtual void Clear(Option_t *opt="")       { TDetectorHit::Clear(opt); fValue = 0; }
    virtual void Print(Option_t *opt="") const { printf("[%i] = %i\n",GetId(),GetValue());}
    virtual void Copy(TObject &obj)      const { TDetectorHit::Copy(obj); ((TS800Channel&)obj).fValue = fValue; }

    virtual int  GetCharge() const { return GetValue(); }

  private:
    //void Init() { Class()->IgnoreTObjectStreamer(true); }
    short fValue;

  ClassDef(TS800Channel,1);
};

class TTrigger : public TDetectorHit {
  public:
    TTrigger();
    ~TTrigger();

    void SetRegistr(unsigned short reg) { fregistr=reg; }
    void SetS800Source(short sou)       { fs800source=sou; }
    void SetExternalSource1(short sou)  { fexternalsource1=sou; }
    void SetExternalSource2(short sou)  { fexternalsource2=sou; }
    void SetSecondarySource(short sou)  { fsecondarysource=sou; }

    unsigned short GetRegistr() const { return fregistr; }
    short GetS800Source()       const { return fs800source; }
    short GetExternalSource1()  const { return fexternalsource1; }
    short GetExternalSource2()  const { return fexternalsource2; }
    short GetSecondarySource()  const { return fsecondarysource; }

    virtual void Copy(TObject &)         const;
    virtual void Print(Option_t *opt="") const;
    virtual void Clear(Option_t *opt="");

  private:
    virtual int GetCharge() const { return 0; }

    unsigned short fregistr;
    short fs800source;
    short fexternalsource1;
    short fexternalsource2;
    short fsecondarysource;

  ClassDef(TTrigger,1)
};

class TCrdc : public TDetectorHit {
  public:
    TCrdc();
    ~TCrdc();

    short GetId()    { return fId;   }
    short GetAnode() { return anode; }
    short GetTime()  { return time;  }
    float GetTimeRand()  { return ((float)(time)+gRandom->Uniform());  }

    int  GetSize()        const { return channel.size(); }
    int  GetNSamples() const { return sample.size(); }

    void SetId(short id)    { fId = id;  }
    void SetAnode(short an) {anode = an; }
    void SetTime(short ti)  {time = ti;  }


    int  GetAddress(int i) const { return TDetectorHit::GetAddress() + channel.at(i); }
    void AddPoint(int chan,int samp,int dat) { channel.push_back(chan);
                                               sample.push_back(samp);
                                               data.push_back(dat);    }
    int GetChannel(int i) const    { if(i>=GetSize()) return -1; return channel.at(i);    }
    int GetSample(int i)  const    { if(i>=GetSize()) return -1; return sample.at(i);     }
    int GetData(int i)    const    { if(i>=GetSize()) return -1; return data.at(i);       }

    int GetWidth();

    float GetDispersiveX() const;
    float GetNonDispersiveY();
    float GetDispersiveX(int) const;
    float GetNonDispersiveY(int);

    int GetMaxPad() const;
    int GetMaxPadFast() const;
    int GetMaxPadSum() const;

    double GetDrifttime() const { return fDrifttime; }
    void SetDrifttime(double drifttime) { if(fDrifttime < 10) fDrifttime = drifttime; }

    virtual void Copy(TObject&) const;
    virtual void Print(Option_t *opt="") const;
    virtual void Clear(Option_t *opt="");

    int Sum() const { int result = 0; for(unsigned int x=0;x<data.size();x++) result +=data[x]; return result; }

  private:
    virtual int GetCharge() const { return 0; }

    bool IsGoodSample(int i) const;

    short fId;
    std::vector<int> channel;
    std::vector<int> sample;
    std::vector<int> data;

    unsigned short anode;
    unsigned short time;
    mutable bool has_cached_dispersive_x; //!
    mutable double cached_dispersive_x; //!
    double fDrifttime;

  ClassDef(TCrdc,1)
};

class TIonChamber : public TDetectorHit {
  public:
    TIonChamber();
    ~TIonChamber();

    void Set(int ch, int data); // { fChan.push_back(ch); fData.push_back(data); }

    int GetChannel(int i) const { if(i>=GetSize()) return -1; return fChan.at(i); }
    int GetData(int i)    const { if(i>=GetSize()) return -1; return fData.at(i); }
    int GetSize() const { return fChan.size(); }
    float GetdE(TCrdc *);
    float GetdE(double crdc_1_x, double crdc_1_y);
    float GetSum() const;
    float GetAve();

    int  GetAddress(int i) const { return TDetectorHit::GetAddress() + GetChannel(i); }


    float GetCalData(int i) const {
      TChannel *c = TChannel::GetChannel(GetAddress(i));
      if(c){
        return c->CalEnergy(GetData(i));
      }else{
       return (float)GetData(i);
      }
    }

    virtual void Copy(TObject&) const;
    virtual void Print(Option_t *opt="") const;
    virtual void Clear(Option_t *opt="");
    int GetCharge() const { int sum=0;for(int i=0;i<GetSize();i++)sum+=GetData(i);return sum;}
    //int GetSum() const { int sum=0;for(int i=0;i<GetSize();i++)sum+=GetData(i);return sum;}
  private:


    std::vector<int> fChan;
    std::vector<int> fData;
    //float            fdE;   //!          Average of all fData.
    ClassDef(TIonChamber,1)
};

class THodoHit : public TDetectorHit{
  public:
    THodoHit() { Clear(); }
    THodoHit(const THodoHit&);

    virtual void Copy(TObject &) const;
    virtual void Print(Option_t *opt="") const;
    virtual void Clear(Option_t *opt="");


    void SetChannel(int channel) { fChannel = channel;}

    int GetCharge() const { return GetCharge(); }
    int GetChannel() const { return fChannel; }
  private:
    int fChannel;

    ClassDef(THodoHit, 1);
};

class THodoscope : public TDetectorHit {
  public:
    THodoscope();
    virtual ~THodoscope();
    THodoscope(const THodoscope&);

    virtual void Copy(TObject&) const;
    virtual void Print(Option_t *opt="") const;
    virtual void Clear(Option_t *opt="");
    void InsertHit(const TDetectorHit&);

    int  GetAddress(int i) const { return TDetectorHit::GetAddress() + GetChannel(i); }

    int GetChannel(int i) const { return GetHodoHit(i).GetChannel();}
    TDetectorHit& GetHit(int i);
    THodoHit& GetHodoHit(int i);
    const THodoHit& GetHodoHit(int i) const;
    std::vector<THodoHit> GetHodoHits() const { return hodo_hits; } ;

    size_t GetSize() const { return hodo_hits.size(); }

    std::vector<THodoHit> hodo_hits;

    ClassDef(THodoscope, 1);

};

class TMTof : public TDetectorHit {
  public:
    TMTof();
    ~TMTof();
    TMTof(const TMTof&);

    virtual void Copy(TObject&) const;
    virtual void Print(Option_t *opt="") const;
    virtual void Clear(Option_t *opt="");

    int E1UpSize()       const { return fE1Up.size();   }
    int E1DownSize()     const { return fE1Down.size();  }
    int XfpSize()        const { return fXfp.size();      }
    int ObjSize()        const { return fObj.size();       }
    int RfSize()         const { return fRf.size();        }
    int Crdc1AnodeSize() const { return fCrdc1Anode.size();}
    int Crdc2AnodeSize() const { return fCrdc2Anode.size();}
    int HodoSize()       const { return fHodoscope.size(); }
    int RefSize()        const { return fRef.size(); }

    //Generic correlated Time of flight functions, for multihit TDC's, can
    double GetCorrelatedTof(int ch1, int ch2, double target, double shift) const;
    double GetCorrelatedTof(std::string ch1, std::string ch2, double target, double shift) const;
    //Selected MTDC vectors (fE1Up, fOb  etc.)
    std::vector<unsigned short> GetMTofVectorFromString(std::string) const;
    std::vector<unsigned short> GetMTofVector(int i = 15) const;

  private:

  public:
    std::vector<unsigned short> fE1Up;         // Channel 0
    std::vector<unsigned short> fE1Down;       // Channel 1
    std::vector<unsigned short> fXfp;          // Channel 2
    std::vector<unsigned short> fObj;          // Channel 3
    std::vector<unsigned short> fRf;           // Channel 5
    std::vector<unsigned short> fCrdc1Anode;   // Channel 6
    std::vector<unsigned short> fCrdc2Anode;   // Channel 7
    std::vector<unsigned short> fHodoscope;    // Channel 12
    std::vector<unsigned short> fRef;          // Channel 15, same as E1Up (different cable.)
    virtual Int_t GetCharge() const  {return 0;}
  ClassDef(TMTof,1)
};

class TS800Hit : public TDetectorHit {
  public:
    TS800Hit()   {  }
    ~TS800Hit()  {  }
    virtual Int_t Charge() const { return 0; }

    virtual void Clear(Option_t *opt ="")       { TDetectorHit::Clear(opt); }
    virtual void Print(Option_t *opt ="") const {  }
    virtual void Copy(TObject &obj)       const { TDetectorHit::Copy(obj);  }

  private:

  ClassDef(TS800Hit,1);
};

#endif
