#ifndef TGRETINAHIT_H
#define TGRETINAHIT_H

#include <TObject.h>
#include <Rtypes.h>
#include <TVector3.h>
#include <TMath.h>
#include <TChain.h>
#include <TRandom.h>
#include <cmath>

#include "TDetectorHit.h"

#define MAXHPGESEGMENTS 36

class TSmartBuffer;

/*
class TInteractionPoint {
  public:
    TInteractionPoint() { }
    TInteractionPoint(const TInteractionPoint &IP);
    TInteractionPoint(int seg,float eng,float frac,TVector3 lpos) :
                      fSegNum(seg),fEng(eng),fDecompEng(frac),fLPosition(lpos) { }
    virtual ~TInteractionPoint() { }

    virtual void Copy(const TInteractionPoint&);
    void Add(TInteractionPoint&);

    virtual int   GetSegNum()              const { return fSegNum;    }
    virtual float GetPreampE()             const { return fEng;       }
    virtual float GetDecompE()             const { return fDecompEng; }
    virtual float GetAssignE()             const { return fAssignedEng; }
    virtual int   GetOrder()               const { return fOrder; }
    virtual TVector3 GetPosition(int xtal) const;
    TVector3 GetLocalPosition()            const { return fLPosition; }
    void SetOrder(int o)     { fOrder=o; }
    void SetAssignE(float e) { fAssignedEng = e; }

    void Print(Option_t *opt="") const;
    void Clear(Option_t *opt="");

    void SetSegNum(int seg) { fSegNum = seg; }

    int Wedge() const { return ((GetSegNum()-1)%6); }

    bool operator == (const TInteractionPoint &rhs) const {
      return (GetDecompE() == rhs.GetDecompE());
    }
    bool operator < (const TInteractionPoint &rhs) const{
      return (GetDecompE() > rhs.GetDecompE());
    }

  private:
    int   fSegNum;
    float fEng;          // energy as recorded by preamp.  energy in mode2 format
    float fDecompEng;    // energy as assigned by decomp.  fraction in mode2 format
    float fAssignedEng;  // percent eng assigned by decomp scaled to the core. not in mode2 format
    int   fOrder;        // interaction order
    //TVector3 fPosition;  // in global coordinates
    TVector3 fLPosition; // in local coordinates.
  ClassDef(TInteractionPoint,1)

};
*/
class interaction_point {
  public:
  interaction_point():fSeg(-1),fX(sqrt(-1)),fY(sqrt(-1)),fZ(sqrt(-1)),fEng(sqrt(-1)),fFrac(sqrt(-1)) { }
  interaction_point(int seg,float x,float y,float z,float energy,float fraction=100.0)
    : fSeg(seg),fX(x),fY(y),fZ(z),fEng(energy),fFrac(fraction) { }
  virtual ~interaction_point() { }
  int   fSeg;
  float fX;
  float fY;
  float fZ;
  float fEng;
  float fFrac;

  bool operator<(const interaction_point &other) const {
    if(fEng!=other.fEng) {
      return fEng>other.fEng;
    }
    if(fSeg==other.fSeg) {
      return fFrac>other.fFrac;
    }
    return fSeg<other.fSeg;
  }

  void Print(Option_t *opt="") const {
    printf("Seg[%02i]\tWedge[%i]\tEng: % 4.2f / % 4.2f  \t(X,Y,Z) % 3.2f % 3.2f % 3.2f\n",
            fSeg,(fSeg%6),fEng,fFrac,fX,fY,fZ);
  };
  ClassDef(interaction_point,1)
};

#ifndef __CINT__

struct intpnt_compare {
  bool operator()(const interaction_point &p1,const interaction_point &p2) const {
    return p1.fSeg < p2.fSeg;
  }
};

struct intpnt_compare_wedge {
  bool operator()(const interaction_point &p1,const interaction_point &p2) const {
    return (p1.fSeg%6) < (p2.fSeg%6);
  }
};

#endif

class TGretaHit : public TDetectorHit {

public:
  TGretaHit();
  TGretaHit(const TGretaHit& hit){ hit.Copy(*this); };
  ~TGretaHit();

  void Copy(TObject& obj) const;

  void BuildFrom(TSmartBuffer& raw);


  Double_t GetTime()             const { return ((double)GetTimestamp() + (double)fT0)*10.; } //Return Time in ns
  Float_t  GetT0()               const { return fT0; }
  Int_t    GetCrystalId()        const { return fCrystalId;      }
  Int_t    GetHoleNumber()       const { return fCrystalId/4-1;  }
  Int_t    GetCrystalNumber()    const { return fCrystalId%4;    }
  Int_t    GetSpecOrder()        const;
  Int_t    GetInverseSpecOrder() const;
  Float_t  GetBaseline()	 const { return fBaseline; 	  }
  Float_t  GetPoststep()	 const { return fPoststep;	  }
  Float_t  GetPrestep()		 const { return fPrestep;	  }
  Float_t  GetChi2()		 const { return fChisq;           }
  Float_t  GetNormChi2()	 const { return fNormChisq;       }
  Int_t GetPad()                 const { return fPad; }

  const char *GetName() const;

  void  Print(Option_t *opt="") const;
  void  Clear(Option_t *opt="");

  double GetX() const { return GetPosition().X(); }
  double GetY() const { return GetPosition().Y(); }
  double GetZ() const { return GetPosition().Z(); }

  int GetABDepth() const { return fAB; }

  double GetPhi() const {
    double phi = GetPosition().Phi();
    if(phi<0) {
      return TMath::TwoPi()+phi;
    } else {
      return phi;
    }
  }
  double GetTheta()    const { return GetPosition().Theta(); }
  double GetPhiDeg()   const { return GetPhi()*TMath::RadToDeg(); }
  double GetThetaDeg() const { return GetTheta()*TMath::RadToDeg(); }

  Int_t Compare(const TObject *obj) const {
    TGretaHit *other = (TGretaHit*)obj;
    if(this->GetEnergy()>other->GetEnergy())
      return -1;
    else if(this->GetEnergy() < other->GetEnergy())
      return 1;  //sort largest to smallest.
    return 0;
  }

  bool HasInteractions() { return !fSegments.empty(); }

  bool operator<(const TGretaHit &rhs) const { return GetEnergy() > rhs.GetEnergy(); }

  double GetDoppler(double beta,const TVector3 *vec=0) const;

  double GetDopplerYta(double beta , double yta, const TVector3 *vec=0) const;
  double GetDopplerYta(double beta , double yta, double target_x_shift, double target_y_shift, double target_z_shift, const TVector3 *vec=0) const;

  Int_t    GetNSegments()                const { return (int)fSegments.size(); }
  Int_t    GetSegmentId(int i=-1)        const { if(i>=GetNSegments()||GetNSegments()==0) return -1;
                                                 if(i==-1) return fSegments.at(0).fSeg;
                                                 else	   return fSegments.at(i).fSeg;  }
  Float_t  GetSegmentEng(const int &i)   const { return fSegments.at(i).fEng;  }

  TVector3 GetIntPosition(unsigned int i)   const;  // position of the ith segment, Global coor.
  TVector3 GetLocalPosition(unsigned int i) const;  // position of the ith segment, Local coor.
  TVector3 GetPosition()                    const;
  TVector3 GetLastPosition()                const;

  //Set Smearing flag for position Used for UCGretina data
  void SetSmearWidth(double smearsig)	    const { fSmearWidth = smearsig; }
  void SetSmear(bool smear)	    const { fSmear = smear; }
  TVector3 GetCrystalPosition()           const;

  void NNAdd(const TGretaHit& other);
  void SetABDepth(int ab) const { fAB = ab; }
  void TrimSegments(int type); // 0: drop multiple ident int pnts.  1: make into wedge "data"
  bool IsClean() const { return !fPad; }

private:
/***********************************************************/
/* All possible decomp information and  where is is stored */
/* Int_t     type ******************************************/
/* Int_t     crystal_id  -> TGreintaHit.fCrystalId *********/
/* Int_t     num         -> TGreintaHit.fSegments.size() ***/
/* Int_t     core_e[4]   -> TGretaHit.fCoreCharge[4] *******/
/* Long_t    timestamp   -> TDetectorHit.fTimestamp ********/
/* Long_t    trig_time *************************************/
/* Float_t   t0          -> TGretaHit.fT0 ******************/
/* Float_t   cfd         -> TDetectorHit.fCfd **************/
/* Float_t   chisq       -> TGretaHit.fChisq ***************/
/* Float_t   norm_chisq ************************************/
/* Float_t   baseline    -> TGretaHit.fBaseline ************/
/* Float_t   prestep     -> TGretaHit.fPrestep *************/
/* Float_t   poststep    -> TGretaHit.fPoststep ************/
/* Int_t     pad         -> TGretaHit.fPad *****************/
/***********************************************************/

  Int_t           fCrystalId;
  Int_t           fCoreCharge[4];
  Int_t   	  fPad;
  Float_t         fT0;
  Float_t	  fChisq;
  Float_t	  fNormChisq;
  Float_t	  fBaseline;
  Float_t	  fPoststep;
  Float_t	  fPrestep;
  mutable Int_t   fAB;

  std::vector<TGretaHit> fSingles;
  bool fSetFirstSingles = false;
  std::vector<interaction_point> fSegments;
  static double fSmearWidth;
  static bool fSmear;
  ClassDef(TGretaHit,5)
};


#endif
