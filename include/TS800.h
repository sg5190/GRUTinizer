#ifndef _TS800HUNDRED_H_
#define _TS800HUNDRED_H_

#include "TClonesArray.h"
#include <iostream>
#include <vector>
#include "TDetector.h"

#include "GValue.h"
#include "TS800Hit.h"
//class TTrigger;

class TChain;

class TS800 : public TDetector {
public:
  TS800();
  virtual ~TS800();

  //////////////////////////////////////////////
  virtual void InsertHit(const TDetectorHit&);
  virtual TDetectorHit& GetHit(int i);
  /////////////////////////////////////////////

  void SetEventCounter(Long_t event) { fEventCounter = event; }

  Long_t GetEventCounter() { return fEventCounter;}
  Long_t GetTimestamp()    { return Timestamp(); }

  TVector3 CRDCTrack();  // not a finished method

  //S800 track determined from inverse map with shifts sata, sbta for ata and
  //bta, respectively. Note these shifts are applied in addition to the shifts
  //from GValues ATA_SHIFT and BTA_SHIFT. Also, the flip in direction for the
  //non-dispersive direction to be in the GRETINA reference frame is already
  //applied in this function.
  TVector3 Track(float Ata, float Bta, double sata=0.000,double sbta=0.000) const;
  TVector3 Track(double sata=0.000,double sbta=0.000) const;

  float GetXFP(int i=0) const; // x position in the first(second) CRDC (mm)
  float GetYFP(int i=0) const; // y position in the first(second) CRDC (mm)
  float GetAFP() const; // x-angle in the focal plane (rad)
  float GetBFP() const; // y-angle in the focal plane (rad)

  float GetAFP(float, float) const; // x-angle in the focal plane (rad)
  float GetBFP(float, float) const; // y-angle in the focal plane (rad)

  Float_t GetAta(int i=6) const; // x-angle at the target (rad)
  Float_t GetYta(int i=6) const; // y-offset at the target (mm)
  Float_t GetBta(int i=6) const; // y-angle at the target (rad)
  Float_t GetDta(int i=6) const; // dE/E of outgoing particle, relative to the central b-rho

  Float_t GetAta(float, float, float, float, int i=6) const; // x-angle at the target (rad)
  Float_t GetYta(float, float, float, float, int i=6) const; // y-offset at the target (mm)
  Float_t GetBta(float, float, float, float, int i=6) const; // y-angle at the target (rad)
  Float_t GetDta(float, float, float, float, int i=6) const; // dE/E of outgoing particle, relative to the central b-rho

  Float_t GetScatteringAngle(int i=6)	  	    const { return GetScatteringAngle(GetAta(i),GetBta(i));}
  Float_t GetScatteringAngle(float, float)  const;	//Projectile Scattering Angle

  Float_t GetAzita(int i=6)			    const {return GetAzita(GetAta(i), GetBta(i));}
  Float_t GetAzita(float, float)   const;	//Phi polar angle at target, in radians

  float AdjustedBeta(float beta)                    const {return AdjustedBeta(beta, GetDta());}
  float AdjustedBeta(float, float) const;

  virtual void Copy(TObject& obj)        const;
  virtual void Print(Option_t *opt = "") const {;}
  virtual void Clear(Option_t* opt = "");

  TCrdc         &GetCrdc(int x=0)  const { if(x==0) return (TCrdc&)crdc1;
                                           else return (TCrdc&)crdc2;}
  TTof          &GetTof()          const { return (TTof&)tof;        }
  TMTof         &GetMTof()         const { return (TMTof&)mtof;        }
  TIonChamber   &GetIonChamber()   const { return (TIonChamber&)ion; }
  TScintillator &GetScint(int x=0) const { return (TScintillator&)scint[x]; }
  //Added 1/4/2016 for getting trigbit - BAE
  TTrigger      &GetTrigger()      const { return (TTrigger&)trigger;}
  THodoscope &GetHodoscope() const { return (THodoscope&)hodo;}


  float GetIonSum() const { return ion.GetSum(); }

  //Note c1 is the AFP correction and c2 is the XFP correction
  float GetTofE1_TAC(float c1=0.00,float c2=0.00)  const;
  float GetTofXFPE1_TAC(float c1=0.00,float c2=0.00)  const;
  float GetTofE1_TDC(float c1=0.00,float c2=0.00)  const;
  float GetTofXFP_E1_TDC(float c1=0.00,float c2=0.00)  const;
  float GetTofE1_MTDC(float c1=0.00,float c2=0.00,int i=0) const;

  float GetCorrTOF_OBJTAC() const;
  float GetCorrTOF_XFPTAC() const;
  float GetOBJRaw_TAC() const;
  float GetXFRaw_TAC() const;

  float GetCorrTOF_OBJ() const;
  float GetCorrTOF_XFP() const;
  float GetOBJ_E1Raw() const;
  float GetXF_E1Raw() const;
  //===================================================

  float GetCorrTOF_OBJ_MESY(int i=0) const;

  float GetOBJ_E1Raw_MESY(int i=0) const;
  float GetOBJ_E1Raw_MESY_Ch15(int i=0) const;

  float GetXF_E1Raw_MESY(int i=0) const;
  float GetXF_E1Raw_MESY_Ch15(int i=0) const;

  float GetRawOBJ_MESY(unsigned int i=0) const;

  float GetRawE1_MESY(unsigned int i=0) const;
  float GetRawE1_MESY_Ch15(unsigned int i=0) const;

  float GetRawXF_MESY(unsigned int i=0) const;

  // Return the correlated gvalue corrected time-of-flight obj to e1.
  // If no corrections are passed, GValues OBJ_MTOF_CORR_AFP and
  // OBJ_MTOF_CORR_XFP are used
  double GetMTofObjE1() const ;
  double GetMTofObjE1(double afp_cor, double xfp_cor) const ;
  // Return the correlated gvalue corrected time-of-flight xfp to e1.
  // If no corrections are passed, GValues XFP_MTOF_CORR_AFP and
  // XFP_MTOF_CORR_XFP are used
  double GetMTofXfpE1() const ;
  double GetMTofXfpE1(double afp_cor, double xfp_cor) const ;

  double GetMTofCorr(double correlatedtof, double afp, double xfp, double afp_cor, double xfp_cor) const ;

  unsigned short GetME1Up(int i)       const {if(i<GetME1Size())   return mtof.fE1Up.at(i);  return sqrt(-1); }
  unsigned short GetME1Down(int i)     const {if(i<GetME1Size())   return mtof.fE1Down.at(i);return sqrt(-1); }
           int   GetME1Size()          const { return mtof.fE1Up.size();  }
  unsigned short GetMXfp(int i)        const {if(i<GetMXfpSize())    return mtof.fXfp.at(i);   return sqrt(-1); }
           int   GetMXfpSize()         const { return mtof.fXfp.size();  }
  unsigned short GetMObj(int i)        const {if(i<GetMObjSize())    return mtof.fObj.at(i);   return sqrt(-1); }
           int   GetMObjSize()         const { return mtof.fObj.size();  }
  unsigned short GetMRf(int i)         const {if(i<GetMRfSize())     return mtof.fRf.at(i);    return sqrt(-1);  }
           int   GetMRfSize()          const { return mtof.fRf.size();  }
  unsigned short GetMCrdc1Anode(int i) const {if(i<GetMCrdc1AnodeSize()) return mtof.fCrdc1Anode.at(i);return sqrt(-1); }
           int   GetMCrdc1AnodeSize()  const { return mtof.fCrdc1Anode.size();  }
  unsigned short GetMCrdc2Anode(int i) const {if(i<GetMCrdc2AnodeSize()) return mtof.fCrdc2Anode.at(i);return sqrt(-1); }
           int   GetMCrdc2AnodeSize()  const { return mtof.fCrdc2Anode.size();  }
  unsigned short GetMHodoscope(int i)  const {if(i<GetMHodoscopeSize()) return mtof.fHodoscope.at(i);  return sqrt(-1);  }
           int   GetMHodoscopeSize()   const { return mtof.fHodoscope.size();  }
  unsigned short GetMRef(int i)        const {if(i<GetMRefSize())       return mtof.fRef.at(i);        return sqrt(-1);  }        
           int   GetMRefSize()         const { return mtof.fRef.size();  }


  int GetReg() const { return trigger.GetRegistr(); }
  float GetCorrIonSum() const;


public:

private:
  virtual int  BuildHits(std::vector<TRawEvent>& raw_data);

  bool HandleTrigPacket(unsigned short*,int);     //!
  bool HandleTOFPacket(unsigned short*,int);      //!
  bool HandleScintPacket(unsigned short*,int);    //!
  bool HandleIonCPacket(unsigned short*,int);     //!
  bool HandleCRDCPacket(unsigned short*,int);     //!
  bool HandleMTDCPacket(unsigned short*,int);     //!
  bool HandleHodoPacket(unsigned short*,int);     //!

  TScintillator scint[3];
  TTrigger     trigger;
  TTof         tof;
  TMTof        mtof;
  TIonChamber  ion;
  TCrdc        crdc1;
  TCrdc        crdc2;
  THodoscope   hodo;
  //TMultiHitTof multi_tof;

  Long_t fEventCounter;
  static bool fGlobalReset; //!


  public:

  void SetGlobalReset(bool flag=true) { fGlobalReset=flag; }

  ClassDef(TS800,3);
};

class TS800Track {
  public:
    TS800Track();
    virtual ~TS800Track();
    void Clear();
    //Uses the CRDC X and Y postion and inverse map to calculate all Focal plane and target angles
    void CalculateTracking(const TS800 *, int i = 6);
  public:
    float xfp[2];   // Focal plane dispersive direction position, in m
    float yfp[2];   // FP non-dispersive direction, in m
    float afp;      // FP dispersive angle, in radians
    float bfp;      // FP non-dispersive angle, in radians
    float ata;      // Dispersive angle at the target, in radians
    float yta;      // Non-dispersive position at the target, in m
    float bta;      // Non-dispersive angle at the target, in radians
    float dta;      // Fractional energy at the target, (E-E0)/E0, in parts
    float azita;    // Phi polar angle at target, in radians
    float scatter;  // Theta polar angle at target, in mradians
    TVector3 track; // S800 Track from inverse map
  ClassDef(TS800Track,1)
};

#endif
