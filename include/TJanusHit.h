#ifndef _TJANUSDETECTORHIT_H_
#define _TJANUSDETECTORHIT_H_

#include "TDetectorHit.h"

class TJanusHit : public TDetectorHit {
public:
  TJanusHit() { }
  TJanusHit(const TJanusHit& hit);

  TJanusHit& operator=(const TJanusHit& hit);

  virtual Int_t Charge() const;

  void Clear(Option_t* opt = "");
  void Print(Option_t* opt = "") const;
  void Copy(TObject& obj) const;

  void SetRingNumber(int ring) { fRing = ring; }
  void SetSectorNumber(int sector) { fSector = sector; }
  void SetBackCharge(int charge) { fBackCharge = charge; }

  int GetDetnum() const;
  int GetRing() const { return fRing; }
  int GetSector() const { return fSector; }
  Double_t BackCharge() const;
  double GetDefaultDistance() const;

  bool IsDownstream() const {return (GetDetnum() > 0); }

  TVector3 GetPosition(bool apply_array_offset = true) const;
  TVector3 GetPosition(bool smear, bool apply_array_offset) const;

private:

  int fRing;
  int fSector;
  double fBackCharge;
  ClassDef(TJanusHit,4);
};

#endif /* _TJANUSDETECTORHIT_H_ */
