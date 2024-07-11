#ifndef _TSUNHIT_H_
#define _TSUNHIT_H_

#include "TDetectorHit.h"
#include "TMath.h"

class TSunHit : public TDetectorHit {
public:
  TSunHit();

  virtual void Copy(TObject&) const;
  virtual void Clear(Option_t *opt = "");
  virtual void Print(Option_t *opt = "") const;

  std::vector<unsigned short>* GetTrace()     { return &fTrace; }

  void SetTime(double time)            { fTime = time; }
  void SetTrace(unsigned int trace_length, const unsigned short* trace);
  void SetExtTime(long timestamp)  { fExtTime = timestamp; }

  virtual Int_t Charge() const;

  int GetSlot() const;
  int GetCrate() const;
  int GetChannel() const;
  long GetExtTime() const      { return fExtTime; }

  protected:
    std::vector<unsigned short> fTrace;
    long fExtTime;
  ClassDef(TSunHit,4);
};

#endif /* _TSUNHIT_H_ */
