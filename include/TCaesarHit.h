#ifndef _TCAESARHIT_H_
#define _TCAESARHIT_H_

#include "TDetectorHit.h"
#include "TMath.h"


class TCaesarHit : public TDetectorHit {
public:
  TCaesarHit();

  virtual void Copy(TObject&) const;
  virtual void Clear(Option_t *opt = "");
  virtual void Print(Option_t *opt = "") const;

  void SetTime(unsigned short time)            { fTimes.push_back(time); }
  void SetChannelNum(int chan)	       	       { fChan = chan;   }

  int GetChannel() const { return fChan; };
  int GetRing() const;
  int GetDetectorNumber() const;
  bool IsOverflow() const { return fOverflow; };
  bool IsUnderflow() const { return fUnderflow; };

  TVector3 GetPosition(bool apply_array_offset = false, TVector3 arrayOffset = TVector3(sqrt(-1.0), sqrt(-1.0), sqrt(-1.0))) const;
  unsigned short GetTime0() const { if(fTimes.size()) return fTimes.at(0);
			            else return -1;
			          }
  unsigned short GetMTime(int i) const { if((int)fTimes.size() > i) return fTimes.at(i);
			                 else return - 1;
			               }
  size_t GetNTimes() const { return fTimes.size(); }
  int GetABDepth() const { return fAB; }
  void SetABDepth(int ab) const { fAB = ab; }
  void SetOverflow(bool over) { fOverflow = over; }
  void SetUnderflow(bool under) { fUnderflow = under; }

  virtual Int_t Charge() const;
  void NNAdd(const TCaesarHit& other);

  protected:
  int fChan;
  bool fOverflow;
  bool fUnderflow;
  mutable int fAB;

  std::vector<unsigned short> fTimes;
  ClassDef(TCaesarHit,4);
};

#endif /* _TCAESARHIT_H_ */
