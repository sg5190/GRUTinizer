#ifndef _TSEGAHIT_H_
#define _TSEGAHIT_H_

#include "TDetectorHit.h"
#include "TSegaSegmentHit.h"

#include "TMath.h"

#define MAX_TRACE_LENGTH 100

class TSegaHit : public TDetectorHit {
public:
  TSegaHit();

  virtual void Copy(TObject&) const;
  virtual void Clear(Option_t *opt = "");
  virtual void Print(Option_t *opt = "") const;

  virtual Int_t Charge() const;

  Int_t GetDetnum() const;
  int GetMainSegnum() const;

  bool HasCore() const;

  unsigned int GetNumSegments() const { return fSegments.size(); }
  TSegaSegmentHit& GetSegment(int i) { return fSegments.at(i); }
  unsigned long GetSegmentTimestamp() {
    if(fSegments.size()){
      return fSegments[0].Timestamp();
    } else {
      return -1;
    }
  }

  std::vector<unsigned short>* GetTrace(int segnum=0);

  void SetCFDFail(int cfdbit) { fCfdFail = cfdbit; }
  int GetCFDFail() const { return fCfdFail; }

//  void SetTime(double time)            { fTime = time; }
  void SetTrace(unsigned int trace_length, const unsigned short* trace);

  void SetEnergySumBool(bool sum = false) {fEsum = sum; }
  bool HasEnergySum() const {return fEsum; }

  void SetEnergySum(int sum, int energy_sum) { fEnSum[sum] = energy_sum; }
  int GetEnergySum(int sum) const { return fEnSum[sum]; }


  TSegaSegmentHit& MakeSegmentByAddress(unsigned int address);

  int GetSlot() const;
  int GetCrate() const;
  int GetChannel() const;

  TVector3 GetPosition(bool apply_array_offset = false,
                       TVector3 array_offset = TVector3(sqrt(-1.0),sqrt(-1.0),sqrt(-1.0))) const;

  double GetDoppler(double beta,
                    const TVector3& particle_vec = TVector3(0,0,1),
                    const TVector3& sega_offset = TVector3(sqrt(-1.0),sqrt(-1.0),sqrt(-1.0))) const;

  double GetTraceHeight() const;
  double GetTraceHeightDoppler(double beta,const TVector3& vec = TVector3(0,0,1)) const;


private:
  std::vector<unsigned short> fTrace;
  std::vector<TSegaSegmentHit> fSegments;
  int fEnSum[4];
  bool fEsum;
  int fCfdFail;
  ClassDef(TSegaHit,4);
};

#endif /* _TSEGAHIT_H_ */
