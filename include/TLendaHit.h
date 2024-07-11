#ifndef _TLENDAHIT_H_
#define _TLENDAHIT_H_

#include <TVector3.h>
#include <TDetectorHit.h>
#include <TGraph.h>
#include <TF1.h>
#include <TFitResult.h>

class TLendaHit : public TDetectorHit {
  public:
    TLendaHit() { }
    ~TLendaHit() { }

    virtual void Copy(TObject& obj) const        { TDetectorHit::Copy(obj); }
    virtual void Clear(Option_t *opt = "");//    { TDetectorHit::Clear(opt); }
    virtual void Print(Option_t *opt = "") const { TDetectorHit::Print(opt); }

    std::vector<unsigned short>* GetTrace()   	{ return &fTrace; }
//    std::vector<Double_t>* GetFF()  		{ return &fFastFilter; }

    void SetTrace(unsigned int trace_length, const unsigned short* trace);

    int Size() { return fTrace.size(); }

    void SetPosArray(int arr_pos) { fArrayPosition = arr_pos; }
    void SetPosSeg(int seg_pos) { fSegmentPosition = seg_pos; }

    void SetCFDFail(int cfdbit) { fCfdFail = cfdbit; }
    int GetCFDFail() const { return fCfdFail; }
    void SetExtTimestamp(long timestamp) { fExtTimestamp = timestamp; }
    long GetExtTimestamp() const      { return fExtTimestamp; }

    int GetSlot()    const;
    int GetCrate()   const;
    int GetChannel() const;
    Int_t GetDetnum() const { return fDetector; }
    TVector3 GetPosition() { return TVector3(0,0,1); }

    int GetPosArray() const { return fArrayPosition;   }
    int GetPosSeg  () const { return fSegmentPosition; }

    Double_t numOfBadFits;

  private:
    std::vector<unsigned short> fTrace;
    int fDetector;
    int fArrayPosition;
    int fSegmentPosition;
    int fCfdFail;
    long fExtTimestamp;
  ClassDef(TLendaHit,2)
};

#endif
