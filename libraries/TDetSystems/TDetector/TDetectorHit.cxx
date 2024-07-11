#include "TDetectorHit.h"

#include <iostream>
#include <cmath>

#include <TClass.h>
#include <TBuffer.h>
#include "TRandom.h"

ClassImp(TDetectorHit)

/*******************************************************************************/
/* TDetectorHit ****************************************************************/
/* Base class for all detector hits  *******************************************/
/*******************************************************************************/
const TVector3 TDetectorHit::BeamUnitVec(0,0,1);

TDetectorHit::TDetectorHit() {
  Class()->CanIgnoreTObjectStreamer();
  Clear();
}

TDetectorHit::~TDetectorHit() { }

/*******************************************************************************/
/* Clears hit  *****************************************************************/
/*******************************************************************************/
void TDetectorHit::Clear(Option_t *opt) {
  TObject::Clear(opt);
  fAddress = -1;
  fCharge = -1;
  fTime = -1;
  fTimestamp = -1;
  fFlags = 0;
}

/*******************************************************************************/
/* Virtual Print Function ******************************************************/
/* Print function should be defined in derived classes *************************/
/*******************************************************************************/
void TDetectorHit::Print(Option_t *opt) const { }

/*******************************************************************************/
/* Copy hit ********************************************************************/
/*******************************************************************************/
void TDetectorHit::Copy(TObject& obj) const {
  TObject::Copy(obj);

  TDetectorHit& hit = (TDetectorHit&)obj;
  hit.fAddress = fAddress;
  hit.fCharge = fCharge;
  hit.fTime = fTime;
  hit.fTimestamp = fTimestamp;
  hit.fFlags = fFlags;
}

/*******************************************************************************/
/* Sets detector charge + Energy ***********************************************/
/*******************************************************************************/
void TDetectorHit::SetCharge(int charge) {
  fCharge = charge + gRandom->Uniform();
  fFlags &= ~kIsEnergy;
}

void TDetectorHit::SetEnergy(double energy) {
  fCharge = energy;
  fFlags |= kIsEnergy;
}

Int_t  TDetectorHit::Charge() const {
  if(fFlags & kIsEnergy) {
    return 0;
  } else {
    return fCharge;
  }
}

Int_t TDetectorHit::GetDetnum() const {
  TChannel* chan = TChannel::GetChannel(fAddress);
  int output = -1;
  if(chan && fAddress !=- 1){
    output = chan->GetArrayPosition();
  } else {
    output = -1;
  }
  return output;
}

double TDetectorHit::GetEnergy() const {
  if(fFlags & kIsEnergy) {
    return fCharge;
  } else {
    TChannel* chan = TChannel::GetChannel(fAddress);
    if(!chan){
      return fCharge;
    } else {
      return chan->CalEnergy(fCharge, fTimestamp);
    }
  }
}

void TDetectorHit::AddEnergy(double eng) {
  SetEnergy(eng + GetEnergy());
}

double TDetectorHit::GetTime() const {
  TChannel* chan = TChannel::GetChannel(fAddress);
  if(!chan){
    return Time() + gRandom->Uniform();
  }
  return chan->CalTime(Time(), fTimestamp);
}

Int_t TDetectorHit::Compare(const TObject *obj) const {
  const TDetectorHit* other = (const TDetectorHit*)obj;
  if(GetEnergy() < other->GetEnergy()) {
    return -1;
  } else if (GetEnergy() > other->GetEnergy()) {
    return +1;
  } else {
    return 0;
  }
}

const char* TDetectorHit::GetName() const {
  TChannel* chan = TChannel::GetChannel(fAddress);
  if(chan){
    return chan->GetName();
  } else {
    return "";
  }
}

void TDetectorHit::Streamer(TBuffer &r_b) {
  if(r_b.IsReading()) {
    r_b.ReadClassBuffer(TDetectorHit::Class(),this);
  } else {
    r_b.WriteClassBuffer(TDetectorHit::Class(),this);
  }
}

int   TDetectorHit::GetNumber() const {
  TChannel* chan = TChannel::GetChannel(fAddress);
  if(chan){
    return chan->GetNumber();
  } else {
    return -1;
  }
}

const char* TDetectorHit::GetInfo()   const {
  TChannel* chan = TChannel::GetChannel(fAddress);
  if(chan){
    return chan->GetInfo();
  } else {
    return "";
  }
}

const char* TDetectorHit::GetSystem() const {
  TChannel* chan = TChannel::GetChannel(fAddress);
  if(chan){
    return chan->GetSystem();
  } else {
    return "";
  }
}
