#include "TSunHit.h"
#include "GValue.h"

#include <algorithm>
#include <iostream>
#include "TString.h"


TSunHit::TSunHit() {
  Clear();
}

/*******************************************************************************/
/* Copies hit ******************************************************************/
/*******************************************************************************/
void TSunHit::Copy(TObject& obj) const{
  TDetectorHit::Copy(obj);
  TSunHit& sun = (TSunHit&)obj;
  sun.fTrace = fTrace;
}

/*******************************************************************************/
/* Clear hit *******************************************************************/
/*******************************************************************************/
void TSunHit::Clear(Option_t *opt) {
  TDetectorHit::Clear(opt);
}

/*******************************************************************************/
/* Basic Print Function ********************************************************/
/*******************************************************************************/
void TSunHit::Print(Option_t *opt) const {
  std::cout << "TSunHit:\n" << "\tChannel: " << GetChannel() << "\n" << "\tCharge: " << Charge() << "\n" << std::flush;
}

/*******************************************************************************/
/* Sets trace information if present in data ***********************************/
/*******************************************************************************/
void TSunHit::SetTrace(unsigned int trace_length, const unsigned short* trace) {
  fTrace.clear();
  fTrace.reserve(trace_length);
  copy(trace,trace+trace_length,back_inserter(fTrace));
}

/*******************************************************************************/
/* Returns DDAS crate/slot channel number **************************************/
/*******************************************************************************/
int TSunHit::GetCrate() const {
  return (fAddress&0x00ff0000)>>16;
}

int TSunHit::GetSlot() const {
  return (fAddress&0x0000ff00)>>8;
}

int TSunHit::GetChannel() const {
  return (fAddress&0x000000ff)>>0;
}

/*******************************************************************************/
/* Uncalibrated DDAS energies **************************************************/
/*******************************************************************************/
Int_t TSunHit::Charge() const {
  if(fCharge > 30000) {
    return fCharge - 32768;
  } else {
    return fCharge;
  }
}
