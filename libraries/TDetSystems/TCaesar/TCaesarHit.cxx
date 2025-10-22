#include "TCaesarHit.h"
#include "TCaesar.h"
#include "GValue.h"

#include <algorithm>
#include <iostream>
#include "TString.h"


TCaesarHit::TCaesarHit() {
  Clear();
}

/*******************************************************************************/
/* Copies hit ******************************************************************/
/*******************************************************************************/
void TCaesarHit::Copy(TObject& obj) const{
  TDetectorHit::Copy(obj);
  TCaesarHit& caesar_hits = (TCaesarHit&)obj;
  caesar_hits.fChan = fChan;
  caesar_hits.fTimes = fTimes;
  caesar_hits.fOverflow = fOverflow;
  caesar_hits.fUnderflow = fUnderflow;
}

/*******************************************************************************/
/* Clear hit *******************************************************************/
/*******************************************************************************/
void TCaesarHit::Clear(Option_t *opt) {
  TDetectorHit::Clear(opt);
  fChan = - 1;
  fAB = - 1;
  fOverflow = 0;
  fUnderflow = 0;
  fTimes.clear();
}

/*******************************************************************************/
/* Basic Print Function ********************************************************/
/*******************************************************************************/
void TCaesarHit::Print(Option_t *opt) const {
  std::cout << "TCaesarHit:\n" << "\tChannel: " << GetChannel() << "\n" << "\tCharge: " << Charge() << "\n" << std::flush;
}

/*******************************************************************************/
/* Returns position of Detector ************************************************/
/*******************************************************************************/
TVector3 TCaesarHit::GetPosition(bool apply_array_offset, TVector3 array_offset) const {
  TVector3 array_pos = TCaesar::GetPosition(GetChannel());
  if(apply_array_offset){
    if(std::isnan(array_offset.X()) && std::isnan(array_offset.Y()) && std::isnan(array_offset.Z())) {
      if(!std::isnan(GValue::Value("CAESAR_X_offset")) && !std::isnan(GValue::Value("CAESAR_Y_offset")) && !std::isnan(GValue::Value("CAESAR_Z_offset"))) {
        array_offset = TVector3(GValue::Value("CAESAR_X_offset"), GValue::Value("CAESAR_Y_offset"), GValue::Value("CAESAR_Z_offset"));
      } else array_offset = TVector3(0,0,0);
    }
    array_pos += array_offset;
  }
  return array_pos;
}

/*******************************************************************************/
/* Returns ring of detector ****************************************************/
/*******************************************************************************/
int TCaesarHit::GetRing() const {
  return TCaesar::GetRing(GetChannel());
}

/*******************************************************************************/
/* Returns detector number within ring of detector *****************************/
/*******************************************************************************/
int TCaesarHit::GetDetectorNumber() const {
  return TCaesar::GetDetectorNumber(GetChannel());
}

/*******************************************************************************/
/* Used by Nearest Neighbour Addback code to add energies **********************/
/*******************************************************************************/
void TCaesarHit::NNAdd(const TCaesarHit& rhs) {
  double energy = GetEnergy() + rhs.GetEnergy();
  SetEnergy(energy, true);
}

/*******************************************************************************/
/* Uncalibrated DDAS energies **************************************************/
/*******************************************************************************/
Int_t TCaesarHit::Charge() const {
  if(fCharge > 30000) {
    return fCharge - 32768;
  } else {
    return fCharge;
  }
}
