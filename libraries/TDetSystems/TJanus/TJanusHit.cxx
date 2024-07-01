#include "TJanusHit.h"

#include <mutex>

#include "GValue.h"
#include "TJanus.h"


/*******************************************************************************/
/* Copies hit ******************************************************************/
/*******************************************************************************/
TJanusHit::TJanusHit(const TJanusHit& hit) {
  hit.Copy(*this);
}

TJanusHit& TJanusHit::operator=(const TJanusHit& hit) {
  hit.Copy(*this);
  return *this;
}

void TJanusHit::Copy(TObject& obj) const {
  TDetectorHit::Copy(obj);
  TJanusHit& janus  = (TJanusHit&)obj;
  janus.fRing       = fRing;
  janus.fSector     = fSector;
  janus.fBackCharge = fBackCharge;
}

/*******************************************************************************/
/* Clear hit *******************************************************************/
/*******************************************************************************/
void TJanusHit::Clear(Option_t* opt) {
  TDetectorHit::Clear(opt);
}

/*******************************************************************************/
/* Empty Print function ********************************************************/
/*******************************************************************************/
void TJanusHit::Print(Option_t *opt) const {
  std::cout << "Some Function" << std::endl;
}

/*******************************************************************************/
/* Returns detector number based on channels.cal file definition ***************/
/*******************************************************************************/
int TJanusHit::GetDetnum() const {
  TChannel* chan = TChannel::GetChannel(fAddress);
  int output = -1;
  if(chan && fAddress !=- 1){
    output = chan->GetArrayPosition();
  } else {
    // std::cout << "Unknown address: " << std::hex << fAddress << std::dec << std::endl;
    output = -1;
  }
  return output;
}

/*******************************************************************************/
/* Returns Hit Position based on ring/sector number ****************************/
/* If smear is selected the position is smeared to get a uniform distribuition */
/* apply_array_offset will use offset values defined in values file ************/
/*******************************************************************************/
TVector3 TJanusHit::GetPosition(bool smear, bool apply_array_offset) const {
  bool secdown = false;
  if(GetDetnum() == 0) secdown = true;
  TVector3 output = TJanus::GetPosition(GetRing() - 1, GetSector() - 1, GetDefaultDistance(), secdown, smear);
  if(apply_array_offset) {
    if(std::isnan(GValue::Value("Janus_X_offset")) && std::isnan(GValue::Value("Janus_Y_offset")) && std::isnan(GValue::Value("Janus_Z_offset"))) {
      output += TVector3(GValue::Value("Janus_X_offset"), GValue::Value("Janus_Y_offset"), GValue::Value("Janus_Z_offset"));
    }
  }
  return output;
}

/*******************************************************************************/
/* Return Detector Z position **************************************************/
/*******************************************************************************/
double TJanusHit::GetDefaultDistance() const {
  if(GetDetnum() == 0) return -3.2;
  else return 2.8;
}

/*******************************************************************************/
/* Uncalibrated Energies - ignores overflows by making them < 0 ****************/
/*******************************************************************************/
Int_t TJanusHit::Charge() const {
  if(fCharge > 60000) {
    return fCharge - 70000;
  } else {
    return fCharge;
  }
}

/*******************************************************************************/
/* Uncalibrated "Back Energy" only valid if hit has ring + sector information **/
/*******************************************************************************/
Double_t TJanusHit::BackCharge() const {
  if(fBackCharge > 60000) {
    return fBackCharge - 70000;
  } else {
    return fBackCharge;
  }
}
