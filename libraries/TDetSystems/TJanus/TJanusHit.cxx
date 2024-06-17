#include "TJanusHit.h"

#include <mutex>

#include "GValue.h"
#include "TJanus.h"


TJanusHit::TJanusHit(const TJanusHit& hit) {
  hit.Copy(*this);
}

TJanusHit& TJanusHit::operator=(const TJanusHit& hit) {
  hit.Copy(*this);
  return *this;
}

void TJanusHit::Copy(TObject& obj) const {
  TDetectorHit::Copy(obj);
  TJanusHit& janus = (TJanusHit&)obj;
  janus.fRing      = fRing;
  janus.fSector    = fSector;
}

void TJanusHit::Clear(Option_t* opt) {
  TDetectorHit::Clear(opt);
}


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


TVector3 TJanusHit::GetPosition(bool apply_array_offset) const {
  TVector3 output = TJanus::GetPosition(GetDetnum(), GetRing(), GetSector());
  if(apply_array_offset) {
    output += TVector3(GValue::Value("Janus_X_offset"),
                       GValue::Value("Janus_Y_offset"),
                       GValue::Value("Janus_Z_offset"));
  }
  return output;
}

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

double TJanusHit::GetDefaultDistance() const {
  if(GetDetnum() == 0) return -3.2;
  else return 2.8;
}
