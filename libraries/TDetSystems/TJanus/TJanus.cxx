#include "TJanus.h"

#include <cassert>
#include <iostream>
#include <memory>

#include "TMath.h"
#include "TRandom.h"
#include "DDASDataFormat.h"
#include "TRawEvent.h"

int TJanus::NRing = 24;
int TJanus::NSector = 32;

double TJanus::PhiOffset = 1.5*TMath::Pi();
double TJanus::OuterDiameter  = 7.0;
double TJanus::InnerDiameter  = 2.2;
double TJanus::TargetDistance = 3.1;

double TJanus::TDiff = 900;
double TJanus::EWin = 0.9;
double TJanus::FrontBackOffset = 0;

bool TJanus::multhit = false;

TJanus::TJanus() {
  Clear();
}

TJanus::~TJanus(){ }

void TJanus::Copy(TObject& obj) const {
  TDetector::Copy(obj);

  TJanus& janus = (TJanus&)obj;
  janus.janus_hits = janus_hits;
  janus.ring_hits = ring_hits;
  janus.sector_hits = sector_hits;
}

void TJanus::Clear(Option_t* opt){
  TDetector::Clear(opt);

  janus_hits.clear();
  ring_hits.clear();
  sector_hits.clear();
}

int TJanus::BuildHits(std::vector<TRawEvent>& raw_data){
  long int smallest_timestamp = 0x7fffffffffffffff;
  for(auto& event : raw_data){
    //Unpack raw data into a DDAS Event
    TSmartBuffer buf = event.GetPayloadBuffer();
    TDDASEvent<DDASHeader> ddas(buf);
    unsigned int address = ( (5<<24) + (ddas.GetCrateID()<<16) + (ddas.GetSlotID()<<8) + ddas.GetChannelID() );
    //If channel not found in channels.cal file skip and do nothing
    TChannel* chan = TChannel::GetChannel(address);
    static int lines_displayed = 0;
    if(!chan){
      if(lines_displayed < 10 && ddas.GetCrateID()) {
        std::cout << "Unknown Janus (crate, slot, channel): (" << ddas.GetCrateID() << ", " << ddas.GetSlotID() << ", " <<  ddas.GetChannelID() << endl;
      }
      lines_displayed++;
      continue;
    }
    TJanusHit hit;
    hit.SetAddress(address);
    hit.SetTimestamp(ddas.GetTimestamp()); // Timestamp in ns
    if(hit.Timestamp()<smallest_timestamp) { smallest_timestamp = hit.Timestamp(); }
    hit.SetTime(ddas.GetTime()); // Timestamp + CFD
    hit.SetCFDTime(ddas.GetCFDTime()); // CFD ONLY
    hit.SetCharge(ddas.GetEnergy());
    if(*chan->GetArraySubposition() == 'F'){
      hit.SetRingNumber(chan->GetSegment());
      hit.SetSectorNumber(-1);
      ring_hits.push_back(std::move(hit));
    } else {
      hit.SetRingNumber(-1);
      hit.SetSectorNumber(chan->GetSegment());
      sector_hits.push_back(std::move(hit));
    }
  }
  return ring_hits.size() + sector_hits.size();
}

TDetectorHit& TJanus::GetHit(int i){
  return janus_hits.at(i);
}


TVector3 TJanus::GetPosition(int detnum, int ring_num, int sector_num){
  if(detnum<0 || detnum>1 ||
     ring_num<1 || ring_num>24 ||
     sector_num<1 || sector_num>32){
    return TVector3(std::sqrt(-1),std::sqrt(-1),std::sqrt(-1));
  }

  TVector3 origin = TVector3(0,0,3);
  // Phi of sector 1 of downstream detector
  double phi_offset = 2*3.1415926535*(0.25);

  // Winding direction of sectors.
  bool clockwise = true;

  double janus_outer_radius = 3.5;
  double janus_inner_radius = 1.1;

  TVector3 position(1,0,0);  // Not (0,0,0), because otherwise SetPerp and SetPhi fail.
  double rad_slope = (janus_outer_radius - janus_inner_radius) /24;
  double rad_offset = janus_inner_radius;
  double perp_num = ring_num - 0.5; // Shift to 0.5-23.5
  position.SetPerp(perp_num*rad_slope + rad_offset);
  double phi_num = sector_num;
  double phi =
    phi_offset +
    (clockwise ? -1 : 1) * 2*3.1415926/32 * (phi_num - 1);
  position.SetPhi(phi);

  position += origin;

  if(detnum==0){
    position.RotateY(TMath::Pi());
  }

  return position;
}

Int_t TJanus::GetJanusSize() {
  BuildJanusHit();
  return janus_hits.size();
}

void TJanus::BuildJanusHit() {
  //Constructs front/back coincidences, energy and time differences can be changed with
  //SetFrontBackEn and SetFrontBackTime
  janus_hits.clear();
  //std::cout << "Begin Hit" << std::endl;
  if(ring_hits.size() == 0 || sector_hits.size() == 0) return;

  //For storage of energies
  std::vector<double> EnR, EnS;
  std::vector<bool> UR, US;

  for(size_t i = 0; i < ring_hits.size(); i++) {
    EnR.push_back(ring_hits.at(i).GetEnergy());
    UR.push_back(false);
  }

  for(size_t i = 0; i < sector_hits.size(); i++) {
    EnS.push_back(sector_hits.at(i).GetEnergy());
    US.push_back(false);
  }


  for(size_t i = 0; i < ring_hits.size(); i++) {
    for(size_t j = 0; j < sector_hits.size(); j++) {
      if(ring_hits.at(i).GetDetnum() == sector_hits.at(j).GetDetnum()) {
	//Check time between events is good
	if(abs(ring_hits.at(i).Timestamp() - sector_hits.at(j).Timestamp()) < TDiff) {
	  //Check energy is good
	  if((EnR.at(i)*EWin) < EnS.at(j) && (EnS.at(j)*EWin) < EnR.at(i)) {
	    //Use rings for all data, sector for position only
	    TJanusHit dhit = ring_hits.at(i);
	    dhit.SetSectorNumber(sector_hits.at(j).GetSector());
	    janus_hits.push_back(dhit);
	    //For multihit events
	    UR.at(i) = true;
	    US.at(j) = true;
	  }
	}
      }
    }
  }
  if(multhit) {
    int ringCount = 0;
    int secCount = 0;
    //Check for hits without segments, charge may be shared between multiple strips
    for(size_t i = 0; i < UR.size(); i++) {
      if(!UR.at(i)) ringCount++;
    }
    for(size_t i = 0; i < US.size(); i++) {
      if(!US.at(i)) secCount++;
    }
    // Check for single ring hit + multiple sectors
    for(size_t i = 0; i < UR.size(); i++) {
      if(UR.at(i)) continue;
      for(size_t j = 0; j < US.size(); j++) {
        if(US.at(j)) continue;
        if(ring_hits.at(i).GetDetnum() != sector_hits.at(j).GetDetnum()) continue;
        for(size_t k = j+1; k < US.size(); k++) {
          if(US.at(k)) continue;
          if(ring_hits.at(i).GetDetnum() != sector_hits.at(k).GetDetnum()) continue;
	  //Check time of hits
  	  if(abs(ring_hits.at(i).Timestamp() - sector_hits.at(j).Timestamp()) < TDiff &&
             abs(ring_hits.at(i).Timestamp() - sector_hits.at(k).Timestamp()) < TDiff) {
	    //check ring energy is equal to sum of sector energy
	    if((EnR.at(i)*EWin) < EnS.at(j)+EnS.at(k) && ((EnS.at(j)+EnS.at(k))*EWin) < EnR.at(i)) {
	      int dSec = abs(sector_hits.at(j).GetSector() - sector_hits.at(k).GetSector());
	      if( dSec == 1 || dSec == NSector) {  //Check for neighbouring Sectors
		//Likely a charge sharing event, create a new hit
		TJanusHit dhit = ring_hits.at(i);
		//Assign sector number as highest energy strip
	        if(sector_hits.at(j) < sector_hits.at(k)) {
		  dhit.SetSectorNumber(sector_hits.at(k).GetSector());
		} else {
                  dhit.SetSectorNumber(sector_hits.at(j).GetSector());
		}
		janus_hits.push_back(dhit);
	      } else {
		//Two different hits with the same ring
		//Have to use sector information for hit, ring gives position only
		TJanusHit dhit0 = sector_hits.at(j);
		dhit0.SetRingNumber(ring_hits.at(i).GetRing());
		janus_hits.push_back(dhit0);

		TJanusHit dhit1 = sector_hits.at(k);
		dhit1.SetRingNumber(ring_hits.at(i).GetRing());
		janus_hits.push_back(dhit1);
	      } //Sector difference
	      UR.at(i) = true;
	      US.at(j) = true;
	      US.at(k) = true;
	    } //Energy
	  } //Time
	} //2nd sector
      } //1st sector
    }  //Shared Ring


    ringCount = 0;
    secCount = 0;
    //Check for hits without segments, charge may be shared between multiple strips
    for(size_t i = 0; i < UR.size(); i++) {
      if(!UR.at(i)) ringCount++;
    }
    for(size_t i = 0; i < US.size(); i++) {
      if(!US.at(i)) secCount++;
    }
    // Check for single sector hit + multiple rings
    for(size_t i = 0; i < US.size(); i++) {
      if(US.at(i)) continue;
      for(size_t j = 0; j < UR.size(); j++) {
        if(UR.at(j)) continue;
        if(sector_hits.at(i).GetDetnum() != ring_hits.at(j).GetDetnum()) continue;
        for(size_t k = j+1; k < UR.size(); k++) {
          if(UR.at(k)) continue;
          if(sector_hits.at(i).GetDetnum() != ring_hits.at(k).GetDetnum()) continue;
	  //Check time of hits
  	  if(abs(sector_hits.at(i).Timestamp() - ring_hits.at(j).Timestamp()) < TDiff &&
             abs(sector_hits.at(i).Timestamp() - ring_hits.at(k).Timestamp()) < TDiff) {
	    //check sector energy is equal to sum of ring energy
	    if((EnS.at(i)*EWin) < EnR.at(j)+EnR.at(k) && ((EnR.at(j)+EnR.at(k))*EWin) < EnS.at(i)) {
	      int dRing = abs(ring_hits.at(j).GetRing() - ring_hits.at(k).GetRing());
	      if( dRing == 1 || dRing == NRing) {  //Check for neighbouring rings
		//Likely a charge sharing event, create a new hit
		//Use sector information for hit, ring gives position only
		TJanusHit dhit = sector_hits.at(i);
		//Assign ring number as highest energy strip
	        if(ring_hits.at(j) < ring_hits.at(k)) {
		  dhit.SetRingNumber(ring_hits.at(k).GetRing());
		} else {
                  dhit.SetRingNumber(ring_hits.at(j).GetRing());
		}
		janus_hits.push_back(dhit);
	      } else {
		//Two different hits with the same sector
		//Use ring information for hit, sector gives position only
		TJanusHit dhit0 = ring_hits.at(j);
		dhit0.SetSectorNumber(sector_hits.at(i).GetSector());
		janus_hits.push_back(dhit0);

		TJanusHit dhit1 = ring_hits.at(k);
		dhit1.SetSectorNumber(sector_hits.at(i).GetRing());
		janus_hits.push_back(dhit1);
	      } //Ring difference
	      US.at(i) = true;
	      UR.at(j) = true;
	      UR.at(k) = true;
	    } //Energy
	  } //Time
	} //2nd ring
      } //1st ring
    }  //Shared sector
  }  //Multihit
}

TVector3 TJanus::GetPosition(int ring, int sector, double zoffset, bool sectorsdownstream, bool smear) {
  double ring_width = (OuterDiameter - InnerDiameter) * 0.5 / NRing;
  double inner_radius = (InnerDiameter)/2;
  double radius = inner_radius + ring_width * (ring + 0.5);

  double phi_width = 2 * TMath::Pi() / NSector;
  double phi = phi_width * sector;

  if(sectorsdownstream) {
    phi = -phi;
  }
  phi += PhiOffset;

  //Uniform distrubution over S3 pixel
  if(smear) {
    double sep = ring_width * 0.025;
    double r1 = radius - ring_width * 0.5 + sep;
    double r2 = radius + ring_width * 0.5 - sep;
    radius = sqrt(gRandom->Uniform(r1 * r1, r2 * r2));
    double sepphi = sep / radius;
    phi = gRandom->Uniform(phi - phi_width * 0.5 + sepphi, phi + phi_width * 0.5 - sepphi);
  }
  return TVector3(cos(phi) * radius, sin(phi) * radius, zoffset);
}

void TJanus::InsertHit(const TDetectorHit& hit) {
  janus_hits.emplace_back((TJanusHit&)hit);
  fSize++;
}

TJanusHit& TJanus::GetRingHit(int i) {
  return ring_hits.at(i);
}

TJanusHit& TJanus::GetSectorHit(int i) {
  return sector_hits.at(i);
}

TJanusHit& TJanus::GetJanusHit(int i) {
  return janus_hits.at(i);
}

void TJanus::Print(Option_t *opt) const {
  std::cout << "TJanus @ " << Timestamp() << std::endl;
  std::cout << "Size: " << Size() << std::endl;
  for(unsigned int i=0;i<Size();i++) {
    janus_hits.at(i).Print();
  }
  std::cout << "---------------------------" << std::endl;

}
