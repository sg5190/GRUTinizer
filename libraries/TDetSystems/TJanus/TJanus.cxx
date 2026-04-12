 #include "TJanus.h"

#include <cassert>
#include <iostream>
#include <memory>

#include "TMath.h"
#include "TRandom.h"
#include "DDASDataFormat.h"
#include "TRawEvent.h"

/*******************************************************************************/
/* Geometric definitions used to calculate hit positions ***********************/
/*******************************************************************************/
int TJanus::NRing = 24;
int TJanus::NSector = 32;
double TJanus::PhiOffset = 1.5*TMath::Pi();
double TJanus::OuterDiameter  = 70.0;
double TJanus::InnerDiameter  = 22.0;
double TJanus::TargetDistance = 31.0;

/*******************************************************************************/
/* Default Energy/Time conditions and booleans used to build a TJanusHit *******/
/*******************************************************************************/
double TJanus::DThresh = 400;
double TJanus::UThresh = 400;
double TJanus::TDiff = 100;
double TJanus::EWin = 0.9;
bool TJanus::Addback = false;
bool TJanus::Multihit = false;
bool TJanus::EBuild = false;
bool TJanus::Diagnostics = false;
/*******************************************************************************/
/* TJanus **********************************************************************/
/* Main class for unpacking JANUS (2 S3 Detectors) in DDAS electronics *********/
/*******************************************************************************/
TJanus::TJanus() {
  Clear();
}

TJanus::~TJanus(){ }

/*******************************************************************************/
/* Copies hit ******************************************************************/
/*******************************************************************************/
void TJanus::Copy(TObject& obj) const {
  TDetector::Copy(obj);

  TJanus& janus = (TJanus&)obj;
  janus.janus_hits = janus_hits;
  janus.ring_hits = ring_hits;
  janus.sector_hits = sector_hits;
  janus.fMismatch = fMismatch;
}

/*******************************************************************************/
/* Clear hit *******************************************************************/
/*******************************************************************************/
void TJanus::Clear(Option_t* opt){
  TDetector::Clear(opt);

  janus_hits.clear();
  ring_hits.clear();
  sector_hits.clear();
}

/*******************************************************************************/
/* Unpacks raw DDAS data into ring/sector hit vectors **************************/
/*******************************************************************************/
int TJanus::BuildHits(std::vector<TRawEvent>& raw_data){
  long int smallest_timestamp = 0x7fffffffffffffff;
  for(auto& event : raw_data){
    //Get raw data and unpacks it as a DDAS Event (DDASDataFormat.h)
    TSmartBuffer buf = event.GetPayloadBuffer();
    TDDASEvent<DDASHeader> ddas(buf);
    unsigned int address = ( (5<<24) + (ddas.GetCrateID()<<16) + (ddas.GetSlotID()<<8) + ddas.GetChannelID() );
    //If channel not found in calibration file (*.cal) file skip and do nothing
    TChannel* chan = TChannel::GetChannel(address);
    static int lines_displayed = 0;
    if(!chan){
      if(lines_displayed < 10 && ddas.GetCrateID()) {
        std::cout << "Unknown Janus (crate, slot, channel): (" << ddas.GetCrateID() << ", " << ddas.GetSlotID() << ", " <<  ddas.GetChannelID() << endl;
      }
      lines_displayed++;
      continue;
    }
    //This is not strictly events with the overflow flag, but anything
    //over channel 65000 which is never good anyway so don't bother making a hit
    if(ddas.IsOverflow()) continue;
    //Make a new TJanusHit and set values
    TJanusHit * hit = NULL;
    //Check if the hit is a front hit (ring) or back hit (sector)
    if(*chan->GetArraySubposition() == 'F'){
      ring_hits.emplace_back();
      hit = &ring_hits.back();
      hit->SetRingNumber(chan->GetSegment());
      hit->SetSectorNumber(-1);
    } else {
      sector_hits.emplace_back();
      hit = &sector_hits.back();
      hit->SetRingNumber(-1);
      hit->SetSectorNumber(chan->GetSegment());
    }
    hit->SetAddress(address);
    hit->SetTimestamp(ddas.GetTimestamp(), ddas.GetTsUnit());
    if(hit->GetTimestamp()<smallest_timestamp) { smallest_timestamp = hit->GetTimestamp(); }
    hit->SetTime(ddas.GetTime()); // Timestamp + CFD
    hit->SetCFDTime(ddas.GetCFDTime()); // CFD ONLY
    hit->SetCharge(ddas.GetEnergy());
  }

  return ring_hits.size() + sector_hits.size();
}

/*******************************************************************************/
/* Functions to call JANUS hits ************************************************/
/*******************************************************************************/
TDetectorHit& TJanus::GetHit(int i){
  return janus_hits.at(i);
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

/*******************************************************************************/
/* When function is called builds a Janus hit from the Ring/sector information */
/* and returns size of the newly created hit vector ****************************/
/*******************************************************************************/
Int_t TJanus::GetJanusSize() {
  BuildJanusHit();
  return janus_hits.size();
}

/*******************************************************************************/
/* Main function to build a JANUS hits from ring and sector hit vectors ********/
/* Hits are built according to Energy/Time conditions which can be changed *****/
/* The Values which can be adjusted are ****************************************/
/* DThesh   - Set by TJanus::SetDownChargeThresh(double) - Exclude downstream **/
/*            ring/sector hits with charge less than DThesh ********************/
/* UThesh   - Set by TJanus::SetUpChargeThresh(double) - Exclude upstream ******/
/*            ring/sector hits with charge less than UThesh ********************/
/* TDiff    - Set by SetFrontBackTime(int) - janus hits will be built only if **/
/*            Ring/Sector occur within TDiff (in ns) of each other *************/
/* EDiff    - Set by SetFrontBackEnergy(double) - janus hits will be built *****/
/* 	      only if Ring/Sector Charge/Energy ratio > EDiff of each other ****/
/* EBuild   - Set by SetEnergyBuild(bool) - If true uses calibrated energies ***/
/*	      to make janus hits otherwise uses raw charge *********************/
/* Addback  - Set by SetAddback(bool) - If true loop over ring/sector hits to **/
/*            check for charge/sharing between neighbouring ring/sectors and ***/
/*	      perform a simple addback before trying to create a janus hit *****/
/* MultiHit - Set by SetMultiHit(bool) - If true will consider cases where *****/
/*            you have two hits with the same ring and different sectors and ***/
/*            vice versa, this is typically a small correction *****************/
/*******************************************************************************/
void TJanus::BuildJanusHit() {
  //Don't rebuild hits, hits should be cleared with ResetJanus()
  if(janus_hits.size() > 0) return;
  //If either rings or sectors are empty do nothing
  if(ring_hits.size() == 0 || sector_hits.size() == 0) return;

  //Define vector to have easier access to ring/sector energy/charges
  std::vector<double> EnR, EnS;
  EnR.resize(ring_hits.size(), 0.0);
  EnS.resize(sector_hits.size(), 0.0);

  //Define boolean vectors used to ensure the same ring/sector cannot be used
  //to make multiple hits
  std::vector<bool> UR, US;
  UR.resize(ring_hits.size(), false);
  US.resize(sector_hits.size(), false);

  //Define diagnostic array
  int diagnostics[2][9] = {0};

  //Exclude ring hits if the charge is below the charge threshold performs
  //addback if selected by SetAddback()
  for(size_t i = 0; i < ring_hits.size(); i++) {
    diagnostics[ring_hits.at(i).GetDetnum()%2][0]++;
    if( (ring_hits.at(i).GetCharge() < DThresh  && ring_hits.at(i).IsDownstream()) ||
       (ring_hits.at(i).GetCharge() < UThresh  && !ring_hits.at(i).IsDownstream())) {
      UR.at(i) = true;
      continue;
    }
    diagnostics[ring_hits.at(i).GetDetnum()%2][2]++;
    if(UR.at(i)) continue;
    if(Addback) {
      for(size_t j = i + 1; j < ring_hits.size(); j++) {
	if(UR.at(j)) continue;
        if((ring_hits.at(j).GetCharge() < DThresh) && ring_hits.at(j).IsDownstream()) continue;
        if((ring_hits.at(j).GetCharge() < UThresh) && !ring_hits.at(j).IsDownstream()) continue;
	if((ring_hits.at(i).GetDetnum() == ring_hits.at(j).GetDetnum()) && (abs(ring_hits.at(i).GetTimestampNs() - ring_hits.at(j).GetTimestampNs()) < TDiff)) {
	  int dRing = abs(ring_hits.at(i).GetRing() - ring_hits.at(j).GetRing());
	  if(dRing == 1) {
	    //Sum energy/charge and assign the highest energy hit as the ring number
	    if(EBuild) {
	      ring_hits.at(i).SetEnergy(ring_hits.at(i).GetEnergy() + ring_hits.at(j).GetEnergy(), true);
	      if(ring_hits.at(i).GetEnergy() < ring_hits.at(j).GetEnergy()) ring_hits.at(i).SetRingNumber(ring_hits.at(j).GetRing());
	    }
	    else {
	      ring_hits.at(i).SetCharge(ring_hits.at(i).GetCharge() + ring_hits.at(j).GetCharge());
	      if(ring_hits.at(i).GetCharge() < ring_hits.at(j).GetCharge()) ring_hits.at(i).SetRingNumber(ring_hits.at(j).GetRing());
	    }
	    UR.at(j) = true;
	    diagnostics[ring_hits.at(i).GetDetnum()%2][4]++;
	  }
	}
      }
    }
    if(EBuild) {
      EnR.at(i) = ring_hits.at(i).GetEnergy();
    } else {
      EnR.at(i) = ring_hits.at(i).GetCharge();
    }
  }

  //Exclude sectors hits if the charge is below the charge threshold performs
  //addback if selected by SetAddback()
  for(size_t i = 0; i < sector_hits.size(); i++) {
    diagnostics[sector_hits.at(i).GetDetnum()%2][1]++;
    if( (sector_hits.at(i).GetCharge() < DThresh  && sector_hits.at(i).IsDownstream()) ||
        (sector_hits.at(i).GetCharge() < UThresh  && !sector_hits.at(i).IsDownstream())) {
      US.at(i) = true;
      continue;
    }
    diagnostics[sector_hits.at(i).GetDetnum()%2][3]++;
    if(US.at(i)) continue;
    if(Addback) {
      for(size_t j = i + 1; j < sector_hits.size(); j++) {
	if(US.at(j)) continue;
        if(sector_hits.at(j).GetCharge() < DThresh  && sector_hits.at(j).IsDownstream()) continue;
        if(sector_hits.at(j).GetCharge() < UThresh  && !sector_hits.at(j).IsDownstream()) continue;
	if((sector_hits.at(i).GetDetnum() == sector_hits.at(j).GetDetnum()) && (abs(sector_hits.at(i).GetTimestampNs() - sector_hits.at(j).GetTimestampNs()) < TDiff)) {
	  int dSector = abs(sector_hits.at(i).GetSector() - sector_hits.at(j).GetSector());
          //Sum energy/charge and assign the highest energy hit as the sector number
	  if(dSector == 1 || dSector == NSector) {
            if(EBuild) {
              if(sector_hits.at(i).GetEnergy() < sector_hits.at(j).GetEnergy()) sector_hits.at(i).SetSectorNumber(sector_hits.at(j).GetSector());
 	      sector_hits.at(i).SetEnergy(sector_hits.at(i).GetEnergy() + sector_hits.at(j).GetEnergy(), true);
            }
            else {
              if(sector_hits.at(i).GetCharge() < sector_hits.at(j).GetCharge()) sector_hits.at(i).SetSectorNumber(sector_hits.at(j).GetSector());
  	      sector_hits.at(i).SetCharge(sector_hits.at(i).GetCharge() + sector_hits.at(j).GetCharge());
            }
	    US.at(j) = true;
    	    diagnostics[sector_hits.at(i).GetDetnum()%2][5]++;
	  }
	}
      }
    }
    if(EBuild) {
      EnS.at(i) = sector_hits.at(i).GetEnergy();
    } else {
      EnS.at(i) = sector_hits.at(i).GetCharge();
    }
  }

  //Main part of function which creates a hit
  for(size_t i = 0; i < ring_hits.size(); i++) {
    if(UR.at(i)) continue;
    for(size_t j = 0; j < sector_hits.size(); j++) {
      if(UR.at(i) || US.at(j)) continue;
      if(ring_hits.at(i).GetDetnum() == sector_hits.at(j).GetDetnum()) {
	//Check time between events is good
	if(abs(ring_hits.at(i).GetTimestampNs() - sector_hits.at(j).GetTimestampNs()) < TDiff) {
	  //Check energy is good
	  if((EnR.at(i)*EWin) < EnS.at(j) && (EnS.at(j)*EWin) < EnR.at(i)) {
	    //Use rings for all data, sector for position only
	    TJanusHit dhit = ring_hits.at(i);
	    dhit.SetSectorNumber(sector_hits.at(j).GetSector());
	    dhit.SetBackCharge(EnS.at(j));
	    janus_hits.push_back(dhit);
	    //Now that hit is created do not consider this ring/sector again
	    UR.at(i) = true;
	    US.at(j) = true;
    	    diagnostics[ring_hits.at(i).GetDetnum()%2][6]++;
	  }
	}
      }
    }
  }
  //Consider cases where two hits have the same ring but different sectors or
  //the same sector and different rings
  if(Multihit) {
    //Check for cases with a single ring hit + multiple sectors
    for(size_t i = 0; i < UR.size(); i++) {
      if(UR.at(i)) continue;
      for(size_t j = 0; j < US.size(); j++) {
        if(UR.at(i) || US.at(j)) continue;
        if(ring_hits.at(i).GetDetnum() != sector_hits.at(j).GetDetnum()) continue;
        for(size_t k = j+1; k < US.size(); k++) {
          if(UR.at(i) || US.at(k)) continue;
          if(ring_hits.at(i).GetDetnum() != sector_hits.at(k).GetDetnum()) continue;
	  //Check all hits are within TDiff of each other
  	  if(abs(ring_hits.at(i).GetTimestampNs() - sector_hits.at(j).GetTimestampNs()) < TDiff &&
             abs(ring_hits.at(i).GetTimestampNs() - sector_hits.at(k).GetTimestampNs()) < TDiff) {
	    //check ring energy is equal to sum of sector energies
	    if((EnR.at(i)*EWin) < EnS.at(j)+EnS.at(k) && ((EnS.at(j)+EnS.at(k))*EWin) < EnR.at(i)) {
	      int dSec = abs(sector_hits.at(j).GetSector() - sector_hits.at(k).GetSector());
	      //If using Addback this should never be true, if not that you don't want addback anyway so skip
	      if(( dSec == 1 || dSec == NSector)) continue;
	      else {
		//Two different hits within the same ring
		//Make two hits using sector information for hit, ring gives position only
		//For ring (back) charge/energy take it as the ratio of sector charges/energies
		TJanusHit dhit0 = sector_hits.at(j);
		dhit0.SetRingNumber(ring_hits.at(i).GetRing());
     	        dhit0.SetBackCharge(EnR.at(i)*(EnS.at(j)/(EnS.at(j)+EnS.at(k))));
  	        janus_hits.push_back(dhit0);

		TJanusHit dhit1 = sector_hits.at(k);
		dhit1.SetRingNumber(ring_hits.at(i).GetRing());
     	        dhit1.SetBackCharge(EnR.at(i)*(EnS.at(k)/(EnS.at(j)+EnS.at(k))));
		janus_hits.push_back(dhit1);

     	        diagnostics[ring_hits.at(i).GetDetnum()%2][7] += 2;
	      } //Sector difference
	      UR.at(i) = true;
	      US.at(j) = true;
	      US.at(k) = true;
	    } //Energy
	  } //Time
	} //2nd sector
      } //1st sector
    }  //Shared Ring

    // Check for single sector hit + multiple rings
    for(size_t i = 0; i < US.size(); i++) {
      if(US.at(i)) continue;
      for(size_t j = 0; j < UR.size(); j++) {
        if(US.at(i) || UR.at(j)) continue;
        if(sector_hits.at(i).GetDetnum() != ring_hits.at(j).GetDetnum()) continue;
        for(size_t k = j+1; k < UR.size(); k++) {
          if(US.at(i) || UR.at(k)) continue;
          if(sector_hits.at(i).GetDetnum() != ring_hits.at(k).GetDetnum()) continue;
	  //Check time of hits
  	  if(abs(sector_hits.at(i).GetTimestampNs() - ring_hits.at(j).GetTimestampNs()) < TDiff &&
             abs(sector_hits.at(i).GetTimestampNs() - ring_hits.at(k).GetTimestampNs()) < TDiff) {
	    //check sector energy is equal to sum of ring energy
	    if((EnS.at(i)*EWin) < EnR.at(j)+EnR.at(k) && ((EnR.at(j)+EnR.at(k))*EWin) < EnS.at(i)) {
	      int dRing = abs(ring_hits.at(j).GetRing() - ring_hits.at(k).GetRing());
	      //If using Addback this should never be true, if not that you don't want addback anyway so skip
	      if( dRing == 1)  continue;
	      else {
		//Two different hits within the same sector
		//Make two hits using ring information for hit, sector gives position only
		//For sector (back) charge/energy take it as the ratio of ring charges/energies
		TJanusHit dhit0 = ring_hits.at(j);
		dhit0.SetSectorNumber(sector_hits.at(i).GetSector());
     	        dhit0.SetBackCharge(EnS.at(i)*(EnR.at(j)/EnR.at(j)+EnR.at(k)));
		janus_hits.push_back(dhit0);

		TJanusHit dhit1 = ring_hits.at(k);
		dhit1.SetSectorNumber(sector_hits.at(i).GetRing());
     	        dhit1.SetBackCharge(EnS.at(i)*(EnR.at(k)/EnR.at(j)+EnR.at(k)));
		janus_hits.push_back(dhit1);

     	        diagnostics[sector_hits.at(i).GetDetnum()%2][8] += 2;
	      } //Ring difference
  	      //Now that hit is created do not consider this ring/sector again
	      US.at(i) = true;
	      UR.at(j) = true;
	      UR.at(k) = true;
	    } //Energy
	  } //Time
	} //2nd ring
      } //1st ring
    }  //Shared sector
  }  //Multihit

  if(Diagnostics && diagnostics[1][6] > 0) {
    bool mismatch = false;
    string names[9] = {"Rings             ", "Sector            ", "Good Rings        ", "Good Sectors      ", "Addback Rings     ", "Addback Sector    ", "Janus Events      ", "Multihit 2 Sectors", "Multihit 2 Rings"};
    if(diagnostics[1][2] == diagnostics[1][6] && diagnostics[1][3] == diagnostics[1][6] ) mismatch = true;
    if(Addback && diagnostics[1][2]-diagnostics[1][4] == diagnostics[1][6] && diagnostics[1][3]-diagnostics[1][5] == diagnostics[1][6]) mismatch = true;
    if(Multihit && (diagnostics[1][2] == (diagnostics[1][6] + diagnostics[1][7]/2)) && diagnostics[1][3] == (diagnostics[1][6] + diagnostics[1][7])) mismatch = true;
    if(Multihit && (diagnostics[1][2] == (diagnostics[1][6] + diagnostics[1][7])) && (diagnostics[1][3] == (diagnostics[1][6] + diagnostics[1][7]/2))) mismatch = true;
    if(Multihit && Addback && (diagnostics[1][2]-diagnostics[1][4] == (diagnostics[1][6] + diagnostics[1][7]/2)) && (diagnostics[1][3]-diagnostics[1][5]) == (diagnostics[1][6] + diagnostics[1][7])) mismatch = true;
    if(Multihit && Addback && ((diagnostics[1][2]-diagnostics[1][4]) == (diagnostics[1][6] + diagnostics[1][7])) && (diagnostics[1][3]-diagnostics[1][5] == (diagnostics[1][6] + diagnostics[1][7]/2))) mismatch = true;
    fMismatch = mismatch;
    if(!mismatch) {
      std::cout << "Built Janus Event\n\t\tUpstream\tDownstream" << std::endl;
      for(int i = 0; i < 9; i++) std::cout << names[i] << "\t" << diagnostics[0][i] << "\t" << diagnostics[1][i] << endl;
      std::cout << std::endl;
      std::cout << "Ring Energies " << std::endl;
      for(int i = 0; i < (int)ring_hits.size(); i++) if(ring_hits.at(i).GetCharge() > DThresh) std::cout << ring_hits.at(i).GetRing() << " " << ring_hits.at(i).GetCharge() << "\t";
      std::cout << std::endl;
      std::cout << "Sector Energies " << std::endl;
      for(int i = 0; i < (int)sector_hits.size(); i++) if(sector_hits.at(i).GetCharge() > DThresh) std::cout << sector_hits.at(i).GetSector() << " " << sector_hits.at(i).GetCharge() << "\t";
      std::cout << std::endl;
    }
  }
}

/*******************************************************************************/
/* Gets hit position based on which ring and sector are hit ********************/
/* Use smear to get a uniform distribution over the size of the "Pixel" ********/
/*******************************************************************************/
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

/*******************************************************************************/
/* Required by TDetector, unused otherwise *************************************/
/*******************************************************************************/
void TJanus::InsertHit(const TDetectorHit& hit) {
  janus_hits.emplace_back((TJanusHit&)hit);
  fSize++;
}

/*******************************************************************************/
/* Basic print function - Not currently too useful *****************************/
/*******************************************************************************/
void TJanus::Print(Option_t *opt) const {
  std::cout << "TJanus @ " << GetTimestampNs() << std::endl;
  std::cout << "Size: " << GetSize() << std::endl;
  for(unsigned int i=0;i<GetSize();i++) {
    janus_hits.at(i).Print();
  }
  std::cout << "---------------------------" << std::endl;

}
