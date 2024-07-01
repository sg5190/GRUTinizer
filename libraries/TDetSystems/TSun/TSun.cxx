#include "TSun.h"

#include <algorithm>
#include <iostream>
#include <fstream>

#include "DDASDataFormat.h"
#include "DDASBanks.h"
#include "TNSCLEvent.h"
#include "TRawEvent.h"
#include "TChannel.h"
#include "TGRUTOptions.h"

/*******************************************************************************/
/* TSun ************************************************************************/
/*******************************************************************************/
TSun::TSun(){ }

TSun::~TSun(){ }

/*******************************************************************************/
/* Copies hit ******************************************************************/
/*******************************************************************************/
void TSun::Copy(TObject& obj) const {
  TDetector::Copy(obj);
  TSun& sun= (TSun&)obj;
  sun.sun_hits = sun_hits;
}

/*******************************************************************************/
/* Clear hit *******************************************************************/
/*******************************************************************************/
void TSun::Clear(Option_t* opt){
  TDetector::Clear(opt);
  sun_hits.clear();
}

/*******************************************************************************/
/* Functions to call Sun hits *************************************************/
/*******************************************************************************/
TSunHit& TSun::GetSunHit(int i){
  return sun_hits.at(i);
}

TDetectorHit& TSun::GetHit(int i){
  return sun_hits.at(i);
}

/*******************************************************************************/
/* Unpacks DDAS data and builds TSun events ****************************/
/*******************************************************************************/
int TSun::BuildHits(std::vector<TRawEvent>& raw_data) {

  //Loop over raw data
  for(auto& event : raw_data){
    TSmartBuffer buf = event.GetPayloadBuffer();
    TDDASEvent<DDASHeader> ddasevt(buf);
    unsigned int address = ( (70<<24) + (ddasevt.GetCrateID()<<16) + (ddasevt.GetSlotID()<<8) + ddasevt.GetChannelID() );

    //Check for channel address from .cal file otherwise skip this event
    TChannel* chan = TChannel::GetChannel(address);
    static int lines_displayed = 0;
    if(!chan){
      if(lines_displayed < 10) {
        std::cout << "Unknown DDAS (crate, slot, channel): (" << ddasevt.GetCrateID() << ", " << ddasevt.GetSlotID() << ", " << ddasevt.GetChannelID() << ")0x" << std::hex << address << std::dec << std::endl;
      }
      lines_displayed++;
      continue;
    }

    //Create a TSunHit
    TSunHit hit;
    hit.SetAddress(address);
    hit.SetDetectorNumber(chan->GetArrayPosition());
    hit.SetCharge(ddasevt.GetEnergy());
    hit.SetTimestamp(ddasevt.GetTimestamp()); //Timestamp in ns
    hit.SetExtTime(ddasevt.GetExtTimestamp());
    hit.SetTime(ddasevt.GetTime()); // Timestamp + CFD
    hit.SetCFDTime(ddasevt.GetCFDTime()); // CFD ONLY
    hit.SetTrace(ddasevt.GetTraceLength(), ddasevt.trace);
    sun_hits.push_back(hit);
    fSize++;
  }
  return Size();
}

/*******************************************************************************/
/* Required by TDetectorFactory ************************************************/
/*******************************************************************************/
void TSun::InsertHit(const TDetectorHit& hit) {
  sun_hits.emplace_back((TSunHit&)hit);
  fSize++;
}
