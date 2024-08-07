
#include <TLenda.h>

#include <TRawEvent.h>
#include <DDASDataFormat.h>
#include <DDASBanks.h>

ClassImp(TLenda)

int TLenda::TDiff = 1000;

void TLenda::Copy(TObject& obj) const {
  TDetector::Copy(obj);

  TLenda& lenda = (TLenda&)obj;
  lenda.bar_hits = bar_hits;
  lenda.top_hits = top_hits;
  lenda.bottom_hits = bottom_hits;
  lenda.ref_hits = ref_hits;
}

void TLenda::Clear(Option_t* opt){
  TDetector::Clear(opt);
  bar_hits.clear();
  top_hits.clear();
  bottom_hits.clear();
  ref_hits.clear();
}

int TLenda::BuildHits(std::vector<TRawEvent>& raw_data) {
  for(auto& event : raw_data){
    SetTimestamp(event.GetTimestamp());

    TSmartBuffer buf = event.GetPayloadBuffer();

    int total_size = *(int*)buf.GetData();
    const char* buffer_end = buf.GetData() + total_size;
    buf.Advance(sizeof(int));

    while(buf.GetData() < buffer_end) {
      TLendaHit hit;

      // Constructor advances the buffer to end of each channel
      TDDASEvent<DDASGEBHeader> ddas(buf);
      unsigned int address = ( (21 << 24) + (ddas.GetCrateID() << 16) + (ddas.GetSlotID() << 8) + ddas.GetChannelID() );
      TChannel* chan = TChannel::GetChannel(address);
      static int lines_displayed = 0;
      if(!chan){
        if(lines_displayed < 20) {
          std::cout << "Unknown LendaChannel Address 0x" << std::hex << address << std::dec << std::endl;
        }
        lines_displayed++;
        continue;
      }
      hit.SetAddress(address);
      hit.SetCharge(ddas.GetEnergy());
      hit.SetTimestamp(ddas.GetTimestamp());
      hit.SetExtTimestamp(ddas.GetExtTimestamp());
      hit.SetTime(ddas.GetTime()); // Timestamp + CFD
      hit.SetCFDTime(ddas.GetCFDTime()); // CFD ONLY
      hit.SetCFDFail(ddas.GetCFDFailBit());
//std::cout << ddas.GetCFDTime() << "\t" << hit.CFDTime() << std::endl;
//std::cout << ddas.GetCFDTrig() << "\t" << hit.GetCFDFail() << std::endl;
      hit.SetTrace(ddas.GetTraceLength(), ddas.trace);

      hit.SetPosArray(chan->GetArrayPosition());
      hit.SetPosSeg(chan->GetSegment());

      if( strcmp(chan->GetSystem(), "LET") ==0) {
	top_hits.push_back(std::move(hit));
      } else if( strcmp(chan->GetSystem(), "LEB") ==0) {
	bottom_hits.push_back(std::move(hit));
      } else if( strcmp(chan->GetSystem(), "LER") ==0) {
	ref_hits.push_back(std::move(hit));
       }
    }
  }
  return top_hits.size() + bottom_hits.size() + ref_hits.size();
}

//Vector of Top Hits
TLendaHit* TLenda::GetTopHit(const int&i) {
  if(i < GetTopSize()) {
    return &top_hits.at(i);
  }
  return nullptr;
}

//Vector of Top Hits
TLendaHit* TLenda::GetBottomHit(const int&i) {
  if(i < GetBottomSize()) {
    return &bottom_hits.at(i);
  }
  return nullptr;
}

TLendaHit* TLenda::GetReferenceHit(const int&i) {
  if(i < GetReferenceSize()) {
    return &ref_hits.at(i);
  }
  return nullptr;
}


//Size of TLendaBar
int TLenda::GetBarSize() {
  BuildLendaBars();
  return bar_hits.size();
}

void TLenda::BuildLendaBars() {
  bar_hits.clear();
  if(top_hits.size() == 0 || bottom_hits.size() == 0) return;

  for(size_t i = 0; i < top_hits.size(); i++) {
    for(size_t j = 0; j < bottom_hits.size(); j++) {
      //Check the top and bottom hit are the same bar
      if(top_hits.at(i).GetDetnum() == bottom_hits.at(j).GetDetnum()) {
	//Checks for good time
        if(top_hits.at(i).Timestamp() - bottom_hits.at(j).Timestamp() < TDiff) {
	  //Use top hit for all data
	  TLendaHit dhit = top_hits.at(i);
	  bar_hits.push_back(dhit);
	}
      }
    }
  }
}

//Vector of build bar hits, must call GetBarSize/Build Bars first
TLendaHit* TLenda::GetBarHit(const int&i) {
  if(i < GetBarSize()) {
    return &bar_hits.at(i);
  }
  return nullptr;
}
