#include "TGRTEvent.h"

#include "TString.h"
#include <arpa/inet.h>
#include <endian.h>

ClassImp(TGRTEvent)

TGRTEvent::TGRTEvent() { }

TGRTEvent::TGRTEvent(const TRawEvent &raw) {
  raw.Copy(*this);
}

TGRTEvent::~TGRTEvent() { }

Long_t TGRTEvent::GetTimestamp() const {

  return be64toh(*((Long_t*)(TRawEvent::GetBody() + 0)));
}

Int_t TGRTEvent::GetSeqNum() const {
  return ntohs(((GRTroutingHdr*)(&fEventHeader))->seqnum());
}

Long_t TGRTEvent::GetCheckSum() const {
  return be64toh(*((Long_t*)(TRawEvent::GetBody() + 8)));
}

const char* TGRTEvent::GetPayload() const {
  if(fTimestamp != -1) {
    return fBody.GetData();
  } else {
    return fBody.GetData();
  }
}

TSmartBuffer TGRTEvent::GetPayloadBuffer() const {
  if(fTimestamp != -1){
    return fBody;
  } else {
    return fBody.BufferSubset(sizeof(Long_t)+sizeof(Long_t));
  }
}

void TGRTEvent::Clear(Option_t *opt) {
  TRawEvent::Clear(opt);
}

//void TGRTEvent::Print(Option_t *opt) const {
//  TString options(opt);
//  std::cout << BLUE << "Type:     \t" << DYELLOW << GetEventType() << RESET_COLOR << std::endl;
//  std::cout << BLUE << "SubType:  \t" << DYELLOW << GetSubType() << RESET_COLOR << std::endl;
//  std::cout << BLUE << "Size:     \t" << DYELLOW << GetBodySize()  << RESET_COLOR << std::endl;
//  std::cout << BLUE << "Seq Num   :\t" << DYELLOW << GetSeqNum() << RESET_COLOR << std::endl;
//  std::cout << BLUE << std::hex << "Timestamp:\t" << DYELLOW << GetTimestamp() << BLUE << "  tens of ns" << RESET_COLOR << std::endl;
//  std::cout << BLUE << std::hex << "Check Sum:\t" << DYELLOW << GetCheckSum() << RESET_COLOR << std::endl;
//  if(options.Contains("all")) {
//    TString pass_opt("bodyonly");
//    pass_opt.Append(options);
//    TRawEvent::Print(pass_opt.Data());
//  }
//}

