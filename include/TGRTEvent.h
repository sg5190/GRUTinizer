#ifndef _TGRTEVENT_H_
#define _TGRTEVENT_H_

#include "TRawEvent.h"


class TGRTEvent : public TRawEvent {
public:
  TGRTEvent();
  TGRTEvent(const TRawEvent&);
  virtual ~TGRTEvent();


  Int_t GetSeqNum() const;

  long GetTimestamp() const;
  long GetCheckSum() const;

  const char* GetPayload() const;
  TSmartBuffer GetPayloadBuffer() const;

  virtual void Clear(Option_t *opt ="");
//  virtual void Print(Option_t *opt ="") const;

  ClassDef(TGRTEvent, 0);
};

#endif /* _TGRTEVENT_H_ */
