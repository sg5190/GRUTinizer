#ifndef _TSUN_H_
#define _TSUN_H_

#include "TDetector.h"
#include "TSunHit.h"

class TSun : public TDetector {
public:
  TSun();
  virtual ~TSun();

  void Copy(TObject& obj) const;

  virtual void Clear(Option_t* opt = "");
  virtual TSunHit& GetSunHit(int i);
  virtual TDetectorHit& GetHit(int i);
  virtual void InsertHit(const TDetectorHit&);

  // Allows for looping over all hits with for(auto& hit : sun) { }
  std::vector<TSunHit>::iterator begin() { return sun_hits.begin(); }
  std::vector<TSunHit>::iterator end() { return sun_hits.end(); }

private:
  virtual int BuildHits(std::vector<TRawEvent>& raw_data);

  std::vector<TSunHit> sun_hits;
  ClassDef(TSun,2);
};

#endif /* _TSUN_H_ */
