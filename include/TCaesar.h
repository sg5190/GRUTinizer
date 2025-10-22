#ifndef _TCAESAR_H_
#define _TCAESAR_H_

#include "TDetector.h"
#include "TCaesarHit.h"

class TCaesar : public TDetector {
public:
  TCaesar();
  virtual ~TCaesar();

  void Copy(TObject& obj) const;
  virtual void Clear(Option_t* opt = "");

  virtual TCaesarHit& GetCaesarHit(int i);
  virtual TDetectorHit& GetHit(int i);
  const TCaesarHit& GetAddbackHit(int i) const { return addback_hits.at(i); }

  virtual size_t Size() const { return caesar_hits.size(); }
  virtual Int_t AddbackSize(int SortDepth = 4) { BuildAddback(SortDepth); return addback_hits.size(); }
  void ResetAddback() { addback_hits.clear();}

  virtual void InsertHit(const TDetectorHit&);
  TCaesarHit* FindHit(int, bool);

  static TVector3 GetPosition(int chNum);
  static int GetRing(int chNum);
  static int GetDetectorNumber(int chNum);
  static void LoadDetectorPositions();
  int Npack() const { return num_packets;}
  int Nsize() const { return nsize;}

  static bool IsNeighbour(int ID1, int ID2) {SetCaesarNeighbours(); return caesarNeighbour[ID1][ID2];}
  static bool IsNeighbour(const TCaesarHit &a, const TCaesarHit &b, bool timegate=true) {
    bool tmpB = false;
    tmpB = IsNeighbour(a.GetChannel(),b.GetChannel());
    if(!timegate) return tmpB;
    if(a.GetNTimes() > 0 && b.GetNTimes() && abs(a.GetTime0()-b.GetTime0()) < 100) {
      return tmpB;
    } else return false;
  }

  // Allows for looping over all hits with for(auto& hit : caesar) { }
  std::vector<TCaesarHit>::iterator begin() { return caesar_hits.begin(); }
  std::vector<TCaesarHit>::iterator end() { return caesar_hits.end(); }

private:
  virtual int BuildHits(std::vector<TRawEvent>& raw_data);
  void BuildAddback(int SortDepth=4) const;

  int num_packets;
  int nsize;

  std::vector<TCaesarHit> caesar_hits;
  mutable std::vector<TCaesarHit> addback_hits;

  struct coord {
    TVector3 pos;
    int ring;
    int det;
  };
  static std::map<int,coord> detector_positions;
  static void SetCaesarNeighbours();
  static bool fNEIGHBOURSet;
  static bool caesarNeighbour[192][192];

  ClassDef(TCaesar,2);
};

#endif /* _TCAESAR_H_ */
