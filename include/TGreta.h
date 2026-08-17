#ifndef TGRETA_H
#define TGRETA_H

#ifndef __CINT__
#include <functional>
#endif

#include <TObject.h>
#include <TMath.h>

#include "TDetector.h"
#include "TGretaHit.h"

class TGreta : public TDetector {

public:
  TGreta();
  ~TGreta();

  virtual void Copy(TObject& obj) const;
  virtual void Print(Option_t *opt = "") const;
  virtual void Clear(Option_t *opt = "");

  virtual Int_t GetAddbackSize(int SortDepth = 6, bool AllHits = false) { BuildAddback(SortDepth, AllHits); return addback_hits.size(); }
  void ResetAddback() {
    addback_hits.clear();
  }

  virtual void InsertHit(const TDetectorHit& hit);
  virtual TDetectorHit& GetHit(int i)            { return greta_hits.at(i); }

  const TGretaHit& GetGretaHit(int i) const { return greta_hits.at(i); }
        TGretaHit& GetGretaHit(int i)       { return greta_hits.at(i); }
  const TGretaHit& GetAddbackHit(int i) const { return addback_hits.at(i); }


  void PrintHit(int i){ greta_hits.at(i).Print(); }

  static TVector3 CrystalToGlobal(int cryId,Float_t localX=0,Float_t localY=0,Float_t localZ=0);
  static TVector3 GetSegmentPosition(int cryid,int segment); //return the position of the segemnt in the lab system
  static TVector3 GetCrystalPosition(int cryid); //return the position of the crysal in the lab system

  static bool sortcrymap(std::pair<int, TVector3> i, std::pair<int, TVector3> j) {
    if(std::round(i.second.Theta()*1000.)/1000. == std::round(j.second.Theta()*1000.)/1000.) {
      return (i.second.Phi() < 0 ? i.second.Phi() + 2*TMath::Pi():i.second.Phi())  < (j.second.Phi() < 0 ? j.second.Phi() + 2*TMath::Pi():j.second.Phi());
    }
    else return i.second.Theta() < j.second.Theta();
  }

  static int GetSpecId(int cryid, bool inverse); //Return CrystalID ordered by theta depending on channel map file

  static bool IsNeighbour(int ID1, int ID2) {SetGretNeighbours(); return gretNeighbour[ID1][ID2];}
  static bool IsNeighbour(const TGretaHit &a, const TGretaHit &b, bool timegate=true) {
    bool tmpB = false;
    if(a.GetCrystalId() == 92 || b.GetCrystalId() == 92) return false;
    tmpB = IsNeighbour(a.GetCrystalId(),b.GetCrystalId());
    if(!timegate) return tmpB;
    if(abs(a.GetTimestamp()-b.GetTimestamp()) < 100) {
      return tmpB;
    } else return false;
  }

//  static bool IsNeighbour(const TGretaHit &a, const TGretaHit &b, bool timegate=false) {
//    bool tmpB = false;
//    tmpB = IsNeighbour(a.GetCrystalId(),b.GetCrystalId());
//    if(abs(a.GetTimestamp()-b.GetTimestamp()) < 100) {
//      if(timegate) return tmpB;
//      if(a.GetHoleNumber() == b.GetHoleNumber()) return tmpB;
//      else return false;
//    } else return false;
//  }

  const std::vector<TGretaHit> &GetAllHits() const { return greta_hits; }

  void  SortHits();

  float MSegPos(int type,int seg,int coord) const {
    return m_segpos[type][seg][coord];
  }

private:
  void BuildAddback(int SortDepth=6, bool AllHits=false) const;
  virtual int BuildHits(std::vector<TRawEvent>& raw_data);

  std::vector<TGretaHit> greta_hits;
  mutable std::vector<TGretaHit> addback_hits; //!
  static Float_t crmat[32][4][4][4];
  static Float_t m_segpos[2][36][3];
  static void SetGretNeighbours();
  static bool fNEIGHBOURSet;
  static bool gretNeighbour[124][124];
  static void SetCRMAT();
  static void SetSegmentCRMAT();
  static bool fCRMATSet;
  static void SetSpecId();
  static std::map<int, int> fThetaCryMap;
  static bool fSetSpec;
  static bool fFileDetermined;
  static bool fGRT;

  ClassDef(TGreta,3);
};

#endif
