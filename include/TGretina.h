#ifndef TGRETINA_H
#define TGRETINA_H

#ifndef __CINT__
#include <functional>
#endif

#include <TObject.h>
#include <TMath.h>


#include "TDetector.h"
#include "TGretinaHit.h"

class TGretina : public TDetector {

public:
  TGretina();
  ~TGretina();

  virtual void Copy(TObject& obj) const;
  virtual void Print(Option_t *opt = "") const;
  virtual void Clear(Option_t *opt = "");

  virtual size_t Size() const { return gretina_hits.size(); }
  virtual Int_t AddbackSize(int SortDepth = 6, int EngRange = -1) { BuildAddback(SortDepth, EngRange); return addback_hits.size(); }
  void ResetAddback() { addback_hits.clear();}

  virtual void InsertHit(const TDetectorHit& hit);
  virtual TDetectorHit& GetHit(int i)            { return gretina_hits.at(i); }

  const TGretinaHit& GetGretinaHit(int i) const { return gretina_hits.at(i); }
        TGretinaHit& GetGretinaHit(int i)       { return gretina_hits.at(i); }
  const TGretinaHit& GetAddbackHit(int i) const { return addback_hits.at(i); }


  void PrintHit(int i){ gretina_hits.at(i).Print(); }

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
  static bool IsNeighbour(const TGretinaHit &a, const TGretinaHit &b, bool timegate=true) {
    bool tmpB = false;
    tmpB = IsNeighbour(a.GetCrystalId(),b.GetCrystalId());
    if(!timegate) return tmpB;
    if(abs(a.GetTime()-b.GetTime()) < 50) {
      return tmpB;
    } else return false;
  }

  const std::vector<TGretinaHit> &GetAllHits() const { return gretina_hits; }

  void  SortHits();

  void  CleanHits(int i=-1) {
    for(auto x=gretina_hits.begin();x!=gretina_hits.end();) { //x++) {
      if(x->GetPad()==0) {
        x->TrimSegments(1);
        if(i>=0) x->SetCoreEnergy(x->GetCoreEnergy(i));
        x++;
      } else {
        x = gretina_hits.erase(x);
      }
    }
  }

  double SumHits(bool clean=true) {
    double sum =0.0;
    for(auto x=gretina_hits.begin();x!=gretina_hits.end();x++) {
      if(clean && x->GetPad()==0) {
        sum += x->GetCoreEnergy();
      } else if(!clean) {
        sum += x->GetCoreEnergy();
      }
    }
    return sum;
  }

  float MSegPos(int type,int seg,int coord) const {
    return m_segpos[type][seg][coord];
  }

private:
  void BuildAddback(int SortDepth=6, int EngRange = -1) const;
  virtual int BuildHits(std::vector<TRawEvent>& raw_data);

  std::vector<TGretinaHit> gretina_hits;
  mutable std::vector<TGretinaHit> addback_hits; //!

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
  ClassDef(TGretina,3);
};





#endif
