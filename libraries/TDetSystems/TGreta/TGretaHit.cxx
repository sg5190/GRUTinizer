#include "TGretaHit.h"

#include <algorithm>
#include <cmath>
#include <set>

#include <TRandom.h>

#include "GValue.h"
#include "TGEBEvent.h"
#include "TGreta.h"
//#include "TS800.h"

TGretaHit::TGretaHit(){ Clear(); }

TGretaHit::~TGretaHit(){ }

bool TGretaHit::fSmear = false;
double TGretaHit::fSmearWidth = 0.2; // 2 mm
/*******************************************************************************/
/* Copies hit ******************************************************************/
/*******************************************************************************/
void TGretaHit::Copy(TObject &rhs) const {
  TDetectorHit::Copy(rhs);
  ((TGretaHit&)rhs).fT0             = fT0;
  ((TGretaHit&)rhs).fPad            = fPad;
  ((TGretaHit&)rhs).fCrystalId      = fCrystalId;
  ((TGretaHit&)rhs).fCoreCharge[0]  = fCoreCharge[0];
  ((TGretaHit&)rhs).fCoreCharge[1]  = fCoreCharge[1];
  ((TGretaHit&)rhs).fCoreCharge[2]  = fCoreCharge[2];
  ((TGretaHit&)rhs).fCoreCharge[3]  = fCoreCharge[3];
  ((TGretaHit&)rhs).fSegments       = fSegments;
  ((TGretaHit&)rhs).fAB	            = fAB;
  ((TGretaHit&)rhs).fChisq          = fChisq;
  ((TGretaHit&)rhs).fNormChisq      = fNormChisq;
  ((TGretaHit&)rhs).fBaseline	    = fBaseline;
  ((TGretaHit&)rhs).fPrestep	    = fPrestep;
  ((TGretaHit&)rhs).fPoststep	    = fPoststep;
}

/*******************************************************************************/
/* Gets mnemonic for channel if a channels.cal file is provided ****************/
/*******************************************************************************/
const char *TGretaHit::GetName() const {
  TChannel *channel = TChannel::GetChannel(GetAddress());
  if(channel) return channel->GetName();
  return "";
}


/*******************************************************************************/
/* Returns position vector for the centre of the crystal ***********************/
/*******************************************************************************/
TVector3 TGretaHit::GetCrystalPosition()  const {
  return TGreta::GetCrystalPosition(fCrystalId);
}

/*******************************************************************************/
/* Build TGreta events from mode2 data ***************************************/
/*******************************************************************************/
void TGretaHit::BuildFrom(TSmartBuffer& buf){
  const TRawEvent::GEBBankType1& raw = *(const TRawEvent::GEBBankType1*)buf.GetData();

  //std::cout << "GretaHit: " << raw << std::endl;
  //std::cout << raw << std::endl;
  Clear();

  SetTimestamp(raw.timestamp, 10.);
  SetEnergy(raw.tot_e, false);
  fT0 = raw.t0;
  fCrystalId = raw.crystal_id;
  fAB = -1;

  int board_id = ((fCrystalId/4) << 8) ;  //hole  number : 0x1f00
  board_id += ((fCrystalId%4) << 6) ;  //x-tal number : 0x00c0
  board_id += 9;                       //chan  number : 0x000f  information not available here(assume core).
  SetAddress((1<<24) + board_id);

  TChannel *channel = TChannel::GetChannel(GetAddress());
  if(!channel) {
    TChannel *newchan = new TChannel(Form("GRR%2.2iXN%2.2iX", fCrystalId/4, fCrystalId%4),GetAddress());
    TChannel::AddChannel(newchan);
  }

  for(int i=0; i<4; i++){
    fCoreCharge[i] = raw.core_e[i];
  }
  fChisq = raw.chisq;
  fNormChisq = raw.norm_chisq;
  fBaseline = raw.baseline;
  fPrestep = raw.prestep;
  fPoststep = raw.poststep;

  fPad = raw.pad;
  for(int i=0; i<raw.num; i++) {
    try {
      interaction_point pnt;
      pnt.fSeg = raw.intpts[i].seg;
      pnt.fX   = raw.intpts[i].x;
      pnt.fY   = raw.intpts[i].y;
      pnt.fZ   = raw.intpts[i].z;
      pnt.fFrac   = raw.intpts[i].e;
      pnt.fEng = raw.intpts[i].seg_ener;
      fSegments.push_back(pnt);
    } catch(...) {
      std::cout << "in try catch block!" << std::endl;
    }
  }
  std::sort(fSegments.begin(),fSegments.end());
  //  Print("all");
}

/*******************************************************************************/
/* Returns position vector based on interaction point **************************/
/*******************************************************************************/
TVector3 TGretaHit::GetIntPosition(unsigned int i) const {
  if(i<fSegments.size()){
    double xoffset = GValue::Value("GRETA_X_OFFSET");
    if(std::isnan(xoffset))
      xoffset=0.00;
    double yoffset = GValue::Value("GRETA_Y_OFFSET");
    if(std::isnan(yoffset))
      yoffset=0.00;
    double zoffset = GValue::Value("GRETA_Z_OFFSET");
    if(std::isnan(zoffset))
      zoffset=0.00;
    TVector3 v = TGreta::CrystalToGlobal(fCrystalId,fSegments.at(i).fX + xoffset,
        fSegments.at(i).fY + yoffset,
        fSegments.at(i).fZ + zoffset);
    return v;
  } else {
    return TDetectorHit::BeamUnitVec;
  }
}
/*******************************************************************************/
/* Returns vector for loacl position within the crystal ************************/
/*******************************************************************************/
TVector3 TGretaHit::GetLocalPosition(unsigned int i) const {
  if(i<fSegments.size()){
    return TVector3(fSegments.at(i).fX,
        fSegments.at(i).fY,
        fSegments.at(i).fZ);
  } else {
    return TVector3(0,0,1);
  }
}

/*******************************************************************************/
/* Returns doppler corrected energies based on highest energy interaction ******/
/* point, particle vectors and beam beta. **************************************/
/*******************************************************************************/
double TGretaHit::GetDoppler(double beta,const TVector3 *vec) const {
  if(GetNSegments()<1)
    return 0.0;
  if(vec==0) {
    vec = &BeamUnitVec;
  }
  double tmp = 0.0;
  double gamma = 1/(sqrt(1-pow(beta,2)));

  TVector3 gret_pos = GetPosition();
  double xoffset = GValue::Value("TARGET_X_OFFSET");
  if(std::isnan(xoffset)) xoffset=0.00;
  double yoffset = GValue::Value("TARGET_Y_OFFSET");
  if(std::isnan(yoffset)) yoffset=0.00;
  double zoffset = GValue::Value("TARGET_Z_OFFSET");
  if(std::isnan(zoffset)) zoffset=0.00;

  gret_pos.SetX(gret_pos.X() - xoffset);
  gret_pos.SetY(gret_pos.Y() - yoffset);
  gret_pos.SetZ(gret_pos.Z() - zoffset);

  tmp = GetEnergy()*gamma *(1 - beta*TMath::Cos(gret_pos.Angle(*vec)));
  return tmp;
}

/*******************************************************************************/
/* Returns doppler corrected energies based on highest energy interaction ******/
/* point, particle vectors and beam beta. Corrects for non-dispersive position */
/* YTA, calculated using TS800 class. ******************************************/
/*******************************************************************************/
double TGretaHit::GetDopplerYta(double beta, double yta, const TVector3 *vec) const {
  if(GetNSegments()<1)
    return 0.0;
  if(vec==0) {
    vec = &BeamUnitVec;
  }
  //Target offsets determine new reference point in lab frame
  double xoffset = GValue::Value("TARGET_X_OFFSET");
  if(std::isnan(xoffset)) xoffset=0.00;
  double yoffset = GValue::Value("TARGET_Y_OFFSET");
  if(std::isnan(yoffset)) yoffset=0.00;
  double zoffset = GValue::Value("TARGET_Z_OFFSET");
  if(std::isnan(zoffset)) zoffset=0.00;
  return GetDopplerYta(beta, yta, xoffset, yoffset, zoffset, vec);
}

double TGretaHit::GetDopplerYta(double beta, double yta, double xoffset, double yoffset, double zoffset, const TVector3 *vec) const {
  if(GetNSegments()<1)
    return 0.0;
  if(vec==0) {
    vec = &BeamUnitVec;
  }
  double tmp = 0.0;
  double gamma = 1./(sqrt(1.-pow(beta,2.)));
  TVector3 gret_pos = GetPosition();

  //Target offsets determine new reference point in lab frame
  gret_pos.SetX(gret_pos.X() - xoffset);
  gret_pos.SetY(gret_pos.Y() - (yoffset - yta));
  gret_pos.SetZ(gret_pos.Z() - zoffset);

  tmp = GetEnergy()*gamma *(1 - beta*TMath::Cos(gret_pos.Angle(*vec)));
  return tmp;
}

/*******************************************************************************/
/* Used by Nearest Neighbour Addback code to add energies **********************/
/* and handle interaction points ***********************************************/
/*******************************************************************************/
void TGretaHit::NNAdd(const TGretaHit& rhs) {
  //copy original hit into singles
  if (!fSetFirstSingles){
    TGretaHit myCopy(*this);
    fSingles.push_back(myCopy);
    fSetFirstSingles = true;
  }
  SetEnergy((GetEnergy() + rhs.GetEnergy()), true);
  fSingles.push_back(rhs);

  // S.G. - Does it make sense to add interactions points from other crystal?
  // Fill all interaction points
  for(unsigned int i=0; i<rhs.fSegments.size(); i++){
    if(GetNSegments() >= MAXHPGESEGMENTS){
      break;
    }
    fSegments.push_back(rhs.fSegments[i]);
  }
}

/*******************************************************************************/
/* Returns position vector based on final interaction point ********************/
/*******************************************************************************/
TVector3 TGretaHit::GetLastPosition() const {
  if(fSegments.size()<1)
    return TDetectorHit::BeamUnitVec;
  return GetIntPosition(fSegments.size()-1);
}

/*******************************************************************************/
/* Print function **************************************************************/
/*******************************************************************************/
void TGretaHit::Print(Option_t *opt) const {

  std::cout << "TGretaHit:" <<  std::endl;
  std::cout << "\tAddress:        \t0x" << std::hex << fAddress << std::dec << std::endl;
  std::cout << "\tHole:           \t" << GetHoleNumber()               << std::endl;
  std::cout << "\tCrystalNum      \t" << GetCrystalNumber()            << std::endl;
  std::cout << "\tCrystalId:      \t" << GetCrystalId()                << std::endl;

  std::cout << "\tLocal Timestamp:\t" << GetTimestamp()                << std::endl;
  std::cout << "\tCorrected time: \t" << GetTime()                     << std::endl;
  std::cout << "\tDecomp Energy:  \t" << GetEnergy()                   << std::endl;
  if(!strcmp(opt,"all")) {
    for(int i=0;i<4;i++)
      std::cout << "\tCharge[" << i << "]:       \t" << fCoreCharge[i] << std::endl;
  }
  std::cout << "\tErrorCode:       \t" << fPad                         << std::endl;
  std::cout << "\tInteractions:    \t" << fSegments.size()             << std::endl;

  if(!strcmp(opt,"all")) {
    for(int i=0;i<GetNSegments();i++) {
      printf("\t\t");
      fSegments.at(i).Print();
    }
  }
  std::cout << "------------------------------"  << std::endl;

}

/*******************************************************************************/
/* Clears hit ******************************************************************/
/*******************************************************************************/
void TGretaHit::Clear(Option_t *opt) {
  TDetectorHit::Clear(opt);
  fT0	          = sqrt(-1);
  fCrystalId      = -1;
  fCoreCharge[0]  = -1;
  fCoreCharge[1]  = -1;
  fCoreCharge[2]  = -1;
  fCoreCharge[3]  = -1;
  fPad            = 0;
  fAB             = -1;
  fSegments.clear();
}


/*******************************************************************************/
/* Reduced fSegment size *******************************************************/
/* 0 - drop multiple ident int pnts ********************************************/
/* 1 - make into wedge data ***************************************************/
/*******************************************************************************/
void TGretaHit::TrimSegments(int type) {
  if(type==0) {
    std::set<interaction_point,intpnt_compare> pset;
    for(auto x=fSegments.begin();x!=fSegments.end();x++) {
      pset.insert(*x);
    }
    fSegments.clear();
    for(auto x=pset.begin();x!=pset.end();x++) {
      fSegments.push_back(*x);
    }
    std::sort(fSegments.begin(),fSegments.end());
  } else if (type==1) {
    std::set<interaction_point,intpnt_compare_wedge> pset;
    for(auto x=fSegments.begin();x!=fSegments.end();x++) {
      pset.insert(*x);
    }
    fSegments.clear();
    for(auto x=pset.begin();x!=pset.end();x++) {
      fSegments.push_back(*x);
      fSegments.back().fSeg = (fSegments.back().fSeg%6);
    }
    std::sort(fSegments.begin(),fSegments.end());
  }
}

/*******************************************************************************/
/* Returns position vector to first interaction point **************************/
/* For simulations the position can be smeared with the functions **************/
/* SetSmear(true/false) and SetSmearWidth(width) width should be in cm *********/
/*******************************************************************************/
TVector3 TGretaHit::GetPosition() const {
  TVector3 gret_pos = GetIntPosition(0);
  if(!fSmear) {
  } else {
    gret_pos.SetX(gret_pos.X() - gRandom->Gaus(0,fSmearWidth));
    gret_pos.SetY(gret_pos.Y() - gRandom->Gaus(0,fSmearWidth));
    gret_pos.SetZ(gret_pos.Z() - gRandom->Gaus(0,fSmearWidth));
  }
  return gret_pos;
}

/*******************************************************************************/
/* Remaps crystal ids in order of Theta/Phi position depending on which ********/
/* crystals are in the data stream. Value is equivalent to SpecTcl Order *******/
/* Useful for ordering crystals by theta to calculare target offsets using *****/
/* find_optimal_shifts.C macro or omparing data to SpecTcl *********************/
/*******************************************************************************/
Int_t TGretaHit::GetSpecOrder() const {
  return TGreta::GetSpecId(fCrystalId, false);
}

/*******************************************************************************/
/* Inverse function of GetSpecOrder() returns GetCrystalId() for a given *******/
/* SpecTcl Order ***************************************************************/
/*******************************************************************************/
Int_t TGretaHit::GetInverseSpecOrder() const {
  return TGreta::GetSpecId(fCrystalId, true);
}

/*******************************************************************************/
/* TInteractionPoint ***********************************************************/
/* Used by TCluster to do tracking *********************************************/
/* Currently not used for FRIB fast beam analysis ******************************/
/*******************************************************************************/
/*
TInteractionPoint::TInteractionPoint(const TInteractionPoint &IP) {
  this->Copy(IP);
}

void TInteractionPoint::Copy(const TInteractionPoint &IP) {
  fSegNum      = IP.GetSegNum();
  fEng         = IP.GetPreampE();
  fDecompEng   = IP.GetDecompE();
  fAssignedEng = IP.GetAssignE();
  fOrder       = IP.GetOrder();
  fLPosition   = IP.GetLocalPosition();
}

TVector3 TInteractionPoint::GetPosition(int xtal) const {
  return TGreta::CrystalToGlobal(xtal, fLPosition.X(), fLPosition.Y(), fLPosition.Z());
}


void TInteractionPoint::Print(Option_t *opt) const {
  TVector3 lv = GetLocalPosition();
  printf("seg[%02i] wedge[%i]    %.1f / %.1f   [ %.1f, %.1f. %.1f] \n",
          GetSegNum(),Wedge(),GetAssignE(),GetPreampE(),lv.X(),lv.Y(),lv.Z());
}


void TInteractionPoint::Clear(Option_t *opt) {
    fSegNum      = -1;
    fEng         = sqrt(-1);
    fDecompEng   = sqrt(-1);
    fAssignedEng = sqrt(-1);
    fOrder       = -1;
    fLPosition.SetXYZ(0,0,1);
}
*/
