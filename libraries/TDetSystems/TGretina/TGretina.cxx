#include <deque>
#include <fstream>
#include <string>
#include <sstream>
#include <set>

#include "TGretina.h"
#include "TGEBEvent.h"

/*******************************************************************************/
/* TGretina ********************************************************************/
/* Main Class for Gretina mode2 (Decomposed) analysis **************************/
/*******************************************************************************/
TGretina::TGretina(){
  Clear();
}

TGretina::~TGretina() {
}

Float_t TGretina::crmat[32][4][4][4];
Float_t TGretina::m_segpos[2][36][3];
bool    TGretina::fCRMATSet = false;
bool    TGretina::fNEIGHBOURSet = false;
bool	TGretina::gretNeighbour[124][124];
std::map<int, int> TGretina::fThetaCryMap;
bool    TGretina::fSetSpec = false;

/*******************************************************************************/
/* Main function for making addback hits according to procedure defined in *****/
/* D. Weisshaar et al., Nucl. Instrum. Methods Phys. Res., Sect. A 847, ********/
/* 18, (2017). Sec 3.3 *********************************************************/
/* SortDepth determines how many combinations of events to consider default=4 **/
/* AddbackDepth classifies events as, n0 (0), n1 (1), n2 (2) or ng (3) *********/
/*******************************************************************************/
//SG TO DO
//How to handle interaction points from added crystals, does it make sense to add them?
void TGretina::BuildAddback(int SortDepth, int EngRange) const {
  //Don't rebuild events
  if( addback_hits.size() > 0 || gretina_hits.size() == 0) {
    return;
  }

  std::vector<TGretinaHit> temp_hits = gretina_hits;
  for(unsigned int i = 0; i < temp_hits.size(); i++) {
    if(EngRange < 0) temp_hits.at(i).SetCoreEnergy(temp_hits.at(i).GetCoreEnergy());
    else temp_hits.at(i).SetCoreEnergy(temp_hits.at(i).GetCoreEnergy(EngRange));
  }

  //sort so that the first hit has the greatest energy
  //this way we can loop through i,j with i < j and know that
  //any hit with higher energy cannot be an addback to one with lower energy
  std::sort(temp_hits.begin(), temp_hits.end(),
    [](const TGretinaHit& a, const TGretinaHit& b) {
    return a.GetCoreEnergy() > b.GetCoreEnergy();
  });

  //vector used to store hit indices when crystals are pairs
  std::vector<int> paired;
  //loop through every hit
  for(unsigned int i=0; i < temp_hits.size(); i++) {
    paired.clear();
    TGretinaHit &current_hit = temp_hits[i];
    //only pick unqiue pairs of hits
    if(SortDepth > 1) {
      for(unsigned int j=i+1; j < temp_hits.size(); j++) {
        if (IsNeighbour(current_hit,temp_hits[j])){
	  paired.push_back(j);
          //check every hit k to see if it is a neighbour with j
	  if(SortDepth > 2) {
            for (unsigned int k=0; k < temp_hits.size(); k++){
	      if( (i == k) || (j == k) ) continue;
              if (IsNeighbour(temp_hits[j],temp_hits[k])){
	        paired.push_back(k);
	        //Edge cases, requires at least 4 crystals in a line
	        if(SortDepth > 3) {
	          for (unsigned int l=0; l < temp_hits.size(); l++){
                    if( (i == l) || (j == l) || (k == l)) continue;
	            if (IsNeighbour(temp_hits[k],temp_hits[l])) {
		      paired.push_back(l);
 	              //Really Edge cases, requires at least 5 crystals in a line
		      if(SortDepth > 4) {
	                for (unsigned int m=0; m < temp_hits.size(); m++){
		          if( (i == m) || (j == m) || (k == m) || (l == m)) continue;
	                  if(IsNeighbour(temp_hits[l],temp_hits[m])) {
		            paired.push_back(m);
 	                    //Extreme Edge cases, requires at least 6 crystals in a line
			    if(SortDepth > 5) {
                              for (unsigned int n=0; n < temp_hits.size(); n++){
                                if( (i == n) || (j == n) || (k == n) || (l == n) || (m == n)) continue;
                                if(IsNeighbour(temp_hits[m],temp_hits[n])) paired.push_back(n);
			      }
		            }
		          }
			}
		      }
	            }
		  }
	        }
              }
	    }
          }
        }
      }
    }

    //Reverse sort paired vector required when removing hits later
    std::sort(paired.rbegin(), paired.rend());
    //Vector can contain multiple version of same hit, removes duplicate values
    paired.erase(unique(paired.begin(), paired.end()), paired.end());
    if (paired.size() == 0) { //n0 events
      addback_hits.push_back(current_hit);
      addback_hits.back().SetABDepth(0);
    } else if (paired.size() == 1) { //n1 Events
      current_hit.NNAdd(temp_hits[paired.at(0)]);
      addback_hits.push_back(current_hit);
      addback_hits.back().SetABDepth(1);
      //Erase hit after adding otherwise event can make a new addback event
      temp_hits.erase(temp_hits.begin() + paired.at(0));
    }
    else if (paired.size() == 2) { //n2 Or Ng
      if(IsNeighbour(current_hit,temp_hits[paired.at(0)]) && IsNeighbour(current_hit,temp_hits[paired.at(1)]) && IsNeighbour(temp_hits[paired.at(0)],temp_hits[paired.at(1)]) ) { //n2
        current_hit.NNAdd(temp_hits[paired.at(1)]);
        current_hit.NNAdd(temp_hits[paired.at(0)]);
        addback_hits.push_back(current_hit);
        addback_hits.back().SetABDepth(2);
        //Erase hit after adding otherwise event can make a new addback event
        temp_hits.erase(temp_hits.begin() + paired.at(0));
        temp_hits.erase(temp_hits.begin() + paired.at(1));
      } else { //Also ng
          addback_hits.push_back(current_hit);
          addback_hits.back().SetABDepth(3);
          for(int p = 0; p < (int)paired.size(); p++) {
          addback_hits.push_back(temp_hits[paired.at(p)]);
          addback_hits.back().SetABDepth(3);
          //Erase hit after adding otherwise event can make a new addback event
          temp_hits.erase(temp_hits.begin() + paired.at(p));
        }
      }
    } else { //ng
      addback_hits.push_back(current_hit);
      addback_hits.back().SetABDepth(3);
      for(int p = 0; p < (int)paired.size(); p++) {
        addback_hits.push_back(temp_hits[paired.at(p)]);
        addback_hits.back().SetABDepth(3);
        //Erase hit after adding otherwise event can make a new addback event
        temp_hits.erase(temp_hits.begin() + paired.at(p));
      }
    }
  }

  return;
}

/*******************************************************************************/
/* Loads files containing rotation matrices of core required for position ******/
/*******************************************************************************/
void TGretina::SetCRMAT() {
  if(fCRMATSet){
    return;
  }

  FILE *fp;
  std::string temp = getenv("GRUTSYS");
  temp.append("/libraries/TDetSystems/TGretina/crmat.dat");
  const char *fn = temp.c_str();
  float f1, f2, f3, f4;
  int pos, xtal;
  int nn = 0;
  char *st, str[256];
  fp = fopen64(fn, "r");
  if (fp == NULL) {
    printf("Could not open \"%s\".\n", fn);
    exit(1);
  }
  //printf("\"%s\" open....", fn);
  /* Read values. */
  nn = 0;
  st = fgets(str, 256, fp);
  while (st != NULL) {
    if (str[0] == 35) {
      /* '#' comment line, do nothing */
    } else if (str[0] == 59) {
      /* ';' comment line, do nothing */
    } else if (str[0] == 10) {
      /* Empty line, do nothing */
    } else {
      sscanf(str, "%i %i", &pos, &xtal);
      for(int i=0; i<4; i++) {
        st = fgets(str, 256, fp);
        sscanf(str, "%f %f %f %f", &f1, &f2, &f3, &f4);
        crmat[pos-1][xtal][i][0] = f1;
        crmat[pos-1][xtal][i][1] = f2;
        crmat[pos-1][xtal][i][2] = f3;
        crmat[pos-1][xtal][i][3] = f4;
      }
      nn++;
    }
    /* Attempt to read the next line. */
    st = fgets(str, 256, fp);
  }
  SetSegmentCRMAT();
  fCRMATSet = true;
  //printf("Read %i rotation matrix coefficients.\n", nn);
  /* Done! */
}

/*******************************************************************************/
/* Loads Segment Positions - Typically unused as position is derived from ******/
/* interaction points **********************************************************/
/*******************************************************************************/
void TGretina::SetSegmentCRMAT() {
  if(fCRMATSet)
    return;
  //FILE *infile;
  int NUMSEG = 36;
  std::string temp = getenv("GRUTSYS");
  temp.append("/libraries/TDetSystems/TGretina/crmat_segpos.dat");
  std::ifstream infile;
  infile.open(temp.c_str());
  if(!infile.is_open()) {
    return;
  }
  // notice: In file type-A is first, then type-B but as type-B are even
  //         det_ids in the data stream we define type=0 for type-B not
  //         type-A.
  std::string line;
  int type = 0;
  int seg  = 0;
  while(getline(infile,line)) {
    if(seg==NUMSEG) {
      seg  = 0;
      type = 1;
    }
    std::stringstream ss(line);
    int segread;
    double x,y,z;
    ss >> segread >> x >> y >> z;
    if((seg+1)!=segread) {
      fprintf(stderr,"%s: seg[%i] read but seg[%i] expected.\n",__PRETTY_FUNCTION__,segread,seg+1);
    }
    m_segpos[(type+1)%2][seg][0] = x;
    m_segpos[(type+1)%2][seg][1] = y;
    m_segpos[(type+1)%2][seg][2] = z;
    seg++;
  }
  infile.close();
}

/*******************************************************************************/
/* Gets Actual Position of Segment *********************************************/
/*******************************************************************************/
TVector3 TGretina::GetSegmentPosition(int cry_id,int segment) {
  SetCRMAT();
  float x = m_segpos[cry_id%2][segment][0];
  float y = m_segpos[cry_id%2][segment][1];
  float z = m_segpos[cry_id%2][segment][2];
  TVector3 v(x,y,z);
  if((cry_id%2)==0) {
    if(cry_id==96 || cry_id==98) {
    } else {
      v.RotateZ(-TMath::Pi()/3.);
    }
    //v.RotateX(TMath::Pi());
  } else {
    //v.RotateY(TMath::Pi());
  }
  //v.RotateX(TMath::Pi());
  //v.RotateY(TMath::Pi());
  //return CrystalToGlobal(cry_id,x,y,z);
  return CrystalToGlobal(cry_id,v.X(),v.Y(),v.Z());
}
/*******************************************************************************/
/* Defines which Gretina crystals are Neighbours used for Addback **************/
/*******************************************************************************/
void TGretina::SetGretNeighbours() {
  if(fNEIGHBOURSet){
    return;
  }

  std::string temp = getenv("GRUTSYS");
  temp.append("/libraries/TDetSystems/TGretina/gretina-pairs.dat");
  std::ifstream infile;
  infile.open(temp.c_str());
  if(!infile.is_open()) {
    std::cout<<"error with file"<<std::endl;
    return;
  }

  int i = 0;
  bool tmpB;
  while(infile.good()) {
    for(int j = 0; j < 124; j++) {
      infile >> tmpB;
      gretNeighbour[i][j] = tmpB;
    }
    i++;
  }

  infile.close();
  fNEIGHBOURSet = true;
  return;
}

/*******************************************************************************/
/* Returns vector for lab position of the centre of the crystal ****************/
/*******************************************************************************/
TVector3 TGretina::GetCrystalPosition(int cry_id) {
  SetCRMAT();

  TVector3 v;
  v.SetXYZ(0.0,0.0,0.0);
  for(int i=30;i<36;i++) {
    TVector3 a = GetSegmentPosition(cry_id,i);
    v.SetXYZ(v.X()+a.X(),v.Y()+a.Y(),v.Z()+a.Z());
  }
  v.SetXYZ(v.X()/6.0,v.Y()/6.0,v.Z()/6.0);
  return v;
}

/*******************************************************************************/
/* Copies hit ******************************************************************/
/*******************************************************************************/
void TGretina::Copy(TObject& obj) const {
  TDetector::Copy(obj);

  TGretina& gretina = (TGretina&)obj;
  gretina.gretina_hits = gretina_hits;
}

/*******************************************************************************/
/* Inserts Hit into TGretinaHit vector Copies hit ******************************/
/*******************************************************************************/
void TGretina::InsertHit(const TDetectorHit& hit){
  gretina_hits.emplace_back((TGretinaHit&)hit);
  fSize++;
}

/*******************************************************************************/
/* Unpacks GEB data and builds TGretinaHit *************************************/
/*******************************************************************************/
int TGretina::BuildHits(std::vector<TRawEvent>& raw_data){
  if(raw_data.size()<1)
    return Size();
  long smallest_time = 0x3fffffffffffffff;

  for(auto& event : raw_data){
    if(event.GetTimestamp()<smallest_time) smallest_time = event.GetTimestamp();
    TGretinaHit hit;
    TSmartBuffer buf = event.GetPayloadBuffer();
    hit.BuildFrom(buf);
    if(hit.GetPad() == 2 || hit.GetPad() == 3 || hit.GetPad() == 4 || hit.GetPad() == 6) continue;
    InsertHit(hit);
  }
  SetTimestamp(smallest_time);
  return Size();
}

/*******************************************************************************/
/* Converts Local Hit Position To Global Position ******************************/
/* x,y,z must be in cm *********************************************************/
/*******************************************************************************/
TVector3 TGretina::CrystalToGlobal(int cryId,Float_t x,Float_t y,Float_t z) {
  SetCRMAT();

  Int_t detectorPosition = cryId/4 - 1;
  Int_t crystalNumber    = cryId%4;

  double xl = ( (crmat[detectorPosition][crystalNumber][0][0] * x) +
                (crmat[detectorPosition][crystalNumber][0][1] * y) +
                (crmat[detectorPosition][crystalNumber][0][2] * z) +
                (crmat[detectorPosition][crystalNumber][0][3]) );

  double yl = ( (crmat[detectorPosition][crystalNumber][1][0] * x) +
                (crmat[detectorPosition][crystalNumber][1][1] * y) +
                (crmat[detectorPosition][crystalNumber][1][2] * z) +
                (crmat[detectorPosition][crystalNumber][1][3]) );

  double zl = ( (crmat[detectorPosition][crystalNumber][2][0] * x) +
                (crmat[detectorPosition][crystalNumber][2][1] * y) +
                (crmat[detectorPosition][crystalNumber][2][2] * z) +
                (crmat[detectorPosition][crystalNumber][2][3]) );

  return TVector3(xl, yl, zl);
}

/*******************************************************************************/
/* Basic Print Function ********************************************************/
/*******************************************************************************/
void TGretina::Print(Option_t *opt) const {
  printf(BLUE "GRETINA: size = %i" RESET_COLOR "\n",(int)Size());
  for(unsigned int x=0;x<Size();x++) {
    printf(DYELLOW);
    GetGretinaHit(x).Print(opt);
    printf(RESET_COLOR);
  }
  printf(BLUE "--------------------------------" RESET_COLOR "\n");
}

/*******************************************************************************/
/* Sorts Hits by energy ********************************************************/
/*******************************************************************************/
void TGretina::SortHits() {
  std::sort(gretina_hits.begin(),gretina_hits.end());
}

/*******************************************************************************/
/* Clears Hits *****************************************************************/
/*******************************************************************************/
void TGretina::Clear(Option_t *opt) {
  TDetector::Clear(opt);
  gretina_hits.clear();
  addback_hits.clear();
}


/*******************************************************************************/
/* Called by TGretinaHit::GetSpecOrder *****************************************/
/* Creates a map of crystals in theta/phi order and returns this order for *****/
/* given crystal id ************************************************************/
/*******************************************************************************/
int TGretina::GetSpecId(int cry_id, bool inverse) {
  SetSpecId();
  if(!inverse) {
    if(fThetaCryMap.find(cry_id) == fThetaCryMap.end()) return -1;
    else return fThetaCryMap[cry_id];
  } else {
    for (const auto& [key, value] : fThetaCryMap) {
      if(value == cry_id) return key;
    }
    return -1;
  }
}

/*******************************************************************************/
/* Uses the channel.cal file to determine which crystals are present and sorts */
/* them in order of Theta/Phi **************************************************/
/*******************************************************************************/
void TGretina::SetSpecId() {
  if(fSetSpec) {
    return;
  }

  std::vector<std::pair<int, TVector3>> crymap;
  for(auto const& it:  TChannel::fChannelMap) {
    if((strcmp(it.second->GetSystem(), "GRR") == 0) && (strcmp(it.second->GetCollectedCharge(), "N") == 0)) {
      int cryid = 4*((it.first & 0x1F00) >> 8) + ((it.first & 0x00C0) >> 6);
      TVector3 vec = TGretina::GetCrystalPosition(cryid);
      crymap.push_back(std::make_pair(cryid, vec));
    }
  }
  std::sort(crymap.begin(), crymap.end(), sortcrymap);
  int size = fThetaCryMap.size();
  for(int i = 0; i < (int)crymap.size(); i++) {
    fThetaCryMap.insert(std::make_pair(crymap.at(i).first, size));
    if(size < (int) fThetaCryMap.size()) size = (int)fThetaCryMap.size();
  }
  fSetSpec = true;
}
