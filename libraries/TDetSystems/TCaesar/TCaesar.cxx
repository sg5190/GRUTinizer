#include "TCaesar.h"

#include <algorithm>
#include <iostream>
#include <fstream>

#include "DDASDataFormat.h"
#include "DDASBanks.h"
#include "TNSCLEvent.h"
#include "TRawEvent.h"
#include "TRawBanks.h"
#include "TChannel.h"
#include "TGRUTOptions.h"
#include "TVector3.h"


#include "chrono"
using namespace std::chrono;

std::map<int,TCaesar::coord> TCaesar::detector_positions;
bool TCaesar::fNEIGHBOURSet = false;
bool TCaesar::caesarNeighbour[192][192];

/*******************************************************************************/
/* TCaesar ****************************************************************/
/*******************************************************************************/
TCaesar::TCaesar(){ }

TCaesar::~TCaesar(){ }

/*******************************************************************************/
/* Main function for making addback hits based on GRETINA algorithm ************/
/* SortDepth determines how many combinations of events to consider default=4 **/
/* AddbackDepth classifies events as, n0 (0), n1 (1), n2 (2) or ng (3) *********/
/*******************************************************************************/
void TCaesar::BuildAddback(int SortDepth) const {
  //Don't rebuild events
  if( addback_hits.size() > 0 || caesar_hits.size() == 0) {
    return;
  }

  std::vector<TCaesarHit> temp_hits;
  //Get rid of low energy, overlow and underflow events
  for(unsigned int i = 0; i < caesar_hits.size(); i++) {
    if(caesar_hits.at(i).GetEnergy() > 50 && !caesar_hits.at(i).IsOverflow() && !caesar_hits.at(i).IsUnderflow()) temp_hits.push_back(caesar_hits.at(i));
  }

  //sort so that the first hit has the greatest energy
  //this way we can loop through i,j with i < j and know that
  //any hit with higher energy cannot be an addback to one with lower energy
  std::sort(temp_hits.begin(), temp_hits.end(),
    [](const TCaesarHit& a, const TCaesarHit& b) {
    return a.GetEnergy() > b.GetEnergy();
  });

  //vector used to store hit indices when crystals are pairs
  std::vector<int> paired;
  //loop through every hit
  for(unsigned int i=0; i < temp_hits.size(); i++) {
    paired.clear();
    TCaesarHit &current_hit = temp_hits[i];
    //only pick unqiue pairs of hits
    if(SortDepth > 1) {
      for(unsigned int j=i+1; j < temp_hits.size(); j++) {
        if (IsNeighbour(current_hit,temp_hits[j])){
          paired.push_back(j);
          //check every hit k to see if it is a neighbour with j
          if(SortDepth > 2) {
            for (unsigned int k=i; k < temp_hits.size(); k++){
              if( (i == k) || (j == k) ) continue;
              if (IsNeighbour(temp_hits[j],temp_hits[k])){
                paired.push_back(k);
                //Edge cases, requires at least 4 crystals in a line
                if(SortDepth > 3) {
                  for (unsigned int l=i; l < temp_hits.size(); l++){
                    if( (i == l) || (j == l) || (k == l)) continue;
                    if (IsNeighbour(temp_hits[k],temp_hits[l])) {
                      paired.push_back(l);
                      //Really Edge cases, requires at least 5 crystals in a line
                      if(SortDepth > 4) {
                        for (unsigned int m=i; m < temp_hits.size(); m++){
                          if( (i == m) || (j == m) || (k == m) || (l == m)) continue;
                          if(IsNeighbour(temp_hits[l],temp_hits[m])) {
                            paired.push_back(m);
                            //Extreme Edge cases, requires at least 6 crystals in a line
                            if(SortDepth > 5) {
                              for (unsigned int n=i; n < temp_hits.size(); n++){
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
/* Copies hit ******************************************************************/
/*******************************************************************************/
void TCaesar::Copy(TObject& obj) const {
  TDetector::Copy(obj);

  TCaesar& caesar = (TCaesar&)obj;
  caesar.caesar_hits = caesar_hits;
 }

/*******************************************************************************/
/* Clear hit *******************************************************************/
/*******************************************************************************/
void TCaesar::Clear(Option_t* opt){
  TDetector::Clear(opt);
  caesar_hits.clear();
}

/*******************************************************************************/
/* Functions to call DDAS hits *************************************************/
/*******************************************************************************/
TCaesarHit& TCaesar::GetCaesarHit(int i){
  return caesar_hits.at(i);
}

TDetectorHit& TCaesar::GetHit(int i){
  return caesar_hits.at(i);
}

/*******************************************************************************/
/* Unpacks  data and builds TCaesar events *************************************/
/*******************************************************************************/
int TCaesar::BuildHits(std::vector<TRawEvent>& raw_data) {
  //Loop over raw data
  for(auto& event : raw_data){
    TSmartBuffer buf = event.GetPayloadBuffer();
    nsize = buf.GetSize();
    const char *data = buf.GetData();

    const TRawEvent::VMUSB_Header* vmusb_header = (TRawEvent::VMUSB_Header*)data;
    data += sizeof(TRawEvent::VMUSB_Header);
    data += 2;
    const  TRawEvent::VME_Timestamp* vme_timestamp = (TRawEvent::VME_Timestamp*)data;
    data += sizeof(TRawEvent::VME_Timestamp);
    long timestamp = vme_timestamp->ts2() * 8.7184750e-7 * 1e9;
    num_packets = vmusb_header->size()/2 - 3;

    int adcN = 0;
    int tdcN = 0;
    for(int i = 0; i < num_packets; i++) {
      TCaesarHit *hit = NULL;
      if (tdcN < 6) {
        const TRawEvent::Mesy_Word* Mesy_Word = (TRawEvent::Mesy_Word*) data;
        data += sizeof(TRawEvent::Mesy_Word);
        if(Mesy_Word->isData()) {
          const TRawEvent::M_TDC_Data* Mt = (TRawEvent::M_TDC_Data*) Mesy_Word;
          if(Mt->isTrig()) continue;
          unsigned int address = ((80 << 24) + (tdcN << 8) + (Mt->Chan()));
          hit = FindHit(address, true);
	  hit->SetTime(Mt->Time());
          hit->SetChannelNum(tdcN*32 + Mt->Chan());
	  hit->SetTimestamp(timestamp);
        }
        if(Mesy_Word->isALLF()) {
          tdcN++;
        }
      } else if(adcN < 6) {
        const TRawEvent::CAEN_DataPacket* Caen_Word = (TRawEvent::CAEN_DataPacket*) data;
        data += sizeof(TRawEvent::CAEN_DataPacket);
        if(Caen_Word->isData()) {
	  //As per CAEN V785 manual sec 2.2 with sliding scale enabled adc values
          //of 3841 - 4095 are not correct so skip these
//	  if(Caen_Word->adcvalue() > 3840) continue;
          unsigned int address = ((80 << 24) + (adcN << 8) + (Caen_Word->channel_num()));
	  hit = FindHit(address, false);
	  hit->SetCharge(Caen_Word->adcvalue());
	  hit->SetChannelNum(adcN*32 + Caen_Word->channel_num());
	  hit->SetABDepth(-1);
	  hit->SetTimestamp(timestamp);
	  hit->SetOverflow(Caen_Word->overflow());
	  hit->SetUnderflow(Caen_Word->underflow());
        }
        if(Caen_Word->isEnd()) adcN++;;
      }
    }
    fSize++;
  }
  return GetSize();
}

/*******************************************************************************/
/* Find Hit by Address used to match ADC/TDC ***********************************/
/* Creates a new CAESAR hit for ADCs and when called by TDC unpacker searches **/
/* for existing ADC hits if not is found a new hit is created ******************/
/*******************************************************************************/
TCaesarHit *TCaesar::FindHit(int address, bool search) {
  if(search) {
    for(unsigned int i=0;i<caesar_hits.size();i++) {
      if(caesar_hits.at(i).GetAddress()==address) {
        return &caesar_hits.at(i);
      }
    }
  }
  caesar_hits.emplace_back();
  TCaesarHit* hit = &caesar_hits.back();
  hit->SetAddress(address);
  return hit;
}

/*******************************************************************************/
/* Returns real position (in mm) from channel number ***************************/
/*******************************************************************************/
TVector3 TCaesar::GetPosition(int chNum) {
  LoadDetectorPositions();
  static int lines_displayed = 0;
  if(chNum < 0 || chNum > 192) {
    if(lines_displayed < 100) {
      std::cout << "No positions loaded for CAESAR Ch " << chNum << std::endl;
    }
    lines_displayed++;
    return TVector3(std::sqrt(-1),std::sqrt(-1),std::sqrt(-1));
  }
  coord &cc = detector_positions[chNum];
  return cc.pos;
}


/*******************************************************************************/
/* Returns ring number from channel number *************************************/
/*******************************************************************************/
int TCaesar::GetRing(int chNum) {
  LoadDetectorPositions();
  static int lines_displayed = 0;
  if(chNum < 0 || chNum > 192) {
    if(lines_displayed < 1000) {
      std::cout << "No ring positions loaded for CAESAR Ch " << chNum << std::endl;
    }
    lines_displayed++;
    return -1;
  }
  coord &cc = detector_positions[chNum];
  return cc.ring;
}

/*******************************************************************************/
/* Returns detector position in Ring from channel number ***********************/
/*******************************************************************************/
int TCaesar::GetDetectorNumber(int chNum) {
  LoadDetectorPositions();
  static int lines_displayed = 0;
  if(chNum < 0 || chNum > 192) {
    if(lines_displayed < 100) {
      std::cout << "No detector positions loaded for CAESAR Ch " << chNum << std::endl;
    }
    lines_displayed++;
    return -1;
  }
  coord &cc = detector_positions[chNum];
  return cc.det;
}

/*******************************************************************************/
/* Defines which Gretina crystals are Neighbours used for Addback **************/
/*******************************************************************************/
void TCaesar::SetCaesarNeighbours() {
  if(fNEIGHBOURSet){
    return;
  }

  std::string temp = getenv("GRUTSYS");
  std::string filename = std::string(getenv("GRUTSYS")) + "/config/CAESAR_Neighbours.txt";
  std::ifstream infile;
  infile.open(filename.c_str());
  if(!infile.is_open()) {
    std::cout<<"error with file"<<std::endl;
    return;
  }

  int i = 0;
  bool tmpB;
  while(infile.good()) {
    for(int j = 0; j < 192; j++) {
      infile >> tmpB;
      caesarNeighbour[i][j] = tmpB;
    }
    i++;
  }

  infile.close();
  fNEIGHBOURSet = true;
  return;
}

/*******************************************************************************/
/* Load in Detector Positions to calculate positions for doppler corections ****/
/*******************************************************************************/
void TCaesar::LoadDetectorPositions() {
  static bool loaded = false;
  if(loaded){
    return;
  }
  loaded = true;
  std::string filename = std::string(getenv("GRUTSYS")) + "/config/CAESAR_DetPos.txt";

  //Read the locations from file.
  std::ifstream infile(filename);
  if(!infile){  std::cout << "CAESAR detector position file \"" << filename << "\"" << " does not exist, skipping" << std::endl;
    return;
  }

  std::string line;
  while (std::getline(infile,line)) {
    //Parse the line
    int ringN, detN, adcN, chN;
    double x,y,z;
    int extracted = sscanf(line.c_str(),"%i\t%i\t%i\t%i\t%lf\t%lf\t%lf", &ringN, &detN, &adcN, &chN, &x, &y, &z);
    if (extracted!=7) {
      continue;
    }
    int tmpCh = adcN * 32 + chN;
    detector_positions[tmpCh].ring = ringN;
    detector_positions[tmpCh].det = detN;
    TVector3 pvec = TVector3(x, y, z);
    detector_positions[tmpCh].pos = pvec;
  }
}

/*******************************************************************************/
/* Required by TDetectorFactory ************************************************/
/*******************************************************************************/
void TCaesar::InsertHit(const TDetectorHit& hit) {
  caesar_hits.emplace_back((TCaesarHit&)hit);
  fSize++;
}
