#include "Globals.h"
#include "TS800Hit.h"

#include "chrono"
using namespace std::chrono;

/*******************************************************************************/
/* TTrigger ********************************************************************/
/* Stores data from ULM trigger module *****************************************/
/* Used To identift trigger source  purposes ***********************************/
/*******************************************************************************/
TTrigger::TTrigger() {
  Clear();
}

TTrigger::~TTrigger() {
}

void TTrigger::Copy(TObject &obj) const {
  TDetectorHit::Copy(obj);
  TTrigger &trig =(TTrigger&)obj;
  trig.fregistr  = fregistr;
  trig.fs800source      = fs800source;
  trig.fexternalsource1 = fexternalsource1;
  trig.fexternalsource2 = fexternalsource2;
  trig.fsecondarysource = fsecondarysource;
}

void TTrigger::Clear(Option_t *opt) {
  fregistr         = -1;
  fs800source      = -1;
  fexternalsource1 = -1;
  fexternalsource2 = -1;
  fsecondarysource = -1;
}

void TTrigger::Print(Option_t *opt) const {  }

/*******************************************************************************/
/* TIonChamber *****************************************************************/
/* Used to store data from S800 Focal Plane IC used for PID ********************/
/*******************************************************************************/
TIonChamber::TIonChamber() {
  Clear();
}

TIonChamber::~TIonChamber() {
}

/*******************************************************************************/
/* Fills IC Hit ****************************************************************/
/*******************************************************************************/
void TIonChamber::Set(int ch, int data){
  fChan.push_back(ch);
  fData.push_back(data);
}

/*******************************************************************************/
/* Gets average charge in IC calibrates data if .cal file is provided **********/
/*******************************************************************************/
float TIonChamber::GetAve(){
  float temp = 0.0;
  for(unsigned int x=0;x<fData.size();x++) {
    TChannel *c = TChannel::GetChannel(GetAddress(x));
    if (c){
      temp += c->CalEnergy(fData.at(x));
    } else{
      temp += fData.at(x);
    }
  }
  if(temp > 0) temp = temp/((float)fData.size());
  return temp;
}

/*******************************************************************************/
/* Gets sum of charge in IC calibrates data if .cal file is provided ***********/
/*******************************************************************************/
float TIonChamber::GetSum() const {
  float temp = 0.0;
  for(int x = 0; x < GetSize(); x++){
    TChannel *c = TChannel::GetChannel(GetAddress(x));
    if(c){
      temp += c->CalEnergy(GetData(x));
    } else{
      temp += GetData(x);
    }
  }
  return temp;
}

/*******************************************************************************/
/* Gets Energy Loss in IC uses CRDC posistion to account for entrance angle ****/
/*******************************************************************************/
 float TIonChamber::GetdE(TCrdc *crdc){
  float x   = crdc->GetDispersiveX();
  float y   = crdc->GetNonDispersiveY();

  return GetdE(x,y);
}

float TIonChamber::GetdE(double crdc_1_x, double crdc_1_y){
  float sum = GetAve();

  float xtilt  = GValue::Value("IC_DE_XTILT");
  float ytilt  = GValue::Value("IC_DE_YTILT");
  float x0tilt = GValue::Value("IC_DE_X0TILT");

  if (isnan(xtilt) || isnan(ytilt) || isnan(x0tilt)){
    std::cout << "Define  IC_DE_XTILT, IC_DE_YTILT, and IC_DE_X0TILT before using TIonChamber::GetdE()!\n";
  }
  sum += sum * ytilt * crdc_1_y;
  if (crdc_1_x < x0tilt){
    sum *= TMath::Exp(xtilt*(x0tilt-crdc_1_x));
  }
  return sum;
}

/*******************************************************************************/
/* Copies hit ******************************************************************/
/*******************************************************************************/
void TIonChamber::Copy(TObject &obj) const {
  TDetectorHit::Copy(obj);
  TIonChamber &ic =(TIonChamber&)obj;
  ic.fChan = fChan;
  ic.fData = fData;
}

/*******************************************************************************/
/* Clears hit ******************************************************************/
/*******************************************************************************/
void TIonChamber::Clear(Option_t *opt) {
  TDetectorHit::Clear(opt);
  fChan.clear();
  fData.clear();
}

/*******************************************************************************/
/* Blank Print Function ********************************************************/
/*******************************************************************************/
void TIonChamber::Print(Option_t *opt) const { }

/*******************************************************************************/
/* TCrdc ***********************************************************************/
/* Used to store data from the CRDCS *******************************************/
/* Gets CRDC1/2 X/y positions to calculate angles in S800 **********************/
/*******************************************************************************/
TCrdc::TCrdc() {
  Clear();
}

TCrdc::~TCrdc() {
}

/*******************************************************************************/
/* Find the pad number with the largest measured charge ************************/
/* Calibrates data if .cal file is provided ************************************/
/*******************************************************************************/
int TCrdc::GetMaxPad() const {
//  auto start = high_resolution_clock::now();
  if(!data.size()) return -1.0;
  int maxp = 0;
  float maxd = 0;
  for(unsigned int i = 0; i < data.size(); i++){
    if(!IsGoodSample(i)) {
std::cout << "BAD SAMPLE" << std::endl;
      continue;
    }
    TChannel *c = TChannel::GetChannel(GetAddress(i));
    double cal_data;
    if(c){
      cal_data = c->CalEnergy(data.at(i));
    } else{
      cal_data = (double)data.at(i);
    }
    if(cal_data > maxd) {
      maxp = channel.at(i);
      maxd = cal_data;
    }
  }
//  auto stop = high_resolution_clock::now();
//  auto duration = duration_cast<microseconds>(stop - start);
//  std::cout << "Time taken by MaxPad: " << duration.count() << " microseconds" << std::endl;
  return maxp;
}

/*******************************************************************************/
/* Fast Function Find the pad number with the largest measured charge **********/
/* WIP - Do not use ************************************************************/
/*******************************************************************************/
int TCrdc::GetMaxPadFast() const {
//  auto start = high_resolution_clock::now();
  //WIP
  //Attempting to speed up function
  if(!data.size()) return -1.0;
  int dis = std::distance(data.begin(), std::max_element(data.begin(), data.end()));
  int maxp = channel.at(dis);

//  auto stop = high_resolution_clock::now();
//  auto duration = duration_cast<microseconds>(stop - start);
//  std::cout << "Time taken by Fast function: " << duration.count() << " microseconds" << std::endl;
  return maxp;
}

/*******************************************************************************/
/* Equivalent to GetMaxPad() but returns the charge not pad number *************/
/* Calibrates data if .cal file is provided ************************************/
/*******************************************************************************/
int TCrdc::GetMaxPadSum() const{
  if(!data.size()) return -1.0;
  float maxd = 0;
  for(unsigned int i = 0; i < data.size(); i++){
    if(!IsGoodSample(i)) {
      continue;
    }
    TChannel *c = TChannel::GetChannel(GetAddress(i));
    double cal_data;
    if(c){
      cal_data = c->CalEnergy(data.at(i));
    } else{
      cal_data = (double)data.at(i);
    }
    if(cal_data > maxd) {
      maxd = cal_data;
    }
  }
  return maxd;
}
/*******************************************************************************/
/* Not Really sure the purpose of this function ********************************/
/*******************************************************************************/
int TCrdc::GetWidth() {
  if(GetSize()<2)
    return 0;
  return sample.back()-sample.front()+1;
}

/*******************************************************************************/
/* Copies hit ******************************************************************/
/*******************************************************************************/
void TCrdc::Copy(TObject &obj) const {
  TDetectorHit::Copy(obj);
  TCrdc &c = (TCrdc&)obj;
  c.fId      = fId;
  c.channel  = channel;
  c.sample   = sample;
  c.data     = data;
  c.anode    = anode;
  c.time     = time;
  c.fDrifttime = fDrifttime;
}

/*******************************************************************************/
/* Clears hit ******************************************************************/
/*******************************************************************************/
void TCrdc::Clear(Option_t *opt) {
  TDetectorHit::Clear(opt);
  fId   = -1;
  anode = 0;
  time  = 0;
  fDrifttime = 0;
  channel.clear();
  sample.clear();
  data.clear();
  has_cached_dispersive_x = false;
  cached_dispersive_x = std::sqrt(-1);
}

/*******************************************************************************/
/* Checks if good sample from digitizer used by GetMaxPad/GetDispersiveX *******/
/*******************************************************************************/
bool TCrdc::IsGoodSample(int i) const {
  return true;
  if (i == 0 && data.size()>1) {
    return channel.at(i) == channel.at(i+1);
  } else if(i == (int)data.size()-1 && data.size()>2) {
    return channel.at(i) == channel.at(i-1);
  } else if(data.size()>2) {
    return (channel.at(i) == channel.at(i-1) ||
	    (channel.at(i) == channel.at(i+1)));
  }
  return false;
}

/*******************************************************************************/
/* Calculate X position of CRDC - Uses the Mean of GetMaxPad +/- 7 strips ******/
/*******************************************************************************/
float TCrdc::GetDispersiveX() const{
  int maxpad = GetMaxPad();

  if (maxpad ==-1){
    has_cached_dispersive_x = true;
    return sqrt(-1);
  }
  float x_slope = sqrt(-1);
  float x_offset = sqrt(-1);
  if(fId == 0) {
    x_slope = GValue::Value("CRDC1_X_SLOPE");
    x_offset = GValue::Value("CRDC1_X_OFFSET");
  } else {
    x_slope = GValue::Value("CRDC2_X_SLOPE");
    x_offset = GValue::Value("CRDC2_X_OFFSET");
  }
  if(std::isnan(x_slope))
    x_slope = 1.0;
  if(std::isnan(x_offset))
    x_offset = 0.0;

  const int GRAVITY_WIDTH = 14;//determines how many pads one uses in averaging
  const int NUM_PADS = 224;
//  const int NUM_PADS = 480;
  int lowpad = maxpad - GRAVITY_WIDTH/2;
  int highpad = lowpad + GRAVITY_WIDTH;
  if (lowpad < 0){
    lowpad = 0;
  }
  if (highpad >= NUM_PADS){
    highpad = NUM_PADS-1;
  }

  double datasum = 0;
  double weighted_sum = 0;
  for(int i=0;i<GetSize();i++) {
    if((channel.at(i) < lowpad) || (channel.at(i)>highpad) || !IsGoodSample(i)) {
      continue;
    }

    TChannel *c = TChannel::GetChannel(GetAddress(i));
    double cal_data;
    if (c){
      cal_data = c->CalEnergy(GetData(i));
    } else{
      cal_data = (double)GetData(i);
    }
    datasum += cal_data;
    weighted_sum += channel.at(i)*cal_data;
  }

  // + 0.5 so that we take the middle of the pad, not the left edge.
  double mean_chan = weighted_sum/datasum + 0.5;
  double output = (mean_chan*x_slope+x_offset);

  has_cached_dispersive_x = true;
  cached_dispersive_x = output;
  return output;
}

/*******************************************************************************/
/* Calculate Y position of CRDC ************************************************/
/*******************************************************************************/
float TCrdc::GetNonDispersiveY() {
  if (GetMaxPad() == -1){
    return sqrt(-1);
  }

  float y_slope = sqrt(-1);
  float y_offset = sqrt(-1);
  if(fId==0) {
    y_slope = GValue::Value("CRDC1_Y_SLOPE");
    y_offset = GValue::Value("CRDC1_Y_OFFSET");
  } else if(fId==1) {
    y_slope = GValue::Value("CRDC2_Y_SLOPE");
    y_offset = GValue::Value("CRDC2_Y_OFFSET");
  }

  if(std::isnan(y_slope)) y_slope = 1.0;
  if(std::isnan(y_offset)) y_offset = 0.0;
  float tmp = ((float)time) * y_slope + y_offset;
  return tmp;
}

/*******************************************************************************/
/* Calculate X position of CRDC - Uses the Mean of GetMaxPad +/- 7 strips ******/
/* Provides maxpad as argument to reduce function calls ************************/
/*******************************************************************************/
float TCrdc::GetDispersiveX(int maxpad) const{
  if (maxpad ==-1){
    has_cached_dispersive_x = true;
    return sqrt(-1);
  }
  float x_slope = sqrt(-1);
  float x_offset = sqrt(-1);
  if(fId == 0) {
    x_slope = GValue::Value("CRDC1_X_SLOPE");
    x_offset = GValue::Value("CRDC1_X_OFFSET");
  } else {
    x_slope = GValue::Value("CRDC2_X_SLOPE");
    x_offset = GValue::Value("CRDC2_X_OFFSET");
  }
  if(std::isnan(x_slope))
    x_slope = 1.0;
  if(std::isnan(x_offset))
    x_offset = 0.0;

  const int GRAVITY_WIDTH = 14;//determines how many pads one uses in averaging
  const int NUM_PADS = 224;
  int lowpad = maxpad - GRAVITY_WIDTH/2;
  int highpad = lowpad + GRAVITY_WIDTH;
  if (lowpad < 0){
    lowpad = 0;
  }
  if (highpad >= NUM_PADS){
    highpad = NUM_PADS-1;
  }

  double datasum = 0;
  double weighted_sum = 0;
  for(int i = 0; i < GetSize(); i++) {
    if((channel.at(i) < lowpad) || (channel.at(i)>highpad) || !IsGoodSample(i)) {
      continue;
    }

    TChannel *c = TChannel::GetChannel(GetAddress(i));
    double cal_data;
    if (c){
      cal_data = c->CalEnergy(GetData(i));
    } else{
      cal_data = (double)GetData(i);
    }
    datasum += cal_data;
    weighted_sum += channel.at(i)*cal_data;
  }

  // + 0.5 so that we take the middle of the pad, not the left edge.
  double mean_chan = weighted_sum/datasum + 0.5;
  double output = (mean_chan*x_slope+x_offset);

  has_cached_dispersive_x = true;
  cached_dispersive_x = output;
  return output;
}

/*******************************************************************************/
/* Calculate Y position of CRDC ************************************************/
/* Provides maxpad as argument to reduce function calls ************************/
/*******************************************************************************/
float TCrdc::GetNonDispersiveY(int maxpad) {
  if (maxpad ==-1){
    return sqrt(-1);
  }

  float y_slope = sqrt(-1);
  float y_offset = sqrt(-1);
  if(fId==0) {
    y_slope = GValue::Value("CRDC1_Y_SLOPE");
    y_offset = GValue::Value("CRDC1_Y_OFFSET");
  } else if(fId==1) {
    y_slope = GValue::Value("CRDC2_Y_SLOPE");
    y_offset = GValue::Value("CRDC2_Y_OFFSET");
  }

  if(std::isnan(y_slope)) y_slope = 1.0;
  if(std::isnan(y_offset)) y_offset = 0.0;
  float tmp = ((float)time)*y_slope+y_offset;
  return tmp;
}

/*******************************************************************************/
/* Basic Print Function ********************************************************/
/*******************************************************************************/
void TCrdc::Print(Option_t *opt) const {

  if(!strcmp(opt,"all")) {
    std::cout << "Channel \t Sample \t Data" << std::endl;
    for(unsigned int i = 0; i < data.size(); i++){
      std::cout << GetChannel(i) << " \t " << GetSample(i) << " \t " << GetData(i) << std::endl;
    }
  }
}

/*******************************************************************************/
/* THodoscope ******************************************************************/
/* Used to store data from the Hodoscope ***************************************/
/*******************************************************************************/
THodoscope::THodoscope() { Clear(); }
THodoscope::~THodoscope() { Clear(); }
THodoscope::THodoscope(const THodoscope &hodoscope) {
  hodoscope.Copy(*this);
}

/*******************************************************************************/
/* Copies hit ******************************************************************/
/*******************************************************************************/
void THodoscope::Copy(TObject &obj) const {
  TDetectorHit::Copy(obj);
  THodoscope &hodo = (THodoscope&)obj;
  hodo.hodo_hits = hodo_hits;
}

/*******************************************************************************/
/* Clears hit ******************************************************************/
/*******************************************************************************/
void THodoscope::Clear(Option_t *opt){
  TDetectorHit::Clear(opt);
  hodo_hits.clear();
}

/*******************************************************************************/
/* Basic Print Function ********************************************************/
/*******************************************************************************/
void THodoscope::Print(Option_t *opt) const {
  printf("Number of Hits: %zu\n",GetSize());
  for(unsigned int i=0;i<GetSize();i++) {
    printf("\t");
    GetHodoHit(i).Print();
  }
  printf("---------------------------------------\n");
}

/*******************************************************************************/
/* Inserts hits ****************************************************************/
/*******************************************************************************/
void THodoscope::InsertHit(const TDetectorHit& hit){
  hodo_hits.emplace_back((THodoHit&)hit);
}

/*******************************************************************************/
/* Functions to call Hodoscope hits ********************************************/
/*******************************************************************************/
TDetectorHit& THodoscope::GetHit(int i){
  return hodo_hits.at(i);
}

const THodoHit& THodoscope::GetHodoHit(int i) const {
  return hodo_hits.at(i);
}

THodoHit& THodoscope::GetHodoHit(int i) {
  return hodo_hits.at(i);
}

/*******************************************************************************/
/* Declare hodoscope hit vector ************************************************/
/*******************************************************************************/
THodoHit::THodoHit(const THodoHit &hit){
  hit.Copy(*this);
}

/*******************************************************************************/
/* Copies hit vector ***********************************************************/
/*******************************************************************************/
void THodoHit::Copy(TObject &obj) const {
  TDetectorHit::Copy(obj);
  THodoHit &hit = (THodoHit&)obj;
  hit.fChannel = fChannel;
}

/*******************************************************************************/
/* Basic Print Function ********************************************************/
/*******************************************************************************/
void THodoHit::Print(Option_t *opt) const {
  printf("Channel: %d Charge: %d\n", GetChannel(), GetCharge());
}

/*******************************************************************************/
/* Clears hodoscope hits *******************************************************/
/*******************************************************************************/
void THodoHit::Clear(Option_t *opt){
  TDetectorHit::Clear(opt);
  fChannel = sqrt(-1);
}

/*******************************************************************************/
/* TMToF ***********************************************************************/
/* Used to store data from Mestec MTDC-32 **************************************/
/* Each channel of MTDC is stored in an individual vector fE1Up, fObj etc. *****/
/* Used primarily for ToF for PID purposes *************************************/
/*******************************************************************************/
TMTof::TMTof() { Clear(); }

TMTof::~TMTof() { }

TMTof::TMTof(const TMTof &mtof) {
  mtof.Copy(*this);
}


/*******************************************************************************/
/* Copy TMTof hit **************************************************************/
/*******************************************************************************/
void TMTof::Copy(TObject &obj) const {
  TDetectorHit::Copy(obj);

  TMTof &mtof = ((TMTof&)obj);

  mtof.fE1Up      = fE1Up;
  mtof.fE1Down    = fE1Down;
  mtof.fXfp       = fXfp;
  mtof.fObj       = fObj;
  mtof.fCrdc1Anode = fCrdc1Anode;
  mtof.fCrdc2Anode = fCrdc2Anode;
  mtof.fRf        = fRf;
  mtof.fHodoscope = fHodoscope;
  mtof.fRef       = fRef;
}

/*******************************************************************************/
/* Clears MTDC vectors *********************************************************/
/*******************************************************************************/
void TMTof::Clear(Option_t *opt) {

  fE1Up.clear();
  fE1Down.clear();
  fXfp.clear();
  fObj.clear();
  fCrdc1Anode.clear();
  fCrdc2Anode.clear();
  fRf.clear();
  fHodoscope.clear();
  fRef.clear();
}


/*******************************************************************************/
/* Generic correlated Time of flight functions between any two MTDC's channels */
/* ch1/ch2 - channels numbers of MTDC for time difference **********************/
/* target - used to line up select specific MToF of peak ***********************/
/* shift - used to shift MToF peak in the event peaks shift run to run *********/
/* If no arguments for "target" are given returns difference between the first */
/* value of MTDC vectors selected by ch1/ch2 ***********************************/
/* Replaces GetCorrelatedbjE1(), GetCorrelatedXFPE1() **************************/
/*******************************************************************************/
double TMTof::GetCorrelatedTof(int ch1, int ch2, double target, double shift) const{
  double fCorr = -1;
  std::vector<unsigned short> refvec1 = GetMTofVector(ch1);
  std::vector<unsigned short> refvec2 = GetMTofVector(ch2);
  if(refvec1.size() && refvec2.size()){
    if(std::isnan(target)) return (refvec1.at(0) - refvec2.at(0));
    fCorr = std::numeric_limits<double>::max();
    for(size_t i = 0; i < refvec1.size(); i++) {
      for (size_t j = 0; j < refvec2.size(); j++){
        double newvalue = refvec1.at(i) - refvec2.at(j);
	if(!std::isnan(shift)) {
	  newvalue += shift;
        }
        if(std::abs(target - newvalue) < std::abs(target - fCorr)) {
          fCorr = newvalue;
        }
      }
    }
  } else return sqrt(-1);
  return fCorr;
}

/*******************************************************************************/
/* Returns one of the MTDC vectors based on channel number *********************/
/* channel - Channel number of MTDC ********************************************/
/* If no arguments are given returns E1Up **************************************/
/*******************************************************************************/
std::vector<unsigned short> TMTof::GetMTofVector(int channel) const {
  switch(channel) {
    case 0  : //E1 up
      return  fE1Up;
      break;
    case 1  : //E1 down
      return  fE1Down;
      break;
    case 2  : //Xfp
      return  fXfp;
      break;
    case 3  : //Obj
      return  fObj;
      break;
    case 5  : //Rf
      return  fRf;
      break;
    case 6  : //CRDC1 anode
      return  fCrdc1Anode;
      break;
    case 7  : //CRDC2 anode
      return  fCrdc2Anode;
      break;
    case 12 : //Hodoscope
      return  fHodoscope;
      break;
    case 15 : //Refrence
      return  fRef;
      break;
    default : //E1 Up
      return  fE1Up;
      break;
  }
  return fE1Up;
}

/*******************************************************************************/
/* Generic correlated Time of flight functions between any two MTDC's channels */
/* ch1/ch2 - Detector name in MTDC for time difference *************************/
/* target used to line up select specific MToF of peak *************************/
/* shift used to shift MToF peak in the event peaks shift run to run ***********/
/* If no arguments for "target" are given returns difference between the first */
/* value of MTDC vectors selected by ch1/ch2 ***********************************/
/* Replaces GetCorrelatedbjE1(), GetCorrelatedXFPE1() **************************/
/*******************************************************************************/
double TMTof::GetCorrelatedTof(std::string ch1, std::string ch2, double target, double shift) const{
  double fCorr = -1;
  std::vector<unsigned short> refvec1 = GetMTofVectorFromString(ch1);
  std::vector<unsigned short> refvec2 = GetMTofVectorFromString(ch2);

  if(refvec1.size() && refvec2.size()){
    if(std::isnan(target)) return (refvec1.at(0) - refvec2.at(0));
    fCorr = std::numeric_limits<double>::max();
    for(size_t i = 0; i < refvec1.size(); i++) {
      for (size_t j = 0; j < refvec2.size(); j++){
        double newvalue = refvec1.at(i) - refvec2.at(j);
	if(!std::isnan(shift)) {
	  newvalue += shift;
        }
        if(std::abs(target - newvalue) < std::abs(target - fCorr)) {
          fCorr = newvalue;
        }
      }
    }
  } else return sqrt(-1);
  return fCorr;
}

/*******************************************************************************/
/* Returns one of the MTDC vectors based on channel name ***********************/
/* channel - Detector name in MTDC *********************************************/
/* If no arguments are given or recognised returns E1Up ************************/
/*******************************************************************************/
std::vector<unsigned short> TMTof::GetMTofVectorFromString(std::string vecname) const {
  if(vecname == "E1Up")        return fE1Up;
  else if(vecname == "E1Down") return  fE1Down;
  else if(vecname == "Xfp")    return  fXfp;
  else if(vecname == "Obj")    return  fObj;
  else if(vecname == "Rf")    return  fRf;
  else if(vecname == "Crdc1")  return  fCrdc1Anode;
  else if(vecname == "Crdc2")  return  fCrdc2Anode;
  else if(vecname == "Hodo")   return  fHodoscope;
  else if(vecname == "Ref")    return  fRef;
  else {
    std::cout << "Vector Name Not Recognised using fE1Up" << std::endl;
    std::cout << "Recognised Names: E1Up E1Down Xfp Obj Rf Crdc1 Crdc2 Hodo Ref" << std::endl;
    return  fE1Up;
  }
}

/*******************************************************************************/
/* Print function - no functionality *******************************************/
/*******************************************************************************/
void TMTof::Print(Option_t *opt) const {    }
