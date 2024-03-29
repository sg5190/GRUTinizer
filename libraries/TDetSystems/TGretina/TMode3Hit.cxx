#include "TMode3Hit.h"
#include "TRandom.h"
#include "TGEBEvent.h"
#include "TGRUTOptions.h"
#include "TChannel.h"

ClassImp(TMode3Hit)

bool TMode3Hit::fExtractWaves = true;

TMode3Hit::TMode3Hit(){ }

TMode3Hit::~TMode3Hit() { }

/*******************************************************************************/
/* Builds the TMode3Hit from TSmartBuffer **************************************/
/*******************************************************************************/
void TMode3Hit::BuildFrom(TSmartBuffer& buf){
  Clear();

  bool read_waveform = TGRUTOptions::Get()->ExtractWaves();

  auto header = (TRawEvent::GEBMode3Head*)buf.GetData();
  buf.Advance(sizeof(TRawEvent::GEBMode3Head));

  board_id = header->board_id;
  SetAddress((1<<24) + board_id);

  auto data   = (TRawEvent::GEBMode3Data*)buf.GetData();
  buf.Advance(sizeof(TRawEvent::GEBMode3Data));

  led = data->GetLed();
  SetCharge(data->GetEnergy(*header));

  dt1  = data->GetDeltaT1();
  dt2  = data->GetDeltaT2();

  charge0  = data->GetEnergy0(*header);
  charge1  = data->GetEnergy1(*header);
  charge2  = data->GetEnergy2(*header);

  cfd = data->GetCfd();

  size_t wave_bytes = header->GetLength()*4 - sizeof(*header) + 4 - sizeof(*data);
  if(read_waveform){
    size_t wavesize = wave_bytes/sizeof(short);
    waveform.resize(wavesize);

    memcpy((char*)&waveform[0], buf.GetData(), wave_bytes);

    for(unsigned int i=0; i<wavesize; i+=2){
      short tmp      = TRawEvent::SwapShort(waveform[i+1]);
      waveform[i+1] = TRawEvent::SwapShort(waveform[i]);
      waveform[i]   = tmp;
    }
  }
  buf.Advance(wave_bytes);
}

/*******************************************************************************/
/* Copies hit ******************************************************************/
/*******************************************************************************/
void TMode3Hit::Copy(TObject& obj) const {
  TDetectorHit::Copy(obj);

  TMode3Hit& mode3 = (TMode3Hit&)obj;
  mode3.board_id = board_id;
  mode3.led      = led;
  mode3.cfd      = cfd;
  mode3.waveform = waveform;
}


/*******************************************************************************/
/* Basic Print Function ********************************************************/
/*******************************************************************************/
void TMode3Hit::Print(Option_t *opt) const {
  printf("hole[%03i] xtal[%i] seg[%03i]:   %i\n",GetHole(),GetCrystal(),GetSegmentId(),GetCharge0());
}

/*******************************************************************************/
/* Clears Hit ******************************************************************/
/*******************************************************************************/
void TMode3Hit::Clear(Option_t *opt) {
  TDetectorHit::Clear(opt);
  board_id = -1;
  led      = -1;
  cfd      = -1;
  dt1      =  0xffff;
  dt2      =  0xffff;
  charge0  = -1;
  charge1  = -1;
  charge2  = -1;

  waveform.clear();
}

/*******************************************************************************/
/* Returns the average of mode3 wave *******************************************/
/*******************************************************************************/
double TMode3Hit::AverageWave(int samples) const {
  if(waveform.size() == 0) {
    return 0.0;
  }
  if(samples < 0 ||
     samples > int(waveform.size())) {
    samples = waveform.size();
  }
  double sum = 0.0;
  for(int i=0;i<samples;i++) {
    sum += waveform[i];
  }
  return sum / ((double)samples);
}

/*******************************************************************************/
/* Get Charges from other Core Channels? ***************************************/
/*******************************************************************************/
double TMode3Hit::GetEnergy0() const {
  double energy;
  TChannel* chan = TChannel::GetChannel(Address());
  if(!chan){
    energy = charge0/128 + gRandom->Uniform();
    //return Charge() + gRandom->Uniform();
  } else {
    energy = chan->CalEnergy(charge0/128, fTimestamp);
  }
  return energy;
}

double TMode3Hit::GetEnergy1() const {
  double energy;
  TChannel* chan = TChannel::GetChannel(Address());
  if(!chan){
    energy = charge1/128 + gRandom->Uniform();
    //return Charge() + gRandom->Uniform();
  } else {
    energy = chan->CalEnergy(charge1/128, fTimestamp);
  }
  return energy;
}

double TMode3Hit::GetEnergy2() const {
  double energy;
  TChannel* chan = TChannel::GetChannel(Address());
  if(!chan){
    energy = charge2/128 + gRandom->Uniform();
    //return Charge() + gRandom->Uniform();
  } else {
    energy = chan->CalEnergy(charge2/128, fTimestamp);
  }
  return energy;
}
