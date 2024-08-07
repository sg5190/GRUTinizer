#include "TUnpackingLoop.h"

#include <chrono>
#include <thread>

//#include "RawDataQueue.h"
//#include "RawDataQueue.h"
//#include "TDetectorEnv.h"
#include "TGEBEvent.h"
#include "TGRUTOptions.h"
#include "TNSCLEvent.h"

#include "TDetectorEnv.h"
#include "TUnpackedEvent.h"

#include "chrono"
using namespace std::chrono;

//#include "TMode3.h"

ClassImp(TUnpackingLoop)

TUnpackingLoop *TUnpackingLoop::Get(std::string name) {
  if(name.length() == 0) {
    name = "unpack_loop";
  }
  TUnpackingLoop *loop = (TUnpackingLoop*)StoppableThread::Get(name);
  if(!loop) {
    loop = new TUnpackingLoop(name);
  }
  return loop;
}

TUnpackingLoop::TUnpackingLoop(std::string name)
  : StoppableThread(name),
    fOutputEvent(NULL), fRunStart(0),
    input_queue(std::make_shared<ThreadsafeQueue<std::vector<TRawEvent> > >()),
    output_queue(std::make_shared<ThreadsafeQueue<TUnpackedEvent*> >()) { }

TUnpackingLoop::~TUnpackingLoop() { }

bool TUnpackingLoop::Iteration(){
//  auto start = high_resolution_clock::now();

  std::vector<TRawEvent> event;
  int error = input_queue->Pop(event);
  if(error < 0){
    if(input_queue->IsFinished()) {
      output_queue->SetFinished();
      return false;
    } else {
      std::this_thread::sleep_for(std::chrono::milliseconds(1000));
      return true;
    }
  }
  //int counter=0;
  //printf("here 1\n"); fflush(stdout);
  fOutputEvent = new TUnpackedEvent;
  for(unsigned int i=0;i<event.size();i++) {
    TRawEvent& raw_event = event[i];
    //printf("counter %i\t%i\n",counter++, raw_event.GetFileType()); fflush(stdout);
    switch(raw_event.GetFileType()){
      case kFileType::NSCL_EVT:
      {
        TNSCLEvent nscl_event(raw_event);
        HandleNSCLData(nscl_event);
      }
        break;

      case kFileType::GRETINA_MODE2:
      case kFileType::GRETINA_MODE3:
      {
        TGEBEvent geb_event(raw_event);
        HandleGEBData(geb_event);
      }
        break;

      default:
        //printf("default\n"); fflush(stdout);
        break;
    }
  }

//  auto stop = high_resolution_clock::now();
//  auto duration = duration_cast<nanoseconds>(stop - start);
//  std::cout << event.size() << " Time taken in Build " << duration.count() << " microseconds" << std::endl;
//  auto start = high_resolution_clock::now();

  fOutputEvent->Build();
  fOutputEvent->SetRunStart(fRunStart);

  if(fOutputEvent->GetDetectors().size() != 0){
    output_queue->Push(fOutputEvent);
    fOutputEvent = NULL;
  }
//  auto stop2 = high_resolution_clock::now();
//  auto duration2 = duration_cast<nanoseconds>(stop2 - start);
//  std::cout << "Time taken in Detector " << duration2.count() << " microseconds" << std::endl;

  return true;
}

void TUnpackingLoop::ClearQueue() {
  std::vector<TRawEvent> raw_event;
  while(input_queue->Size()) {
    input_queue->Pop(raw_event);
  }

  while(output_queue->Size()){
    TUnpackedEvent* event = NULL;
    output_queue->Pop(event);
    if(event){

      delete event;
    }
  }
}

void TUnpackingLoop::HandleNSCLData(TNSCLEvent& event) {
//  auto start = high_resolution_clock::now();

  //printf("in handle nscl\t%i\n",event.GetEventType()); fflush(stdout);
  switch(event.GetEventType()) {
    case kNSCLEventType::BEGIN_RUN:            // 0x0001
    {
      TRawEvent::TNSCLBeginRun* begin = (TRawEvent::TNSCLBeginRun*)event.GetPayload();
      fRunStart = begin->unix_time;
    }
      break;
    case kNSCLEventType::END_RUN:              // 0x0002
    case kNSCLEventType::PAUSE_RUN:            // 0x0003
    case kNSCLEventType::RESUME_RUN:           // 0x0004
    case kNSCLEventType::PACKET_TYPES:         // 0x000A
    case kNSCLEventType::MONITORED_VARIABLES:  // 0x000B
    case kNSCLEventType::RING_FORMAT:          // 0x000C
    case kNSCLEventType::PHYSICS_EVENT_COUNT:  // 0x001F
    case kNSCLEventType::EVB_FRAGMENT:         // 0x0028
    case kNSCLEventType::EVB_UNKNOWN_PAYLOAD:  // 0x0029
    case kNSCLEventType::EVB_GLOM_INFO:        // 0x002A
    case kNSCLEventType::FIRST_USER_ITEM_CODE: // 0x8000
      break;
    case kNSCLEventType::PERIODIC_SCALERS:     // 0x0014
      //HandleNSCLPeriodicScalers(event);
      break;
    case kNSCLEventType::PHYSICS_EVENT:        // 0x001E
      HandleBuiltNSCLData(event);

/*      if(event.IsBuiltData()){
        HandleBuiltNSCLData(event);
      } else {
        HandleUnbuiltNSCLData(event);
      }
*/
      break;
  }
//  auto stop = high_resolution_clock::now();
//  auto duration = duration_cast<nanoseconds>(stop - start);
//  std::cout << "Time taken in NSCL " << duration.count() << " microseconds" << std::endl;
}


void TUnpackingLoop::HandleBuiltNSCLData(TNSCLEvent& event){
//  auto start = high_resolution_clock::now();
  TNSCLBuiltRingItem built(event);

  //printf("i am being called!!!\n"); fflush(stdout);
  //int counter=0;
  for(unsigned int i=0; i<built.NumFragments(); i++){
    //printf("\tcounter = %i\n",counter++); fflush(stdout);
    TNSCLFragment& fragment = built.GetFragment(i);

    int source_id = fragment.GetFragmentSourceID();
//  auto start = high_resolution_clock::now();
    kDetectorSystems detector = TDetectorEnv::Get().DetermineSystem(source_id);
//  auto stop = high_resolution_clock::now();
//  auto duration = duration_cast<nanoseconds>(stop - start);
//  std::cout << built.NumFragments() << " Time in Built " << duration.count() << " nanoseconds" << std::endl;


    TRawEvent frag_event = fragment.GetNSCLEvent();
    // S800 events in CAESAR/S800 event have no body header.
    // This grabs the timestamp from the fragment header so it can be used later.
/*    if(frag_event.GetTimestamp() == -1){
      frag_event.SetFragmentTimestamp(fragment.GetFragmentTimestamp());
    }
*/
//    fOutputEvent->AddRawData(frag_event, kDetectorSystems::SEGA);
    fOutputEvent->AddRawData(frag_event, detector);
  }
//  auto stop = high_resolution_clock::now();
//  auto duration = duration_cast<nanoseconds>(stop - start);
//  std::cout << built.NumFragments() << " Time in Built " << duration.count() << " nanoseconds" << std::endl;

}

void TUnpackingLoop::HandleUnbuiltNSCLData(TNSCLEvent& event){
  kDetectorSystems detector = TDetectorEnv::Get().DetermineSystem(event);
  fOutputEvent->AddRawData(event, detector);
}

void TUnpackingLoop::HandleGEBMode3(TGEBEvent& event, kDetectorSystems system){
  TGEBMode3Event built(event);
  for(unsigned int i=0; i<built.NumFragments(); i++){
    TGEBEvent& fragment = built.GetFragment(i);
    fOutputEvent->AddRawData(fragment, system);
  }
}

void TUnpackingLoop::HandleNSCLPeriodicScalers(TNSCLEvent& event){
  TUnpackedEvent* scaler_event = new TUnpackedEvent;
  scaler_event->AddRawData(event, kDetectorSystems::NSCLSCALERS);
  scaler_event->Build();
  scaler_event->SetRunStart(fRunStart);
  output_queue->Push(scaler_event);
}

void TUnpackingLoop::HandleS800Scaler(TGEBEvent& event){
  TUnpackedEvent* scaler_event = new TUnpackedEvent;
  scaler_event->AddRawData(event, kDetectorSystems::S800SCALER);
  scaler_event->Build();
  scaler_event->SetRunStart(fRunStart);
  output_queue->Push(scaler_event);
}

void TUnpackingLoop::HandleGEBData(TGEBEvent& event){
  int type = event.GetEventType();

  switch(type) {
    case 1: // Gretina decomp data.
      fOutputEvent->AddRawData(event, kDetectorSystems::GRETINA);
      break;
    case 2: // Gretina Mode3 data.
      HandleGEBMode3(event, kDetectorSystems::MODE3);
      break;
    case 5: // S800 Mode2 equvilant.
      fOutputEvent->AddRawData(event, kDetectorSystems::S800);
      break;
    case 6: // NSCL Non Event Data Typically Run Start/Stop and scalers
      break;
    case 8: // Bank 88 (Formerly Bank29)
      HandleGEBMode3(event, kDetectorSystems::BANK88);
      break;
    case 9: // Simulated S800 data from UCGretina
      fOutputEvent->AddRawData(event, kDetectorSystems::S800_SIM);
      break;
    case 10: // S800 scaler data
      HandleS800Scaler(event);
      break;
    case 11: // UCGretina simulated data.
      fOutputEvent->AddRawData(event,kDetectorSystems::GRETINA_SIM);
      break;
    case 17: // Phoswich Wall Mode2 equivlant.
      break;
    case 19: // GODDESS
      break;
    case 21: // Lenda event format identical to Type25
      fOutputEvent->AddRawData(event, kDetectorSystems::LENDA);
      break;
    case 22: // FastScint - S.G. Unsure of this data, type22 should be GODDESS auxillary from GEBHeaders.h
      fOutputEvent->AddRawData(event, kDetectorSystems::FASTSCINT);
      break;
    case 25: // General NSCL DDAS format identical to Type21
      fOutputEvent->AddRawData(event, kDetectorSystems::DDAS);
      break;
    case 29: // Something - S.G. Unsure of this type
      fOutputEvent->AddRawData(event, kDetectorSystems::BANK88);
      break;
    default:
      std::cout << "Unknown EventType: " << type << std::endl;
      break;
  }
}
