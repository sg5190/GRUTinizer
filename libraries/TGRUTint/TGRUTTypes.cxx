#include "TGRUTTypes.h"

#include "TBank88.h"
#include "TCaesar.h"
#include "TFastScint.h"
#include "TGretSim.h"
#include "TGretina.h"
#include "TJanus.h"
#include "TJanusDDAS.h"
#include "TLenda.h"
#include "TGenericDDAS.h"
#include "TMode3.h"
#include "TNSCLScalers.h"
#include "TS800Scaler.h"
#include "TS800Sim.h"
#include "TSega.h"
#include "TSun.h"
// Map from string to detector enum.
// This is used to parse the DetectorEnvironment.env file.
// This is ONLY used when parsing NSCL data files, not
std::map<std::string, kDetectorSystems> detector_system_map{
  {"Unknown",     kDetectorSystems::UNKNOWN_SYSTEM},
  {"Gretina",     kDetectorSystems::GRETINA},
  {"Mode3",       kDetectorSystems::MODE3},
  {"FastScint",   kDetectorSystems::FASTSCINT},
  {"S800",        kDetectorSystems::S800},
  {"S800_Scaler", kDetectorSystems::S800SCALER},
  {"Bank88",      kDetectorSystems::BANK88},
  {"Gretina_Sim", kDetectorSystems::GRETINA_SIM},
  {"S800_SIM",    kDetectorSystems::S800_SIM},
  {"Lenda",       kDetectorSystems::LENDA},
  {"DDAS",	  kDetectorSystems::DDAS},
  {"Sega",        kDetectorSystems::SEGA},
  {"Janus",       kDetectorSystems::JANUS},
  {"JanusDDAS",   kDetectorSystems::JANUS_DDAS},
  {"Sun",	  kDetectorSystems::SUN},
  {"Caesar",      kDetectorSystems::CAESAR},
  {"NSCL_Scalers",kDetectorSystems::NSCLSCALERS}
};

// Map from detector enum to detector factory.
// This is used to build each detector type as it is needed.
// Adding a detector to this list allows grutinizer to unpack that detector type.
std::map<kDetectorSystems, TDetectorFactoryBase*> detector_factory_map {
  {kDetectorSystems::GRETINA,     new TDetectorFactory<TGretina>() },
  {kDetectorSystems::MODE3,       new TDetectorFactory<TMode3>() },
  {kDetectorSystems::FASTSCINT,   new TDetectorFactory<TFastScint>() },
  {kDetectorSystems::S800,        new TDetectorFactory<TS800>() },
  {kDetectorSystems::S800SCALER,  new TDetectorFactory<TS800Scaler>() },
  {kDetectorSystems::BANK88,      new TDetectorFactory<TBank88>() },
  {kDetectorSystems::GRETINA_SIM, new TDetectorFactory<TGretSim>() },
  {kDetectorSystems::S800_SIM,    new TDetectorFactory<TS800Sim>() },
  {kDetectorSystems::LENDA,       new TDetectorFactory<TLenda>() },
  {kDetectorSystems::DDAS,        new TDetectorFactory<TGenericDDAS>() },
  {kDetectorSystems::SEGA,        new TDetectorFactory<TSega>() },
  {kDetectorSystems::JANUS,       new TDetectorFactory<TJanus>() },
  {kDetectorSystems::JANUS_DDAS,  new TDetectorFactory<TJanusDDAS>() },
  {kDetectorSystems::SUN,         new TDetectorFactory<TSun>() },
  {kDetectorSystems::CAESAR,      new TDetectorFactory<TCaesar>() },
  {kDetectorSystems::NSCLSCALERS, new TDetectorFactory<TNSCLScalers>() }
};

std::map<std::string, EColor> color_system_map {
  {"kWhite",   EColor::kWhite},
  {"kBlack",   EColor::kBlack},
  {"kGray",    EColor::kGray},
  {"kRed",     EColor::kRed},
  {"kGreen",   EColor::kGreen},
  {"kBlue",    EColor::kBlue},
  {"kYellow",  EColor::kYellow},
  {"kMagenta", EColor::kMagenta},
  {"kCyan",    EColor::kCyan},
  {"kOrange",  EColor::kOrange},
  {"kSpring",  EColor::kSpring},
  {"kTeal",    EColor::kTeal},
  {"kAzure",   EColor::kAzure},
  {"kViolet",  EColor::kViolet},
  {"kPink",    EColor::kPink}
};
