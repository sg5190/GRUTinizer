#ifndef _TGRUTTYPES_H_
#define _TGRUTTYPES_H_

#include <map>
#include <string>

#include <Rtypes.h>

#include "TDetectorFactory.h"

enum kDetectorSystems {
  UNKNOWN_SYSTEM = -1,
  GRETINA     = 1,
  MODE3       = 2,
  FASTSCINT   = 4,
  S800        = 5,
  S800SCALER  = 10,
  BANK88      = 8,
  GRETINA_SIM = 11,
  S800_SIM    = 9,
  LENDA       = 21,
  DDAS	      = 25,
  SEGA        = 64,
  JANUS       = 65,
  JANUS_DDAS  = 66,
  SUN	      = 70,
  CAESAR      = 80,
  NSCLSCALERS = 100
};

extern std::map<std::string, kDetectorSystems> detector_system_map;

enum kFileType {
  UNKNOWN_FILETYPE = -1,
  NSCL_EVT = 1,
  GRETINA_MODE2 = 2,
  GRETINA_MODE3 = 3,
  ROOT_DATA = 256,
  ROOT_MACRO = 257,
  CALIBRATED = 512,
  GVALUE     = 513,
  PRESETWINDOW = 514,
  DETECTOR_ENVIRONMENT = 1024,
  GUI_HIST_FILE = 2048,
  COMPILED_HISTOGRAMS = 2049,
  COMPILED_FILTER = 2053,
  CONFIG_FILE = 2050,
  S800_INVMAP = 2051,
  CUTS_FILE = 2052
};

extern std::map<std::string, EColor> color_system_map;

extern std::map<kDetectorSystems, TDetectorFactoryBase*> detector_factory_map;

#endif /* _TGRUTTYPES_H_ */
