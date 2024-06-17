// To make a new filter, copy this file under a new name in the "filter" directory.
// The "FilterCondition" function should return a boolean value.
// The boolean indicates whether the event should be kept or not.

#include "TRuntimeObjects.h"

#include <TGretina.h>
#include <TS800.h>
#include <TGates.h>

bool gatesloaded = false;
bool incoming = false;
bool outgoing = false;

TGates *in = new TGates();
TGates *out = new TGates();

extern "C"
bool FilterCondition(TRuntimeObjects& obj) {


  TS800 *s800 = obj.GetDetector<TS800>();
  if(!s800) return false;

  if(!gatesloaded) {
    incoming = in->LoadPIDFile(Form("%s", GValue::Info("INCOMINGPIDFILE").c_str()));
    outgoing = out->LoadPIDFile(Form("%s", GValue::Info("OUTGOINGPIDFILE").c_str()));
    gatesloaded = true;
  }

  float tmpOBJ = s800->GetMTof().GetCorrelatedTof("Obj", "Ref", GValue::Value("TARGET_MTOF_OBJE1"), GValue::Value("SHIFT_MTOF_OBJE1"));
  float tmpXFP = s800->GetMTof().GetCorrelatedTof("Xfp", "Ref", GValue::Value("TARGET_MTOF_XFPE1"), GValue::Value("SHIFT_MTOF_XFPE1"));
  if(in->GateID(tmpOBJ, tmpXFP) < 0) return false;

  int maxp0 = s800->GetCrdc(0).GetMaxPad();
  int maxp1 = s800->GetCrdc(1).GetMaxPad();
  float xfp0 = s800->GetCrdc(0).GetDispersiveX(maxp0);
  float yfp0 = s800->GetCrdc(0).GetNonDispersiveY(maxp0);
  float xfp1 = s800->GetCrdc(1).GetDispersiveX(maxp1);
  float afp = s800->GetAFP(xfp0, xfp1);

  float tmpOBJE1 = s800->GetMTofCorr(tmpOBJ, afp, xfp0, GValue::Value("OBJ_MTOF_CORR_AFP"),  GValue::Value("OBJ_MTOF_CORR_XFP"));
  float tmpIC = s800->GetIonChamber().GetdE(xfp0, yfp0);
  if(!outgoing) return true;
  if(out->GateID(tmpOBJE1,tmpIC) < 0) return false;

  return true;
}



