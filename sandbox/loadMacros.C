void co60Cal(float p1, float p2) {
  float gain = (1332.5 - 1173.2)/(p2 - p1);
  float off = 1173.2 - gain * p1 ;
  cout << gain << "\t" << off << endl;
}

void na22Cal(float p1, float p2) {
  float gain = (1274.537 - 511.0)/(p2 - p1);
  float off = 511.0 - gain * p1 ;
  cout << gain << "\t" << off << endl;
}

void y88Cal(float p1, float p2) {
  float gain = (1836.061 - 898.042)/(p2 - p1);
  float off = 898.042 - gain * p1 ;
  cout << gain << "\t" << off << endl;
}

void loadMacros() {
  std::string grsipath = std::string(getenv("GRUTSYS"));
  gROOT->ProcessLine(Form(".L %s/sandbox/PIDHelper.C", grsipath.c_str()));
  gROOT->ProcessLine(Form(".L %s/sandbox/fit.C", grsipath.c_str()));
//  gROOT->ProcessLine(Form(".L %s/sandbox/find_optimal_shifts.C", grsipath.c_str()));
}
