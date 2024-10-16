void loadMacros() {
  std::string grsipath = std::string(getenv("GRUTSYS"));
  gROOT->ProcessLine(Form(".L %s/sandbox/PIDHelper.C", grsipath.c_str()));
  gROOT->ProcessLine(Form(".L %s/sandbox/fit.C", grsipath.c_str()));
//  gROOT->ProcessLine(Form(".L %s/sandbox/find_optimal_shifts.C", grsipath.c_str()));
}
