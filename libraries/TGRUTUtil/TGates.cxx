#include "TGates.h"

ClassImp(TGates)

TGates::TGates() {
}

TGates::~TGates() { }

bool TGates::LoadPIDFile(const char * infile) {
  std::ifstream ifile(infile);
  if(!ifile.is_open()) {
    std::cout << "PID File " << infile << " not found" << std::endl;
    return false;
  }
  ifile >> fNPid;
  for(int i = 0; i < fNPid; i++) {
    TCutG *tmpG = new TCutG();
    ifile >> fZ >> fA;
    TNucleus *tmpN = new TNucleus(fZ, fA - fZ);
    ifile >> fNP;
    for(int j = 0; j < fNP; j++) {
      ifile >> fX >> fY;
      tmpG->SetPoint(j, fX, fY);
    }
    tmpG->SetName(Form("%s",tmpN->GetSymbol()));
    GateList.push_back(std::make_pair(tmpN, tmpG));
  }
  return true;
}

void TGates::MakePIDFile(const char * filename) {
  std::vector<std::string> data;
  int NPID = 0;
  while(1) {
    TCutG *cut1 = (TCutG*)gPad->WaitPrimitive("CUTG","CutG");
    cut1->SetName("cut1");
    cut1->Print();
    int n = cut1->GetN();
    double *x = cut1->GetX();
    double *y = cut1->GetY();
    int Z, A;
    std::cout << "Z A?" << std::endl;
    std::cin >> Z >> A;
    data.push_back(Form("%i %i",Z,A));
    data.push_back(Form("%i",n));
    for (int i = 0; i < n; i++) {
      data.push_back(Form("%f %f",*x++,*y++));
    }
    delete cut1;
    std::string con;
    std::cout << "Continue (y/n)" << std::endl;
    std::cin >> con;
    NPID++;
    if(con != "y") break;
  }


  std::ofstream outfile;

  std::string fn = filename;
  if(fn.substr(fn.find_last_of(".") + 1) == "pid") {
    outfile.open(Form("%s",filename));
  } else {
    outfile.open(Form("%s.pid",filename));
  }

  outfile << NPID << std::endl;
  for(unsigned int i = 0; i < data.size();i++) {
    outfile << data.at(i) << std::endl;
  }
  outfile.close();
  return;
}

std::string TGates::GetSymbol(int ID) {
  if(ID > (int)GateList.size() - 1) return "NA";
  else return GateList.at(ID).first->GetSymbol();

}

int TGates::GateID(float x, float y) {
  int i = 0;
  for(auto &iter: GateList) {
    if(iter.second->IsInside(x,y)) return i;
    i++;
  }
  return -1;
}

int TGates::GetA(int ID) {
  if(ID > (int)GateList.size() - 1) return -1;
  else return GateList.at(ID).first->GetA();

}

int TGates::GetZ(int ID) {
  if(ID > (int)GateList.size() - 1) return -1;
  else return GateList.at(ID).first->GetZ();

}

int TGates::GetN(int ID) {
  if(ID > (int)GateList.size() - 1) return -1;
  else return GateList.at(ID).first->GetN();

}

double TGates::GetMass(int ID) {
  if(ID > (int)GateList.size() - 1) return -1;
  else return GateList.at(ID).first->GetMass();
}

bool TGates::LoadGateFile(const char * infile) {
  std::ifstream ifile(infile);
  if(!ifile.is_open()) {
    std::cout << "Gate File " << infile << " not found" << std::endl;
    return false;
  }

  gate2D = new TCutG();
  ifile >> fNP;
  for(int j = 0; j < fNP; j++) {
    ifile >> fX >> fY;
    gate2D->SetPoint(j, fX, fY);
  }
  gate2D->SetName("gate2d");

  return true;
}

void TGates::MakeGateFile(const char * filename) {
  std::vector<std::string> data;
  while(1) {
    TCutG *cut1 = (TCutG*)gPad->WaitPrimitive("CUTG","CutG");
    cut1->SetName("cut1");
    cut1->Print();
    int n = cut1->GetN();
    double *x = cut1->GetX();
    double *y = cut1->GetY();
    data.push_back(Form("%i",n));
    for (int i = 0; i < n; i++) {
      data.push_back(Form("%f %f",*x++,*y++));
    }
    delete cut1;
    std::string con;
    std::cout << "Redraw (y/n)" << std::endl;
    std::cin >> con;
    if(con != "y") break;
    else data.clear();
  }
  std::ofstream outfile;

  std::string fn = filename;
  if(fn.substr(fn.find_last_of(".") + 1) == "gate") {
    outfile.open(Form("%s",filename));
  } else {
    outfile.open(Form("%s.gate",filename));
  }
  for(unsigned int i = 0; i < data.size();i++) {
    outfile << data.at(i) << std::endl;
  }
  outfile.close();
  return;
}

bool TGates::InGate(float X, float Y) {
  if(gate2D->IsInside(X, Y)) return true;
  else return false;
}
