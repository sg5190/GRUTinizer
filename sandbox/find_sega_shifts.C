#include <map>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "TSega.h"
#include "TVector3.h"
#include "TMinuit.h"

using namespace std;
const double SBETA = 0.088338;
const int DetNum = 16;
TSega *sega = new TSega();

std::vector<float> senergies;

double invDoppler(double E, double theta){
  double gamma = 1./(sqrt(1.-pow(SBETA,2.)));
  return E / (gamma *(1 - SBETA*TMath::Cos(theta)));
}

void labFrame(vector<float> &eng){
  TSega *sega = new TSega();
  double x = 1.340096;
  double y = -0.221872;
  double z = 2.8;
  TVector3 beam = TVector3(x,y,z);
  int ncrystals = (int) eng.size();
  for (int c=0; c < ncrystals; c++){
    TVector3 crstl_pos = sega->GetSegmentPosition(DetNum, c + 1);
cout << eng[c] << endl;
    eng[c] = invDoppler(eng[c],crstl_pos.Angle(beam));
  }
  return;
}

double dopplerCorrect(int idx, double xshift, double yshift, double zshift){
  double x = 1.340096;
  double y = -0.221872;
  double z = 2.8;
  TVector3 track = TVector3(x,y,z);
  TVector3 sega_pos = sega->GetSegmentPosition(DetNum, idx);
  sega_pos.SetXYZ(sega_pos.x() + xshift, sega_pos.y() + yshift, sega_pos.z() + zshift);
  double gamma = 1./(sqrt(1.-pow(SBETA,2.)));
  return gamma*(1 - SBETA*TMath::Cos(sega_pos.Angle(track)));
}

void fStdev(int &npar, double *gin, double &f, double *par, int iflag){
  double enAvg = 0;
  int ncrystals = (int) senergies.size();
  vector<double> dop_en;
  //make list and find avg
  for (int c=0; c < ncrystals; c++){
    dop_en.push_back(senergies[c]*dopplerCorrect(c + 1, par[0], par[1], par[2]));
    enAvg += dop_en.back();
  }
  enAvg /= ncrystals;
  //calc stdev
  double stdev = 0;
  for (int c=0; c < ncrystals; c++){
    stdev += (enAvg - dop_en[c])*(enAvg - dop_en[c]);
  }
  stdev /= ncrystals;
  f = TMath::Sqrt(stdev);
}

void printFitEnergies(double *par){
  int ncrystals = (int) senergies.size();
  for (int c=0; c < ncrystals; c++){
    printf("%f\n%f\n",senergies[c],senergies[c]*dopplerCorrect(c + 1, par[0], par[1], par[2]));
  }
  return;
}

void find_optimal_shifts(vector<float> eng){
  if(eng.size()==0) {
    cout << "Energy vector empty" << endl;
    return 0;
  }
  senergies.insert(senergies.end(), eng.begin(), eng.end());
  labFrame(senergies);

  //initialize
  int npars = 3;
  TMinuit *min = new TMinuit(npars);
  min->DefineParameter(0,"x_shift",1,0.001,-10.0,10.0);
  min->DefineParameter(1,"y_shift",1,0.001,-10.0,10.0);
  min->DefineParameter(2,"z_shift",1,0.001,-10.0,10.0);
  min->SetFCN(fStdev);
  //do minimization
  double pars[1], parerrs[1];
  min->Migrad();
  for (int p=0; p < npars; p++){
    min->GetParameter(p,pars[p],parerrs[p]);
    cout << p << " " << pars[p] << " " << parerrs[p] << endl;
  }
  printFitEnergies(pars);
  return 1;
}
