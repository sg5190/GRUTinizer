#include <map>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "TGretina.h"
#include "TVector3.h"
#include "TMinuit.h"

using namespace std;
const double BETA = 0.38541500;
TGretina *gret = new TGretina();

std::vector<float> energies;

double invDoppler(double E, double theta){
  double gamma = 1./(sqrt(1.-pow(BETA,2.)));
  return E / (gamma *(1 - BETA*TMath::Cos(theta)));
}

void labFrame(vector<float> &eng){
  TGretina *gret = new TGretina();
  double x, y = 0.0;
  TVector3 beam = TVector3(x,-y,TMath::Sqrt(1 - x*x - y*y));
  int ncrystals = (int) eng.size();
  for (int c=0; c < ncrystals; c++){
    TVector3 crstl_pos = gret->GetCrystalPosition(gret->GetSpecId(c,1));
    eng[c] = invDoppler(eng[c],crstl_pos.Angle(beam));
  }
  return;
}

double dopplerCorrect(int idx, double ataShift, double btaShift, double xshift, double yshift, double zshift){
  TVector3 track = TVector3(ataShift,-btaShift,sqrt(1-ataShift*ataShift-btaShift*btaShift));
  TVector3 gret_pos = gret->GetCrystalPosition(idx);
  gret_pos.SetXYZ(gret_pos.x() - xshift, gret_pos.y() - yshift, gret_pos.z() - zshift);
  double gamma = 1./(sqrt(1.-pow(BETA,2.)));
  return gamma*(1 - BETA*TMath::Cos(gret_pos.Angle(track)));
}

void fChi2(int &npar, double *gin, double &f, double *par, int iflag){
  double enAvg = 0;
  int ncrystals = (int) energies.size();
  vector<double> dop_en;
  //make list and find avg
  for (int c=0; c < ncrystals; c++){
    dop_en.push_back(energies[c]*dopplerCorrect(gret->GetSpecId(c,1),par[0],par[1],par[2],par[3],par[4]));
    enAvg += dop_en.back();
  }
  enAvg /= ncrystals;

  //calc chi2
  double chi2 = 0;
  for (int c=0; c < ncrystals; c++){
    chi2 += (enAvg - dop_en[c])*(enAvg - dop_en[c])/11; //need to add actual sigmas
  }
  f = chi2;
}

void fStdev(int &npar, double *gin, double &f, double *par, int iflag){
  double enAvg = 0;
  int ncrystals = (int) energies.size();
  vector<double> dop_en;
  //make list and find avg
  for (int c=0; c < ncrystals; c++){
    dop_en.push_back(energies[c]*dopplerCorrect(gret->GetSpecId(c,1),par[0],par[1],par[2],par[3],par[4]));
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
  int ncrystals = (int) energies.size();
  for (int c=0; c < ncrystals; c++){
    printf("%f\n",energies[c]*dopplerCorrect(gret->GetSpecId(c,1),par[0],par[1],par[2],par[3],par[4]));
  }
  return;
}

void find_optimal_shifts(vector<float> eng){
  if(eng.size()==0) {
    cout << "Energy vector empty" << endl;
    return 0;
  }
  energies.insert(energies.end(), eng.begin(), eng.end());
  labFrame(energies);

  //initialize
  int npars = 5;
  TMinuit *min = new TMinuit(npars);
  min->DefineParameter(0,"ata_shift",0.01,1E-3,-0.5,0.5);
  min->DefineParameter(1,"bta_shift",0.01,1E-3,-0.5,0.5);
  min->DefineParameter(2,"x_shift",1,0.01,-10.0,10.0);
  min->DefineParameter(3,"y_shift",1,0.01,-10.0,10.0);
  min->DefineParameter(4,"z_shift",1,0.01,-100.0,100.0);
  min->SetFCN(fStdev);
  //do minimization
  double pars[5], parerrs[5];
  min->Migrad();
  for (int p=0; p < npars; p++){
    min->GetParameter(p,pars[p],parerrs[p]);
    cout << pars[p] << endl;
  }
//  printFitEnergies(pars);
  return 1;
}
