#include <iostream>
#include <iomanip>
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include <numeric>
#include "GPeak.h"
using namespace std;

vector<float> fitCentroids() {
  float fitwidth = 30;
  vector<float> en;
  TH2D *tmpHist = NULL;
  TPad *pad = (TPad*)gPad;
  if(pad == NULL) {
    cout << "Need To Open 2D Histogram of Energy vs Detector Number" << endl;
    return en;
  }

  TList *newlist = pad->GetListOfPrimitives();
  for(auto obj: *newlist) {
    if((strcmp(obj->ClassName(), "TH2F") == 0 ) || (strcmp(obj->ClassName(), "TH2D") == 0) ) {
      tmpHist = (TH2D*) newlist->FindObject(obj->GetName());
      break;
    }
  }

  if(tmpHist == NULL) {
    cout << "Histogram Not Found" << endl;
    return en;
  }

  cout << "Select a Peak to Fit" << endl;
  float xp;
  while(1) {
    TMarker *mark = (TMarker*)gPad->WaitPrimitive("TMarker","TMarker");
    xp = mark->GetY();
    cout << "You Have Selected E\u03B3 = " << xp << " keV" << endl;
    cout << "Reselect (y/n)?" << endl;
    string con;
    cin >> con;
    if(con != "y") break;
  }

  while(1) {
    cout << "Select Fit Range " << xp << " +/- ? " << endl;
    cin >> fitwidth;
    cout << "Fit Range " << xp - fitwidth << " - " << xp + fitwidth << endl;
    cout << "Reselect (y/n)?" << endl;
    string con;
    cin >> con;
    if(con != "y") break;
  }

  cout << "Fitting ..." << endl;
  TGraphErrors *gr = new TGraphErrors();
  gr->SetTitle("Centroid vs DetNum");
  for(int i = 0; i < tmpHist->GetNbinsX(); i++) {
    char hname[192];
    sprintf(hname, "tmp1D%i",i);
    TH1D * tmp1D = new TH1D(hname,"",tmpHist->GetNbinsY(), tmpHist->GetYaxis()->GetXmin(),  tmpHist->GetYaxis()->GetXmax());
    tmpHist->ProjectionY(hname, i+1, i+1);
    if(tmp1D->Integral(0,-1) < 10) {
      //en.push_back(0);
      continue;
    }
    Int_t bin = tmpHist->GetXaxis()->FindBin(xp);
    Double_t yp = tmpHist->GetBinContent(bin);
    GGaus *tmpP = new GGaus(xp - fitwidth, xp + fitwidth, "Q && no-print");
    tmpP->Fit(tmp1D,"Q && no-print");
    en.push_back(tmpP->GetParameter(1));
    gr->SetPoint(i, i, tmpP->GetParameter(1));
    gr->SetPointError(i, 0, tmpP->GetParError(1));
  }
  gr->Draw("A*");
  return en;
}
