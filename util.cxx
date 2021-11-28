#include "util.h"

void FillHistVec1D(TH1D *hist[pi::nIntTypes+1], const double &value, const int &partype){
  //hist[0]->Fill(value);
  double weight = 1.;
  if (partype == 3) // kMuon
    weight = 1.;//59;
  if (partype>=0 && partype < pi::nIntTypes+1){
    if (value<hist[partype]->GetXaxis()->GetXmin()){
      hist[partype]->Fill(hist[partype]->GetXaxis()->GetXmin(), weight);
    }
    else if (value<hist[partype]->GetXaxis()->GetXmax()){
      hist[partype]->Fill(value, weight);
    }
    else{
      hist[partype]->Fill(hist[partype]->GetXaxis()->GetXmax()-0.000001, weight);
    }
  }
}

void FillHistVec2D(TH2D *hist[pi::nIntTypes+1], const double &value1, const double &value2, const int &partype){
  double weight = 1.;
  if (partype == 3) // kMuon
    weight = 1.;//59;
  //hist[0]->Fill(value1, value2);
  if (partype>=0 && partype < pi::nIntTypes+1){
    hist[partype]->Fill(value1, value2, weight);
  }
}

void FillHist1D(TH1D *hist, const double &value, const double &wei){
  if (value<hist->GetXaxis()->GetXmin()){
    hist->Fill(hist->GetXaxis()->GetXmin(), wei);
  }
  else if (value<hist->GetXaxis()->GetXmax()){
    hist->Fill(value, wei);
  }
  else{
    hist->Fill(hist->GetXaxis()->GetXmax()-0.000001, wei);
  }
}

