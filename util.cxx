#include "util.h"

void FillHistVec1D(TH1D *hist[nIntTypes+1], const double &value, const int &partype){
  //hist[0]->Fill(value);
  if (partype>=0 && partype < nIntTypes+1){
    if (value<hist[partype]->GetXaxis()->GetXmin()){
      hist[partype]->Fill(hist[partype]->GetXaxis()->GetXmin());
    }
    else if (value<hist[partype]->GetXaxis()->GetXmax()){
      hist[partype]->Fill(value);
    }
    else{
      hist[partype]->Fill(hist[partype]->GetXaxis()->GetXmax()-0.000001);
    }
  }
}

void FillHistVec2D(TH2D *hist[nIntTypes+1], const double &value1, const double &value2, const int &partype){
  //hist[0]->Fill(value1, value2);
  if (partype>=0 && partype < nIntTypes+1){
    hist[partype]->Fill(value1, value2);
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

