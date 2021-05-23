#include "util.h"

void FillHistVec1D(TH1D *hist[nParTypes+1], const double &value, const int &partype){
  hist[0]->Fill(value);
  if (partype>=1 && partype < nParTypes+1){
    if (value<hist[partype]->GetXaxis()->GetXmax()){
      hist[partype]->Fill(value);
    }
    else{
      hist[partype]->Fill(hist[partype]->GetXaxis()->GetXmax()-0.000001);
    }
  }
}

void FillHistVec2D(TH2D *hist[nParTypes+1], const double &value1, const double &value2, const int &partype){
  hist[0]->Fill(value1, value2);
  if (partype>=1 && partype < nParTypes+1){
    hist[partype]->Fill(value1, value2);
  }
}
