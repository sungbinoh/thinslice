#include "SliceParams.h"
#include "EventType.h"
#include "EventSelection.h"
#include "TemplateFitter.h"
#include "TFile.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include <iostream>
#include <string>

static void show_usage(std::string name)
{
  std::cerr << "Usage: " << name <<" <option(s)>\n"
            << "Options:\n"
            << "\t-h,--help\t\tShow this help message\n"
            << "\t-f,--fakedata\t\tDo fit on fake data instead of data"
            << std::endl;
}


int main(int argc, char* argv[]){

  bool fitfakedata = false;

  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if ((arg == "-h") || (arg == "--help")) {
      show_usage(argv[0]);
      return 0;
    } else if ((arg == "-f") || (arg == "--fakedata")) {
      fitfakedata = true;
    }
  }

  TFile *fmc = TFile::Open("mcprod4a.root");

  TFile *fdata;
  if (fitfakedata){
    fdata = TFile::Open("mcprod4a.root");
  }
  else{
    fdata = TFile::Open("data.root");
  }

  TH1D *hmediandEdxSlice[nthinslices][nCuts][nIntTypes+1];
  TH1D *hdaughter_michel_scoreSlice[nthinslices][nCuts][nIntTypes+1];

  double totaldata = 0;
  double totalmc = 0;

  for (int i = 0; i < nCuts; ++i){
    for (int j = 0; j < nIntTypes+1; ++j){
      for (int k = 0; k < nthinslices; ++k){
        if (j==0){
          hmediandEdxSlice[k][i][j] = (TH1D*)fdata->Get(Form("hmediandEdxSlice_%d_%d_%d",k,i,j));
          hdaughter_michel_scoreSlice[k][i][j] = (TH1D*)fdata->Get(Form("hdaughter_michel_scoreSlice_%d_%d_%d",k,i,j));
          if (i==0) totaldata += hmediandEdxSlice[k][i][j]->Integral();
        }
        else{
          hmediandEdxSlice[k][i][j] = (TH1D*)fmc->Get(Form("hmediandEdxSlice_%d_%d_%d",k,i,j));
          hdaughter_michel_scoreSlice[k][i][j] = (TH1D*)fmc->Get(Form("hdaughter_michel_scoreSlice_%d_%d_%d",k,i,j));
          if (i==0) totalmc += hmediandEdxSlice[k][i][j]->Integral();
        }
      }
    }
  }
  std::cout<<totaldata<<" "<<totalmc<<std::endl;
  TemplateFitter fitter;

  std::vector<double> vslice;
  std::vector<double> vcorrproton;
  std::vector<double> vcorrprotonerr;

  for (int i = 0; i<nthinslices; ++i){
    std::cout<<"Slice "<<i<<std::endl;
    vslice.push_back(i);
    TH1D *h0 = hmediandEdxSlice[i][kAPA3][kData];
    TH1D *h1 = hmediandEdxSlice[i][kAPA3][kPiInel];
    h1->Add(hmediandEdxSlice[i][kAPA3][kPiElas]);
    h1->Add(hmediandEdxSlice[i][kAPA3][kMuon]);
    h1->Add(hmediandEdxSlice[i][kAPA3][kMIDcosmic]);
    h1->Add(hmediandEdxSlice[i][kAPA3][kMIDpi]);
    h1->Add(hmediandEdxSlice[i][kAPA3][kMIDmu]);
    h1->Add(hmediandEdxSlice[i][kAPA3][kMIDeg]);
    h1->Add(hmediandEdxSlice[i][kAPA3][kMIDother]);
    TH1D *h2 = hmediandEdxSlice[i][kAPA3][kMIDp];
    h1->Scale(totaldata/totalmc);
    h2->Scale(totaldata/totalmc);
    fitter.SetHistograms(h0, h1, h2);
    fitter.SetFitRange(5,8);
    fitter.Fit();
    vcorrproton.push_back(fitter.GetPar());
    vcorrprotonerr.push_back(fitter.GetParError());
  }

  TGraphErrors *gr_corr_proton = new TGraphErrors(vslice.size(), &vslice[0], &vcorrproton[0], 0, &vcorrprotonerr[0]);

  std::vector<double> vcorrmuon;
  std::vector<double> vcorrmuonerr;
  
  for (int i = 0; i<nthinslices; ++i){
    std::cout<<"Slice "<<i<<std::endl;
    TH1D *h0 = hdaughter_michel_scoreSlice[i][kAPA3][kData];
    TH1D *h1 = hdaughter_michel_scoreSlice[i][kAPA3][kPiInel];
    h1->Add(hdaughter_michel_scoreSlice[i][kAPA3][kPiElas]);
    h1->Add(hdaughter_michel_scoreSlice[i][kAPA3][kMIDp]);
    h1->Add(hdaughter_michel_scoreSlice[i][kAPA3][kMIDcosmic]);
    h1->Add(hdaughter_michel_scoreSlice[i][kAPA3][kMIDpi]);
    h1->Add(hdaughter_michel_scoreSlice[i][kAPA3][kMIDeg]);
    h1->Add(hdaughter_michel_scoreSlice[i][kAPA3][kMIDother]);
    TH1D *h2 = hdaughter_michel_scoreSlice[i][kAPA3][kMuon];
    h2->Add(hdaughter_michel_scoreSlice[i][kAPA3][kMIDmu]);
    h1->Scale(totaldata/totalmc);
    h2->Scale(totaldata/totalmc);
    fitter.SetHistograms(h0, h1, h2);
    fitter.SetFitRange(5,10);
    fitter.Fit();
    vcorrmuon.push_back(fitter.GetPar());
    vcorrmuonerr.push_back(fitter.GetParError());
  }

  TGraphErrors *gr_corr_muon = new TGraphErrors(vslice.size(), &vslice[0], &vcorrmuon[0], 0, &vcorrmuonerr[0]);


  if (fitfakedata){
    TFile f("backgroundfits_mc.root","recreate");
    gr_corr_proton->Write("gr_corr_proton");
    gr_corr_muon->Write("gr_corr_muon");
    f.Close();
  }
  else{
    TFile f("backgroundfits.root","recreate");
    gr_corr_proton->Write("gr_corr_proton");
    gr_corr_muon->Write("gr_corr_muon");
    f.Close();
  }

  return 0;

}
