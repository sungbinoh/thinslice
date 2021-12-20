#include "SliceParams.h"
#include "EventType.h"
#include "EventSelection.h"
#include "TemplateFitter.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TVectorD.h"
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

  TFile *fmc = TFile::Open("/dune/app/users/yinrui/thinslice/build/mcprod4a_nospicut.root");

  TFile *fdata;
  if (fitfakedata){
    fdata = TFile::Open("/dune/app/users/yinrui/thinslice/build/mcprod4a_nospicut.root");
  }
  else{
    fdata = TFile::Open("/dune/app/users/yinrui/thinslice/build/data_nospicut.root");
  }

  double totaldata = 0;
  double totalmc = 0;
  
  TH1D *hvarSlice[pi::nthinslices][pi::nCuts][pi::nIntTypes+1];
  TH1D *hvar[pi::nCuts][pi::nIntTypes+1];
  const char varname[100] = "costheta"; // "daughter_michel_score" "mediandEdx" "costheta"
  const char particle[100] = "secpion"; // "muon" "proton" "secpion"

  for (int i = 0; i < pi::nCuts; ++i){
    for (int j = 0; j < pi::nIntTypes+1; ++j){
      for (int k = 0; k < pi::nthinslices; ++k){
        if (j==0){ // data
          hvarSlice[k][i][j] = (TH1D*)fdata->Get(Form("h%sSlice_%d_%d_%d",varname,k,i,j));
          if (i==0) totaldata += hvarSlice[k][i][j]->Integral();
        }
        else{ // MC
          hvarSlice[k][i][j] = (TH1D*)fmc->Get(Form("h%sSlice_%d_%d_%d",varname,k,i,j));
          if (i==0) totalmc += hvarSlice[k][i][j]->Integral();
        }
        if (k==0){
          if (j==0){ // data
            hvar[i][j] = (TH1D*)fdata->Get(Form("h%s_%d_%d",varname,i,j));
          }
          else{ // MC
            hvar[i][j] = (TH1D*)fmc->Get(Form("h%s_%d_%d",varname,i,j));
          }
        }
      }
    }
  }
  std::cout<<totaldata<<" "<<totalmc<<std::endl;
  TemplateFitter fitter;

  std::vector<double> vslice;
  std::vector<double> vcorr;
  std::vector<double> vcorrerr;

  // constrain proton
  /*for (int i = 0; i<pi::nthinslices; ++i){
    std::cout<<"##### Slice "<<i<<std::endl;
    vslice.push_back(i);
    TH1D *h0 = hvarSlice[i][6][pi::kData];
    TH1D *h1 = hvarSlice[i][6][pi::kPiInel];
    h1->Add(hvarSlice[i][6][pi::kPiElas]);
    h1->Add(hvarSlice[i][6][pi::kMuon]);
    h1->Add(hvarSlice[i][6][pi::kMIDcosmic]);
    h1->Add(hvarSlice[i][6][pi::kMIDpi]);
    h1->Add(hvarSlice[i][6][pi::kMIDmu]);
    h1->Add(hvarSlice[i][6][pi::kMIDeg]);
    h1->Add(hvarSlice[i][6][pi::kMIDother]);
    // components to be rescaled
    TH1D *h2 = hvarSlice[i][6][pi::kMIDp];
    
    h1->Scale(totaldata/totalmc);
    h2->Scale(totaldata/totalmc);
    fitter.SetHistograms(h0, h1, h2);
    fitter.SetFitRange(4,-1);
    fitter.Fit();
    vcorr.push_back(fitter.GetPar());
    vcorrerr.push_back(fitter.GetParError());
  }
  cout<<"##### Global fit #####"<<endl;
  TH1D *h0 = hvar[6][pi::kData];
  TH1D *h1 = hvar[6][pi::kPiInel];
  h1->Add(hvar[6][pi::kPiElas]);
  h1->Add(hvar[6][pi::kMuon]);
  h1->Add(hvar[6][pi::kMIDcosmic]);
  h1->Add(hvar[6][pi::kMIDpi]);
  h1->Add(hvar[6][pi::kMIDmu]);
  h1->Add(hvar[6][pi::kMIDeg]);
  h1->Add(hvar[6][pi::kMIDother]);
  
  TH1D *h2 = hvar[6][pi::kMIDp];
  h1->Scale(totaldata/totalmc);
  h2->Scale(totaldata/totalmc);
  fitter.SetHistograms(h0, h1, h2);
  fitter.SetFitRange(h0->FindBin(2.5),-1);
  fitter.Fit();
  std::cout<<fitter.GetPar()<<" "<<fitter.GetParError()<<std::endl;*/
  
  // constrain muon
  /*for (int i = 0; i<pi::nthinslices; ++i){
    std::cout<<"##### Slice "<<i<<std::endl;
    vslice.push_back(i);
    TH1D *h0 = hvarSlice[i][6][pi::kData];
    TH1D *h1 = hvarSlice[i][6][pi::kPiInel];
    h1->Add(hvarSlice[i][6][pi::kPiElas]);
    h1->Add(hvarSlice[i][6][pi::kMIDp]);
    h1->Add(hvarSlice[i][6][pi::kMIDcosmic]);
    h1->Add(hvarSlice[i][6][pi::kMIDpi]);
    h1->Add(hvarSlice[i][6][pi::kMIDeg]);
    h1->Add(hvarSlice[i][6][pi::kMIDother]);
    // components to be rescaled
    TH1D *h2 = hvarSlice[i][6][pi::kMuon];
    h2->Add(hvarSlice[i][6][pi::kMIDmu]);
    
    h1->Scale(totaldata/totalmc);
    h2->Scale(totaldata/totalmc);
    fitter.SetHistograms(h0, h1, h2);
    fitter.SetFitRange(6,10);
    fitter.Fit();
    vcorr.push_back(fitter.GetPar());
    vcorrerr.push_back(fitter.GetParError());
  }
  cout<<"##### Global fit #####"<<endl;
  TH1D *h0 = hvar[6][pi::kData];
  TH1D *h1 = hvar[6][pi::kPiInel];
  h1->Add(hvar[6][pi::kPiElas]);
  h1->Add(hvar[6][pi::kMIDp]);
  h1->Add(hvar[6][pi::kMIDcosmic]);
  h1->Add(hvar[6][pi::kMIDpi]);
  h1->Add(hvar[6][pi::kMIDeg]);
  h1->Add(hvar[6][pi::kMIDother]);
  TH1D *h2 = hvar[6][pi::kMuon];
  h2->Add(hvar[6][pi::kMIDmu]);
  h1->Scale(totaldata/totalmc);
  h2->Scale(totaldata/totalmc);
  fitter.SetHistograms(h0, h1, h2);
  fitter.SetFitRange(h0->FindBin(0.55),-1);
  fitter.Fit();
  std::cout<<fitter.GetPar()<<" "<<fitter.GetParError()<<std::endl;*/

  // constrain secondary pion
  for (int i = 0; i<pi::nthinslices; ++i){
    std::cout<<"##### Slice "<<i<<std::endl;
    vslice.push_back(i);
    TH1D *h0 = hvarSlice[i][6][pi::kData];
    TH1D *h1 = hvarSlice[i][6][pi::kPiInel];
    h1->Add(hvarSlice[i][6][pi::kPiElas]);
    h1->Add(hvarSlice[i][6][pi::kMuon]);
    h1->Add(hvarSlice[i][6][pi::kMIDcosmic]);
    h1->Add(hvarSlice[i][6][pi::kMIDp]);
    h1->Add(hvarSlice[i][6][pi::kMIDmu]);
    h1->Add(hvarSlice[i][6][pi::kMIDeg]);
    h1->Add(hvarSlice[i][6][pi::kMIDother]);
    // components to be rescaled
    TH1D *h2 = hvarSlice[i][6][pi::kMIDpi];
    
    h1->Scale(totaldata/totalmc);
    h2->Scale(totaldata/totalmc);
    fitter.SetHistograms(h0, h1, h2);
    fitter.SetFitRange(1, h0->FindBin(0.96));
    fitter.Fit();
    vcorr.push_back(fitter.GetPar());
    vcorrerr.push_back(fitter.GetParError());
  }
  cout<<"##### Global fit #####"<<endl;
  TH1D *h0 = hvar[6][pi::kData];
  TH1D *h1 = hvar[6][pi::kPiInel];
  h1->Add(hvar[6][pi::kPiElas]);
  h1->Add(hvar[6][pi::kMuon]);
  h1->Add(hvar[6][pi::kMIDcosmic]);
  h1->Add(hvar[6][pi::kMIDp]);
  h1->Add(hvar[6][pi::kMIDmu]);
  h1->Add(hvar[6][pi::kMIDeg]);
  h1->Add(hvar[6][pi::kMIDother]);
  
  TH1D *h2 = hvar[6][pi::kMIDpi];
  h1->Scale(totaldata/totalmc);
  h2->Scale(totaldata/totalmc);
  fitter.SetHistograms(h0, h1, h2);
  fitter.SetFitRange(1, h0->FindBin(0.96));
  fitter.Fit();
  std::cout<<fitter.GetPar()<<" "<<fitter.GetParError()<<std::endl;

  TCanvas *c1 = new TCanvas("c1","c1");
  TGraphErrors *gr_corr = new TGraphErrors(vslice.size(), &vslice[0], &vcorr[0], 0, &vcorrerr[0]);
  gr_corr->Draw("ALP");
  if (fitfakedata){
    //TFile f(Form("bkgfitsMC_%s.root", particle), "recreate");
    //gr_corr->Write(Form("gr_corr_%s", particle));
    //f.Close();
    c1->Print(Form("bkgfitsMC_%s.png", particle));
  }
  else{
    //TFile f(Form("bkgfits_%s.root", particle), "recreate");
    //gr_corr->Write(Form("gr_corr_%s", particle));
    //f.Close();
    c1->Print(Form("bkgfits_%s.png", particle));
  }

  return 0;

}
