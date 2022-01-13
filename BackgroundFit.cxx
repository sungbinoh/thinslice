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

void bkgFit_mu(TFile *fmc, TFile *fdata, bool fitfakedata){
  const char varname[50] = "daughter_michel_score";
  const char particle[10] = "mu";
  
  double totaldata = 0;
  double totalmc = 0;
  TH1D *hvarSlice[pi::nthinslices][pi::nCuts][pi::nIntTypes+1];
  TH1D *hvar[pi::nCuts][pi::nIntTypes+1];

  for (int i = 0; i < pi::nCuts; ++i){
    for (int j = 0; j < pi::nIntTypes+1; ++j){
      for (int k = 0; k < pi::nthinslices; ++k){
        if (j==0){ // data
          hvarSlice[k][i][j] = (TH1D*)fdata->Get(Form("h%sSlice_%d_%d_%d",varname,k,i,j));
          if (i==6) totaldata += hvarSlice[k][i][j]->Integral(); // normalization should be done before any cut or after all cuts?
        }
        else{ // MC
          hvarSlice[k][i][j] = (TH1D*)fmc->Get(Form("h%sSlice_%d_%d_%d",varname,k,i,j));
          if (i==6) totalmc += hvarSlice[k][i][j]->Integral();
        }
        if (k==0){
          if (j==0){ // data
            hvar[i][j] = (TH1D*)fdata->Get(Form("h%s_bkg_%d_%d",varname,i,j));
          }
          else{ // MC
            hvar[i][j] = (TH1D*)fmc->Get(Form("h%s_bkg_%d_%d",varname,i,j));
          }
        }
      }
    }
  }
  std::cout<<"#Data "<<totaldata<<"; #MC "<<totalmc<<std::endl;
  TemplateFitter fitter;

  std::vector<double> vslice;
  std::vector<double> vcorr;
  std::vector<double> vcorrerr;

  // muon constraint begins
  for (int i = 0; i<pi::nthinslices; ++i){
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
    fitter.SetFitRange(6, 10);
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
  // components to be rescaled
  TH1D *h2 = hvar[6][pi::kMuon];
  h2->Add(hvar[6][pi::kMIDmu]);
  
  h1->Scale(totaldata/totalmc);
  h2->Scale(totaldata/totalmc);
  fitter.SetHistograms(h0, h1, h2);
  fitter.SetFitRange(h0->FindBin(0.55), -1);
  fitter.Fit();
  std::cout<<fitter.GetPar()<<" "<<fitter.GetParError()<<std::endl;
  // muon constraint ends
  
  TCanvas *c1 = new TCanvas("c1","c1");
  TGraphErrors *gr_corr = new TGraphErrors(vslice.size(), &vslice[0], &vcorr[0], 0, &vcorrerr[0]);
  gr_corr->Draw("ALP");
  gr_corr->SetTitle(Form("%s bkg constraint", particle));
  gr_corr->GetXaxis()->SetTitle("Slice");
  gr_corr->GetYaxis()->SetTitle("Scale factor");
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
}

void bkgFit_p(TFile *fmc, TFile *fdata, bool fitfakedata){
  const char varname[50] = "Chi2_proton"; // "Chi2_proton"/"mediandEdx"
  const char particle[10] = "p";
  
  double totaldata = 0;
  double totalmc = 0;
  TH1D *hvarSlice[pi::nthinslices][pi::nCuts][pi::nIntTypes+1];
  TH1D *hvar[pi::nCuts][pi::nIntTypes+1];

  for (int i = 0; i < pi::nCuts; ++i){
    for (int j = 0; j < pi::nIntTypes+1; ++j){
      for (int k = 0; k < pi::nthinslices; ++k){
        if (j==0){ // data
          hvarSlice[k][i][j] = (TH1D*)fdata->Get(Form("h%sSlice_%d_%d_%d",varname,k,i,j));
          if (i==6) totaldata += hvarSlice[k][i][j]->Integral(); // normalization should be done before any cut or after all cuts?
        }
        else{ // MC
          hvarSlice[k][i][j] = (TH1D*)fmc->Get(Form("h%sSlice_%d_%d_%d",varname,k,i,j));
          if (i==6) totalmc += hvarSlice[k][i][j]->Integral();
        }
        if (k==0){
          if (j==0){ // data
            hvar[i][j] = (TH1D*)fdata->Get(Form("h%s_bkg_%d_%d",varname,i,j));
          }
          else{ // MC
            hvar[i][j] = (TH1D*)fmc->Get(Form("h%s_bkg_%d_%d",varname,i,j));
          }
        }
      }
    }
  }
  std::cout<<"#Data "<<totaldata<<"; #MC "<<totalmc<<std::endl;
  TemplateFitter fitter;

  std::vector<double> vslice;
  std::vector<double> vcorr;
  std::vector<double> vcorrerr;

  // proton constraint begins
  for (int i = 0; i<pi::nthinslices; ++i){
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
    fitter.SetFitRange(2, 8);
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
  fitter.SetFitRange(2, 8);
  fitter.Fit();
  std::cout<<fitter.GetPar()<<" "<<fitter.GetParError()<<std::endl;
  // proton constraint ends
  
  TCanvas *c1 = new TCanvas("c1","c1");
  TGraphErrors *gr_corr = new TGraphErrors(vslice.size(), &vslice[0], &vcorr[0], 0, &vcorrerr[0]);
  gr_corr->Draw("ALP");
  gr_corr->SetTitle(Form("%s bkg constraint", particle));
  gr_corr->GetXaxis()->SetTitle("Slice");
  gr_corr->GetYaxis()->SetTitle("Scale factor");
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
}

void bkgFit_spi(TFile *fmc, TFile *fdata, bool fitfakedata){
  const char varname[50] = "costheta";
  const char particle[10] = "spi";
  
  double totaldata = 0;
  double totalmc = 0;
  TH1D *hvarSlice[pi::nthinslices][pi::nCuts][pi::nIntTypes+1];
  TH1D *hvar[pi::nCuts][pi::nIntTypes+1];

  for (int i = 0; i < pi::nCuts; ++i){
    for (int j = 0; j < pi::nIntTypes+1; ++j){
      for (int k = 0; k < pi::nthinslices; ++k){
        if (j==0){ // data
          hvarSlice[k][i][j] = (TH1D*)fdata->Get(Form("h%sSlice_%d_%d_%d",varname,k,i,j));
          if (i==6) totaldata += hvarSlice[k][i][j]->Integral(); // normalization should be done before any cut or after all cuts?
        }
        else{ // MC
          hvarSlice[k][i][j] = (TH1D*)fmc->Get(Form("h%sSlice_%d_%d_%d",varname,k,i,j));
          if (i==6) totalmc += hvarSlice[k][i][j]->Integral();
        }
        if (k==0){
          if (j==0){ // data
            hvar[i][j] = (TH1D*)fdata->Get(Form("h%s_bkg_%d_%d",varname,i,j));
          }
          else{ // MC
            hvar[i][j] = (TH1D*)fmc->Get(Form("h%s_bkg_%d_%d",varname,i,j));
          }
        }
      }
    }
  }
  std::cout<<"#Data "<<totaldata<<"; #MC "<<totalmc<<std::endl;
  TemplateFitter fitter;

  std::vector<double> vslice;
  std::vector<double> vcorr;
  std::vector<double> vcorrerr;

  // secondary pion constraint begins
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
    fitter.SetFitRange(2, 6); // don't use underflow bin
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
  fitter.SetFitRange(2, h0->FindBin(0.955));
  fitter.Fit();
  std::cout<<fitter.GetPar()<<" "<<fitter.GetParError()<<std::endl;
  // secondary pion constraint ends
  
  TCanvas *c1 = new TCanvas("c1","c1");
  TGraphErrors *gr_corr = new TGraphErrors(vslice.size(), &vslice[0], &vcorr[0], 0, &vcorrerr[0]);
  gr_corr->Draw("ALP");
  gr_corr->SetTitle(Form("%s bkg constraint", particle));
  gr_corr->GetXaxis()->SetTitle("Slice");
  gr_corr->GetYaxis()->SetTitle("Scale factor");
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
  
  TFile *fmc = TFile::Open("/dune/app/users/yinrui/thinslice/build/mcprod4a.root");
  TFile *fdata;
  if (fitfakedata){
    fdata = TFile::Open("/dune/app/users/yinrui/thinslice/build/mcprod4a.root");
  }
  else{
    fdata = TFile::Open("/dune/app/users/yinrui/thinslice/build/data.root");
  }

  bkgFit_mu(fmc, fdata, fitfakedata);
  bkgFit_p(fmc, fdata, fitfakedata);
  bkgFit_spi(fmc, fdata, fitfakedata);

  return 0;
}
