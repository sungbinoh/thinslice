#include "../EventType.h"
#include "../EventSelection.h"
#include "../SliceParams.h"

void plotTemplates(){

  gStyle->SetOptStat(0);

  TFile *fmc = TFile::Open("../install/bin/mcprod4a.root");
  
  TH1D *hmediandEdxSlice[nthinslices][nCuts][nIntTypes+1];
  for (int i = 0; i < nCuts; ++i){
    for (int j = 0; j < nIntTypes+1; ++j){
      for (int k = 0; k < nthinslices; ++k){
        hmediandEdxSlice[k][i][j] = (TH1D*)fmc->Get(Form("hmediandEdxSlice_%d_%d_%d",k,i,j));
      }
    }
  }
  
  TCanvas *c1 = new TCanvas("c1", "c1", 800, 800);
  TLegend *leg1 = new TLegend(0.7,0.5,0.9,0.9);
  leg1->SetFillStyle(0);
  for (int i = 0; i<nthinslices; ++i){
    hmediandEdxSlice[i][kAPA3][kMIDp]->Scale(1./hmediandEdxSlice[i][kAPA3][kMIDp]->Integral());
    if (i<9){
      hmediandEdxSlice[i][kAPA3][kMIDp]->SetLineColor(i+1);
      hmediandEdxSlice[i][kAPA3][kMIDp]->SetMarkerColor(i+1);
    }
    else{
      hmediandEdxSlice[i][kAPA3][kMIDp]->SetLineColor(i+2);
      hmediandEdxSlice[i][kAPA3][kMIDp]->SetMarkerColor(i+2);
    }
    if (i==0) {
      hmediandEdxSlice[i][kAPA3][kMIDp]->Draw();
      hmediandEdxSlice[i][kAPA3][kMIDp]->GetYaxis()->SetRangeUser(0,1.1);
    }
    else{
      hmediandEdxSlice[i][kAPA3][kMIDp]->Draw("same");
    }
    leg1->AddEntry(hmediandEdxSlice[i][kAPA3][kMIDp],Form("slice%d",i),"ple");
  }
  leg1->Draw();

  TH1D *hdaughter_michel_scoreSlice[nthinslices][nCuts][nIntTypes+1];
  for (int i = 0; i < nCuts; ++i){
    for (int j = 0; j < nIntTypes+1; ++j){
      for (int k = 0; k < nthinslices; ++k){
        hdaughter_michel_scoreSlice[k][i][j] = (TH1D*)fmc->Get(Form("hdaughter_michel_scoreSlice_%d_%d_%d",k,i,j));
      }
    }
  }
  
  TCanvas *c2 = new TCanvas("c2", "c2", 800, 800);
  TLegend *leg2 = new TLegend(0.7,0.5,0.9,0.9);
  leg2->SetFillStyle(0);
  for (int i = 0; i<nthinslices; ++i){
    hdaughter_michel_scoreSlice[i][kAPA3][kMuon]->Add(hdaughter_michel_scoreSlice[i][kAPA3][kMIDmu]);
    hdaughter_michel_scoreSlice[i][kAPA3][kMuon]->Scale(1./hdaughter_michel_scoreSlice[i][kAPA3][kMuon]->Integral());
    if (i<9){
      hdaughter_michel_scoreSlice[i][kAPA3][kMuon]->SetLineColor(i+1);
      hdaughter_michel_scoreSlice[i][kAPA3][kMuon]->SetMarkerColor(i+1);
    }
    else{
      hdaughter_michel_scoreSlice[i][kAPA3][kMuon]->SetLineColor(i+2);
      hdaughter_michel_scoreSlice[i][kAPA3][kMuon]->SetMarkerColor(i+2);
    }
    if (i==0) {
      hdaughter_michel_scoreSlice[i][kAPA3][kMuon]->Draw();
      hdaughter_michel_scoreSlice[i][kAPA3][kMuon]->GetYaxis()->SetRangeUser(0,1.1);
    }
    else{
      hdaughter_michel_scoreSlice[i][kAPA3][kMuon]->Draw("same");
    }
    leg2->AddEntry(hdaughter_michel_scoreSlice[i][kAPA3][kMuon],Form("slice%d",i),"ple");
  }
  leg2->Draw();

}
