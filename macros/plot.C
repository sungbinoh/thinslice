#include "../EventType.h"
#include "../EventSelection.h"
#include "TLegend.h"

int nc = 0;

TFile *f;

void PrintEvents(string name){
  
  TH1D *h[nCuts][nParTypes+1];
  for (int i = 0; i < nCuts; ++i){
    for (int j = 0; j < nParTypes+1; ++j){
      h[i][j] = (TH1D*)f->Get(Form("%s_%d_%d",name.c_str(),i,j));
    }
  }

  for (int i = 0; i<nCuts; ++i){
    cout<<cutName[i]<<endl;
    for (int j =0; j<5; ++j){
      cout<<parTypeName[j]<<" "<<h[i][j]->GetEntries()<<endl;
    }
  }
}

void plot1d(string name, int cut, string xtitle, string ytitle){

  TH1D *h[nCuts][nParTypes+1];
  for (int i = 0; i < nCuts; ++i){
    for (int j = 0; j < nParTypes+1; ++j){
      h[i][j] = (TH1D*)f->Get(Form("%s_%d_%d",name.c_str(),i,j));
    }
  }
  TCanvas *can = new TCanvas(Form("can_%d",nc), Form("can_%d",nc));
  THStack *hs = new THStack(Form("hs_%d",nc),"");
  h[cut][3]->SetFillColor(kRed);
  h[cut][3]->SetLineWidth(0);
  hs->Add(h[cut][3]);
  h[cut][2]->SetFillColor(kBlue);
  h[cut][2]->SetLineWidth(0);
  hs->Add(h[cut][2]);
  h[cut][1]->SetFillColor(kGreen);
  h[cut][1]->SetLineWidth(0);
  hs->Add(h[cut][1]);
  h[cut][4]->SetFillColor(6);
  h[cut][4]->SetLineWidth(0);
  hs->Add(h[cut][4]);
  hs->Draw("hist");
  hs->SetTitle(cutName[cut]);
  hs->GetXaxis()->SetTitle(xtitle.c_str());
  hs->GetYaxis()->SetTitle(ytitle.c_str());
  TLegend *leg = new TLegend(0.6,0.6,0.9,0.9);
  leg->SetFillStyle(0);
  leg->AddEntry(h[cut][3],"#pi^{+} inelastic","f");
  leg->AddEntry(h[cut][4],"#pi^{+} elastic","f");
  leg->AddEntry(h[cut][2],"#mu^{+}","f");
  leg->AddEntry(h[cut][1],"Misidentified","f");
  leg->Draw();
  can->Print(Form("plots/can_%s_%s.png",name.c_str(), cutName[cut]));
  can->Print(Form("plots/can_%s_%s.pdf",name.c_str(), cutName[cut]));
  gPad->SetLogy();
  can->Print(Form("plots/canlog_%s_%s.png",name.c_str(), cutName[cut]));
  can->Print(Form("plots/canlog_%s_%s.pdf",name.c_str(), cutName[cut]));
  ++nc;
}

void plot(){

  gStyle->SetOptStat(0);

  f = TFile::Open("../install/bin/mcprod4a.root");

  for (int i = 0; i<nCuts; ++i){
    plot1d("hmediandEdx", i, "Median dE/dx (MeV/cm)", "Events");
    plot1d("hdaughter_michel_score", i, "Daughter Michel Score", "Events");
    plot1d("hreco_beam_endZ_SCE", i, "Reco track end Z (cm)", "Events");
    plot1d("htrue_beam_endZ_SCE", i, "True pion/muon end Z (cm)", "Events");
    plot1d("htrue_sliceID", i, "True Slice ID", "Events");
    plot1d("hreco_sliceID", i, "Reco Slice ID", "Events");
    plot1d("hreco_true_beam_endZ_SCE", i, "Reco - True End Z (cm)", "Events");
  }

  PrintEvents("hdaughter_michel_score");

}
