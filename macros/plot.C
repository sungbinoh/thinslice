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
  THStack *hs = new THStack("hs","");
  h[cut][3]->SetFillColor(kRed);
  h[cut][3]->SetLineColor(kRed);
  hs->Add(h[cut][3]);
  h[cut][2]->SetFillColor(kBlue);
  h[cut][2]->SetLineColor(kBlue);
  hs->Add(h[cut][2]);
  h[cut][1]->SetFillColor(kGreen);
  h[cut][1]->SetLineColor(kGreen);
  hs->Add(h[cut][1]);
  h[cut][4]->SetFillColor(6);
  h[cut][4]->SetLineColor(6);
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
  can->Print(Form("can_%s_%s.png",name.c_str(), cutName[cut]));
  can->Print(Form("can_%s_%s.pdf",name.c_str(), cutName[cut]));
  gPad->SetLogy();
  can->Print(Form("canlog_%s_%s.png",name.c_str(), cutName[cut]));
  can->Print(Form("canlog_%s_%s.pdf",name.c_str(), cutName[cut]));
  ++nc;
}

void plot(){

  gStyle->SetOptStat(0);

  f = TFile::Open("../install/bin/hadana.root");

  for (int i = 0; i<nCuts; ++i){
    plot1d("hmediandEdx", i, "Median dE/dx (MeV/cm)", "Events");
    plot1d("hdaughter_michel_score", i, "Daughter Michel Score", "Events");
  }

  PrintEvents("hdaughter_michel_score");

}
