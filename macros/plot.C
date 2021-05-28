#include "../EventType.h"
#include "../EventSelection.h"
#include "TLegend.h"

int nc = 0;
double totaldata = 0;
double totalmc = 0;

int colors[nParTypes] = {2, 3, 5, 7, 33, 9, 46, 28, 41};

TFile *fmc;
TFile *fdata;

void PrintEvents(string name){

  TH1D *hdata[nCuts];
  TH1D *hmc[nCuts][nParTypes+1];
  for (int i = 0; i < nCuts; ++i){
    for (int j = 0; j < nParTypes+1; ++j){
      hmc[i][j] = (TH1D*)fmc->Get(Form("%s_%d_%d",name.c_str(),i,j));
    }
    hdata[i] = (TH1D*)fdata->Get(Form("%s_%d_%d",name.c_str(),i,0));
  }
  
  for (int i = 0; i<nCuts; ++i){
    cout<<"=========="<<cutName[i]<<"=========="<<endl;
    double totalmc = 0;
    for (int j =0; j<nParTypes+1; ++j){
      if (j==0){
        cout<<parTypeName[j]<<" "<<hdata[i]->Integral()<<endl;
      }
      else{
        cout<<parTypeName[j]<<" "<<hmc[i][j]->Integral()<<endl;
        totalmc += hmc[i][j]->Integral();
      }
    }
    cout<<"Total MC: "<<totalmc<<endl;
  }
}

void plot1d(string name, int cut, string xtitle, string ytitle){

  static bool first = true;

  TH1D *hdata[nCuts];
  TH1D *hmc[nCuts][nParTypes+1];
  for (int i = 0; i < nCuts; ++i){
    for (int j = 0; j < nParTypes+1; ++j){
      hmc[i][j] = (TH1D*)fmc->Get(Form("%s_%d_%d",name.c_str(),i,j));
      if (i==0 && j!=0 && first){
        totalmc += hmc[i][j]->Integral();
      }
    }
    hdata[i] = (TH1D*)fdata->Get(Form("%s_%d_%d",name.c_str(),i,0));
    if (i==0 && first) totaldata = hdata[i]->Integral();
  }
  //std::cout<<totaldata<<" "<<totalmc<<std::endl;
  first = false;
  TCanvas *can = new TCanvas(Form("can_%d",nc), Form("can_%d",nc));
  THStack *hs = new THStack(Form("hs_%d",nc),"");
  for (int i = 0; i<nParTypes; ++i){
    hmc[cut][i+1]->SetFillColor(colors[i]);
    hmc[cut][i+1]->SetLineWidth(0);
    hmc[cut][i+1]->Scale(totaldata/totalmc);
    hs->Add(hmc[cut][i+1]);
  }
//  h[cut][3]->SetFillColor(kRed);
//  h[cut][3]->SetLineWidth(0);
//  hs->Add(h[cut][3]);
//  h[cut][2]->SetFillColor(kBlue);
//  h[cut][2]->SetLineWidth(0);
//  hs->Add(h[cut][2]);
//  h[cut][1]->SetFillColor(kGreen);
//  h[cut][1]->SetLineWidth(0);
//  hs->Add(h[cut][1]);
//  h[cut][4]->SetFillColor(6);
//  h[cut][4]->SetLineWidth(0);
//  hs->Add(h[cut][4]);
  double max = TMath::Max(hs->GetMaximum(), hdata[cut]->GetMaximum());
  hs->SetMaximum(1.1*max);
  hs->Draw("hist");
  hs->SetTitle(cutName[cut]);
  hs->GetXaxis()->SetTitle(xtitle.c_str());
  hs->GetYaxis()->SetTitle(ytitle.c_str());
  hdata[cut]->SetMarkerStyle(20);
  hdata[cut]->Draw("same");
  TLegend *leg;
  if (hmc[cut][0]->GetMaximumBin() < hmc[cut][0]->GetNbinsX()/2){
    leg = new TLegend(0.55,0.7,0.9,0.9);
  }
  else{
    leg = new TLegend(0.15,0.7,0.5,0.9);
  }
  leg->SetFillStyle(0);
  leg->SetNColumns(2);
  leg->AddEntry(hdata[cut],parTypeName[0],"ple");
  for (int i = 0; i<nParTypes; ++i){
    leg->AddEntry(hmc[cut][i+1], parTypeName[i+1],"f");
  }
//  leg->AddEntry(h[cut][3],"#pi^{+} inelastic","f");
//  leg->AddEntry(h[cut][4],"#pi^{+} elastic","f");
//  leg->AddEntry(h[cut][2],"#mu^{+}","f");
//  leg->AddEntry(h[cut][1],"Misidentified","f");
  leg->Draw();
  can->Print(Form("plots/can_%s_%d_%s.png", name.c_str(), cut, cutName[cut]));
  can->Print(Form("plots/can_%s_%d_%s.pdf", name.c_str(), cut, cutName[cut]));
  gPad->SetLogy();
  can->Print(Form("plots/canlog_%s_%d_%s.png", name.c_str(), cut, cutName[cut]));
  can->Print(Form("plots/canlog_%s_%d_%s.pdf", name.c_str(), cut, cutName[cut]));
  ++nc;
}

void plot2d(string name, int cut){

  TH2D *h[nCuts][nParTypes+1];
  for (int i = 0; i < nCuts; ++i){
    for (int j = 0; j < nParTypes+1; ++j){
      h[i][j] = (TH2D*)fmc->Get(Form("%s_%d_%d",name.c_str(),i,j));
    }
  }
  for (int i = 1; i<=4; ++i){
    TCanvas *can = new TCanvas(Form("can_%d",nc), Form("can_%d",nc));
    h[cut][i]->Draw("colz");
    can->Print(Form("plots/can_%s_%s_%s.png",name.c_str(), cutName[cut], parTypeName[i]));
    can->Print(Form("plots/can_%s_%s_%s.pdf",name.c_str(), cutName[cut], parTypeName[i]));
    ++nc;
  }
}
void plot(){

  gStyle->SetOptStat(0);
  gErrorIgnoreLevel = kWarning;

  fmc = TFile::Open("../install/bin/mcprod4a.root");
  fdata = TFile::Open("../install/bin/data.root");

  for (int i = 0; i<nCuts; ++i){
    plot1d("hmediandEdx", i, "Median dE/dx (MeV/cm)", "Events");
    plot1d("hdaughter_michel_score", i, "Daughter Michel Score", "Events");
    plot1d("hreco_beam_endZ_SCE", i, "Reco track end Z (cm)", "Events");
    plot1d("htrue_beam_endZ_SCE", i, "True pion/muon end Z (cm)", "Events");
    plot1d("htrue_sliceID", i, "True Slice ID", "Events");
    plot1d("hreco_sliceID", i, "Reco Slice ID", "Events");
    plot1d("hreco_true_beam_endZ_SCE", i, "Reco - True End Z SCE (cm)", "Events");
    plot2d("hreco_vs_true_beam_endZ_SCE",i);
    plot1d("hdeltax", i, "#Deltax/#sigma_{x}", "Events");
    plot1d("hdeltay", i, "#Deltay/#sigma_{y}", "Events");
    plot1d("hdeltaz", i, "#Deltaz/#sigma_{z}", "Events");
    plot1d("hcostheta", i, "cos#theta", "Events");
    plot1d("htrklen", i, "Track length (cm)", "Events");
    plot1d("hreco_beam_startX_SCE", i, "Reco track start X (cm)", "Events");
    plot1d("hreco_beam_startY_SCE", i, "Reco track start Y (cm)", "Events");
    plot1d("hreco_beam_startZ_SCE", i, "Reco track start Z (cm)", "Events");

    plot1d("hreco_beam_dcosX_SCE", i, "dcosX", "Events");
    plot1d("hreco_beam_dcosY_SCE", i, "dcosY", "Events");
    plot1d("hreco_beam_dcosZ_SCE", i, "dcosZ", "Events");

    plot1d("hreco_beam_angleX_SCE", i, "angleX", "Events");
    plot1d("hreco_beam_angleY_SCE", i, "angleY", "Events");
    plot1d("hreco_beam_angleZ_SCE", i, "angleZ", "Events");
  }

  PrintEvents("htrue_beam_endZ_SCE");

}
