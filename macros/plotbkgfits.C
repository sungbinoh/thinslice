#include "../EventType.h"
#include "../EventSelection.h"
#include "../SliceParams.h"

void plotbkgfits(bool fakedata = false){

  int colors[nIntTypes] = {2, 3, 5, 7, 33, 9, 46, 28, 41};

  TFile *fmc = TFile::Open("../install/bin/mcprod4a.root");
  TFile *fdata;
  TFile *fcorr;
  if (fakedata){
    fdata = TFile::Open("../install/bin/mcprod4a.root");
    fcorr = TFile::Open("../install/bin/backgroundfits_mc.root");
  }
  else{
    fdata = TFile::Open("../install/bin/data.root");
    fcorr = TFile::Open("../install/bin/backgroundfits.root");
  }

  TH1D *hmediandEdxSlice[nthinslices][nCuts][nIntTypes+1];
  TH1D *hdaughter_michel_scoreSlice[nthinslices][nCuts][nIntTypes+1];

  double totaldata = 0;
  double totalmc = 0;
  int nc = 0;

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

  TGraphErrors *corr = (TGraphErrors*)fcorr->Get("gr_corr_proton");

  for (int j = 0; j<nthinslices; ++j){
    TCanvas *can = new TCanvas(Form("can_%d",nc), Form("can_%d",nc), 800, 800);
    TPad *pad1 = new TPad(Form("pad1_%d",nc), Form("pad1_%d",nc), 0, 0.2, 1, 1.);
    pad1->SetBottomMargin(0.03);
    pad1->SetGridx();
    pad1->SetGridy();
    pad1->Draw();
    pad1->cd();
    TH1D *htotmc_orig;
    TH1D *htotmc_bestfit;
    THStack *hs = new THStack(Form("hs_%d",nc),"");
    for (int i = 0; i<nIntTypes; ++i){
      hmediandEdxSlice[j][kAPA3][i+1]->Scale(totaldata/totalmc);
      if (i==0) htotmc_orig = (TH1D*)hmediandEdxSlice[j][kAPA3][1]->Clone(Form("htotmc_orig_%d",nc));
      else htotmc_orig->Add(hmediandEdxSlice[j][kAPA3][i+1]);
      double correction = corr->GetPointY(j);
      if (i+1 == kMIDp) hmediandEdxSlice[j][kAPA3][i+1]->Scale(correction);
      if (i==0) htotmc_bestfit = (TH1D*)hmediandEdxSlice[j][kAPA3][1]->Clone(Form("htotmc_bestfit_%d",nc));
      else htotmc_bestfit->Add(hmediandEdxSlice[j][kAPA3][i+1]);
      hmediandEdxSlice[j][kAPA3][i+1]->SetFillColor(colors[i]);
      hmediandEdxSlice[j][kAPA3][i+1]->SetLineWidth(0);
      hs->Add(hmediandEdxSlice[j][kAPA3][i+1]);
    }
    
    double max = TMath::Max(hs->GetMaximum(), hmediandEdxSlice[j][kAPA3][0]->GetMaximum());
    hs->SetMaximum(1.2*max);
    hs->Draw("hist");
    hs->SetTitle(hmediandEdxSlice[j][kAPA3][0]->GetTitle());
    hs->GetXaxis()->SetTitle("Median dE/dx (MeV/cm)");
    hs->GetXaxis()->SetTitleSize(0);
    hs->GetXaxis()->SetLabelSize(0);
    hs->GetYaxis()->SetTitle("Events");
    hs->GetYaxis()->SetTitleSize(0.04);
    hs->GetYaxis()->SetTitleOffset(1.8);
    hs->GetYaxis()->SetLabelSize(0.04);
    hmediandEdxSlice[j][kAPA3][0]->SetMarkerStyle(20);
    hmediandEdxSlice[j][kAPA3][0]->Draw("same");
    htotmc_orig->SetLineStyle(2);
    htotmc_orig->Draw("hist same");
    TLegend *leg = new TLegend(0.15, 0.75, 0.88, 0.9);
    //  if (hmc[j][cut][0]->GetMaximumBin() < hmc[j][cut][0]->GetNbinsX()/2){
    //    leg = new TLegend(0.55,0.7,0.9,0.9);
    //  }
    //  else{
    //    leg = new TLegend(0.15,0.7,0.5,0.9);
    //  }
    //leg->SetFillStyle(0);
    leg->SetNColumns(3);
    leg->AddEntry(hmediandEdxSlice[j][kAPA3][0],Form("%s %.0f",intTypeName[0],hmediandEdxSlice[j][kAPA3][0]->Integral()),"ple");
    leg->AddEntry((TObject*)0,Form("Original MC %.0f",htotmc_orig->Integral()),"");
    leg->AddEntry((TObject*)0,Form("Best fit MC %.0f",htotmc_bestfit->Integral()),"");
    leg->AddEntry((TObject*)0,Form("Proton corr.: %.1f#pm%.1f",corr->GetPointY(j), corr->GetErrorY(j)),"");
    for (int i = 0; i<nIntTypes; ++i){
      leg->AddEntry(hmediandEdxSlice[j][kAPA3][i+1], Form("%s %.0f",intTypeName[i+1],hmediandEdxSlice[j][kAPA3][i+1]->Integral()), "f");
    }
    leg->Draw();


    if (fakedata){
      can->Print(Form("plots/canfitmc_%s_%d_%s_%02d.png", "hmediandEdxSlice", kAPA3, cutName[kAPA3], j));
      can->Print(Form("plots/canfitmc_%s_%d_%s_%02d.pdf", "hmediandEdxSlice", kAPA3, cutName[kAPA3], j));
      hs->SetMaximum(8*max);
      hs->SetMinimum(1);
      pad1->SetLogy();
      can->Print(Form("plots/canlogfitmc_%s_%d_%s_%02d.png", "hmediandEdxSlice", kAPA3, cutName[kAPA3], j));
      can->Print(Form("plots/canlogfitmc_%s_%d_%s_%02d.pdf", "hmediandEdxSlice", kAPA3, cutName[kAPA3], j));
    }
    else{
      can->Print(Form("plots/canfit_%s_%d_%s_%02d.png", "hmediandEdxSlice", kAPA3, cutName[kAPA3], j));
      can->Print(Form("plots/canfit_%s_%d_%s_%02d.pdf", "hmediandEdxSlice", kAPA3, cutName[kAPA3], j));
      hs->SetMaximum(8*max);
      hs->SetMinimum(1);
      pad1->SetLogy();
      can->Print(Form("plots/canlogfit_%s_%d_%s_%02d.png", "hmediandEdxSlice", kAPA3, cutName[kAPA3], j));
      can->Print(Form("plots/canlogfit_%s_%d_%s_%02d.pdf", "hmediandEdxSlice", kAPA3, cutName[kAPA3], j));
    }
    ++nc;
  }
}
