#include "../EventType.h"
#include "../EventSelection.h"
#include "../SliceParams.h"

int nc = 0;
double totaldata = 0;
double totalmc = 0;
int colors[nIntTypes] = {2, 3, 800, 7, 33, 9, 46, 28, 41};

void plotmediandEdx(TH1D *h[nCuts][nIntTypes+1], int j, double correction, double correction_err, bool fakedata){
  
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
    h[kAPA3][i+1]->Scale(totaldata/totalmc);
    if (i==0) htotmc_orig = (TH1D*)h[kAPA3][1]->Clone(Form("htotmc_orig_%d",nc));
    else htotmc_orig->Add(h[kAPA3][i+1]);
    if (i+1 == kMIDp) h[kAPA3][i+1]->Scale(correction);
    if (i==0) htotmc_bestfit = (TH1D*)h[kAPA3][1]->Clone(Form("htotmc_bestfit_%d",nc));
    else htotmc_bestfit->Add(h[kAPA3][i+1]);
    h[kAPA3][i+1]->SetFillColor(colors[i]);
    h[kAPA3][i+1]->SetLineWidth(0);
    hs->Add(h[kAPA3][i+1]);
  }
  
  double max = TMath::Max(hs->GetMaximum(), h[kAPA3][0]->GetMaximum());
  hs->SetMaximum(1.2*max);
  hs->Draw("hist");
  hs->SetTitle(h[kAPA3][0]->GetTitle());
  hs->GetXaxis()->SetTitle("Median dE/dx (MeV/cm)");
  hs->GetXaxis()->SetTitleSize(0);
  hs->GetXaxis()->SetLabelSize(0);
  hs->GetYaxis()->SetTitle("Events");
  hs->GetYaxis()->SetTitleSize(0.04);
  hs->GetYaxis()->SetTitleOffset(1.8);
  hs->GetYaxis()->SetLabelSize(0.04);
  h[kAPA3][0]->SetMarkerStyle(20);
  h[kAPA3][0]->Draw("same");
  htotmc_orig->SetLineStyle(2);
  htotmc_orig->Draw("hist same");
  TLegend *leg = new TLegend(0.15, 0.75, 0.88, 0.9);
  leg->SetFillStyle(0);
  leg->SetNColumns(3);
  if (fakedata){
    leg->AddEntry(h[kAPA3][0],Form("%s %.0f","fake data",h[kAPA3][0]->Integral()),"ple");
  }
  else{
    leg->AddEntry(h[kAPA3][0],Form("%s %.0f","data",h[kAPA3][0]->Integral()),"ple");
  }
  leg->AddEntry(htotmc_orig,Form("Original MC %.0f",htotmc_orig->Integral()),"l");
  leg->AddEntry((TObject*)0,Form("Best fit MC %.0f",htotmc_bestfit->Integral()),"");
  leg->AddEntry((TObject*)0,Form("Proton corr.: %.1f#pm%.1f",correction, correction_err),"");
  for (int i = 0; i<nIntTypes; ++i){
    leg->AddEntry(h[kAPA3][i+1], Form("%s %.0f",intTypeName[i+1],h[kAPA3][i+1]->Integral()), "f");
  }
  leg->Draw();
  
  can->cd();
  TPad *pad2 = new TPad(Form("pad2_%d",nc), Form("pad2_%d",nc), 0, 0, 1, 0.2);
  pad2->SetTopMargin(0.02);
  pad2->SetBottomMargin(0.25);
  pad2->SetGridx();
  pad2->SetGridy();
  pad2->Draw();
  pad2->cd();
  
  TH1D *hratio_bestfit = (TH1D*)h[kAPA3][0]->Clone(Form("hratio_bestfit_%d",nc));
  hratio_bestfit->SetTitle("");
  hratio_bestfit->Divide(htotmc_bestfit);
  hratio_bestfit->GetYaxis()->SetTitle("Data/MC");
  hratio_bestfit->GetXaxis()->SetLabelSize(0.12);
  hratio_bestfit->GetXaxis()->SetTitleSize(0.12);
  hratio_bestfit->GetXaxis()->SetTitleOffset(1.);
  hratio_bestfit->GetYaxis()->SetLabelSize(0.12);
  hratio_bestfit->GetYaxis()->SetTitleSize(0.12);
  hratio_bestfit->GetYaxis()->SetTitleOffset(.5);
  hratio_bestfit->GetYaxis()->SetRangeUser(0,5);
  hratio_bestfit->GetYaxis()->SetNdivisions(505);
  hratio_bestfit->Draw();
  TH1D *hratio_orig = (TH1D*)h[kAPA3][0]->Clone(Form("hratio_orig_%d",nc));
  hratio_orig->Divide(htotmc_orig);
  hratio_orig->SetLineStyle(2);
  hratio_orig->Draw("hist same");
  TLegend *leg2 = new TLegend(0.15,0.7,0.5,0.95);
  leg2->SetFillStyle(0);
  leg2->AddEntry(hratio_bestfit,"best fit","ple");
  leg2->AddEntry(hratio_orig,"original","l");
  leg2->Draw();
  
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

void plotdaughter_michel_score(TH1D *h[nCuts][nIntTypes+1], int j, double correction, double correction_err, bool fakedata){

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
    h[kAPA3][i+1]->Scale(totaldata/totalmc);
    if (i==0) htotmc_orig = (TH1D*)h[kAPA3][1]->Clone(Form("htotmc_orig_%d",nc));
    else htotmc_orig->Add(h[kAPA3][i+1]);
    if (i+1 == kMIDmu || i+1 == kMuon) h[kAPA3][i+1]->Scale(correction);
    if (i==0) htotmc_bestfit = (TH1D*)h[kAPA3][1]->Clone(Form("htotmc_bestfit_%d",nc));
    else htotmc_bestfit->Add(h[kAPA3][i+1]);
    h[kAPA3][i+1]->SetFillColor(colors[i]);
    h[kAPA3][i+1]->SetLineWidth(0);
    hs->Add(h[kAPA3][i+1]);
  }
  
  double max = TMath::Max(hs->GetMaximum(), h[kAPA3][0]->GetMaximum());
  hs->SetMaximum(1.2*max);
  hs->Draw("hist");
  hs->SetTitle(h[kAPA3][0]->GetTitle());
  hs->GetXaxis()->SetTitle("Daughter Michel Score");
  hs->GetXaxis()->SetTitleSize(0);
  hs->GetXaxis()->SetLabelSize(0);
  hs->GetYaxis()->SetTitle("Events");
  hs->GetYaxis()->SetTitleSize(0.04);
  hs->GetYaxis()->SetTitleOffset(1.8);
  hs->GetYaxis()->SetLabelSize(0.04);
  h[kAPA3][0]->SetMarkerStyle(20);
  h[kAPA3][0]->Draw("same");
  htotmc_orig->SetLineStyle(2);
  htotmc_orig->Draw("hist same");
  TLegend *leg = new TLegend(0.15, 0.75, 0.88, 0.9);
  leg->SetFillStyle(0);
  leg->SetNColumns(3);
  if (fakedata){
    leg->AddEntry(h[kAPA3][0],Form("%s %.0f","fake data",h[kAPA3][0]->Integral()),"ple");
  }
  else{
    leg->AddEntry(h[kAPA3][0],Form("%s %.0f","data",h[kAPA3][0]->Integral()),"ple");
  }
  leg->AddEntry(htotmc_orig,Form("Original MC %.0f",htotmc_orig->Integral()),"l");
  leg->AddEntry((TObject*)0,Form("Best fit MC %.0f",htotmc_bestfit->Integral()),"");
  leg->AddEntry((TObject*)0,Form("Muon corr.: %.1f#pm%.1f",correction, correction_err),"");
  for (int i = 0; i<nIntTypes; ++i){
    leg->AddEntry(h[kAPA3][i+1], Form("%s %.0f",intTypeName[i+1],h[kAPA3][i+1]->Integral()), "f");
  }
  leg->Draw();
  
  can->cd();
  TPad *pad2 = new TPad(Form("pad2_%d",nc), Form("pad2_%d",nc), 0, 0, 1, 0.2);
  pad2->SetTopMargin(0.02);
  pad2->SetBottomMargin(0.25);
  pad2->SetGridx();
  pad2->SetGridy();
  pad2->Draw();
  pad2->cd();
  
  TH1D *hratio_bestfit = (TH1D*)h[kAPA3][0]->Clone(Form("hratio_bestfit_%d",nc));
  hratio_bestfit->SetTitle("");
  hratio_bestfit->Divide(htotmc_bestfit);
  hratio_bestfit->GetXaxis()->SetLabelSize(0.12);
  hratio_bestfit->GetXaxis()->SetTitleSize(0.12);
  hratio_bestfit->GetXaxis()->SetTitleOffset(1.);
  hratio_bestfit->GetYaxis()->SetTitle("Data/MC");
  hratio_bestfit->GetYaxis()->SetLabelSize(0.12);
  hratio_bestfit->GetYaxis()->SetTitleSize(0.12);
  hratio_bestfit->GetYaxis()->SetTitleOffset(.5);
  hratio_bestfit->GetYaxis()->SetRangeUser(0,5);
  hratio_bestfit->GetYaxis()->SetNdivisions(505);
  hratio_bestfit->Draw();
  TH1D *hratio_orig = (TH1D*)h[kAPA3][0]->Clone(Form("hratio_orig_%d",nc));
  hratio_orig->Divide(htotmc_orig);
  hratio_orig->SetLineStyle(2);
  hratio_orig->Draw("hist same");
  TLegend *leg2 = new TLegend(0.15,0.7,0.5,0.95);
  leg2->SetFillStyle(0);
  leg2->AddEntry(hratio_bestfit,"best fit","ple");
  leg2->AddEntry(hratio_orig,"original","l");
  leg2->Draw();
  
  if (fakedata){
    can->Print(Form("plots/canfitmc_%s_%d_%s_%02d.png", "hdaughter_michel_scoreSlice", kAPA3, cutName[kAPA3], j));
    can->Print(Form("plots/canfitmc_%s_%d_%s_%02d.pdf", "hdaughter_michel_scoreSlice", kAPA3, cutName[kAPA3], j));
    hs->SetMaximum(8*max);
    hs->SetMinimum(1);
    pad1->SetLogy();
    can->Print(Form("plots/canlogfitmc_%s_%d_%s_%02d.png", "hdaughter_michel_scoreSlice", kAPA3, cutName[kAPA3], j));
    can->Print(Form("plots/canlogfitmc_%s_%d_%s_%02d.pdf", "hdaughter_michel_scoreSlice", kAPA3, cutName[kAPA3], j));
  }
  else{
    can->Print(Form("plots/canfit_%s_%d_%s_%02d.png", "hdaughter_michel_scoreSlice", kAPA3, cutName[kAPA3], j));
    can->Print(Form("plots/canfit_%s_%d_%s_%02d.pdf", "hdaughter_michel_scoreSlice", kAPA3, cutName[kAPA3], j));
    hs->SetMaximum(8*max);
    hs->SetMinimum(1);
    pad1->SetLogy();
    can->Print(Form("plots/canlogfit_%s_%d_%s_%02d.png", "hdaughter_michel_scoreSlice", kAPA3, cutName[kAPA3], j));
    can->Print(Form("plots/canlogfit_%s_%d_%s_%02d.pdf", "hdaughter_michel_scoreSlice", kAPA3, cutName[kAPA3], j));
  }
  ++nc;
}


void plotbkgfits(bool fakedata = false){
  
  gStyle->SetOptStat(0);
  gErrorIgnoreLevel = kWarning;

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

  TVectorD *fitresults = (TVectorD*)fcorr->Get("fitresults");

  TH1D *hmediandEdxSlice[nthinslices][nCuts][nIntTypes+1];
  TH1D *hdaughter_michel_scoreSlice[nthinslices][nCuts][nIntTypes+1];

  TH1D *hmediandEdx[nCuts][nIntTypes+1];
  TH1D *hdaughter_michel_score[nCuts][nIntTypes+1];

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
        if (k==0){
          if (j==0){
            hmediandEdx[i][j] = (TH1D*)fdata->Get(Form("hmediandEdx_%d_%d",i,j));
            hdaughter_michel_score[i][j] = (TH1D*)fdata->Get(Form("hdaughter_michel_score_%d_%d",i,j));
          }
          else{
            hmediandEdx[i][j] = (TH1D*)fmc->Get(Form("hmediandEdx_%d_%d",i,j));
            hdaughter_michel_score[i][j] = (TH1D*)fmc->Get(Form("hdaughter_michel_score_%d_%d",i,j));
          }
        }
      }
    }
  }
  std::cout<<totaldata<<" "<<totalmc<<std::endl;

  TGraphErrors *corr = (TGraphErrors*)fcorr->Get("gr_corr_proton");

  TCanvas *c1 = new TCanvas("c1","c1");
  corr->GetXaxis()->SetTitle("Slice ID");
  corr->GetYaxis()->SetTitle("Proton correction");
  corr->Draw();
  if (fakedata){
    corr->SetTitle("Fake data");
    c1->Print("plots/corr_proton_mc.png");
    c1->Print("plots/corr_proton_mc.pdf");
  }
  else{
    corr->SetTitle("Data");
    c1->Print("plots/corr_proton.png");
    c1->Print("plots/corr_proton.pdf");
  }

  for (int j = 0; j<nthinslices; ++j){

    double correction = corr->GetPointY(j);
    double correction_err = corr->GetErrorY(j);
    plotmediandEdx(hmediandEdxSlice[j], j, correction, correction_err, fakedata);
  }


  corr = (TGraphErrors*)fcorr->Get("gr_corr_muon");

  TCanvas *c2 = new TCanvas("c2","c2");
  corr->GetXaxis()->SetTitle("Slice ID");
  corr->GetYaxis()->SetTitle("Muon correction");
  corr->Draw();
  if (fakedata){
    corr->SetTitle("Fake data");
    c2->Print("plots/corr_muon_mc.png");
    c2->Print("plots/corr_muon_mc.pdf");
  }
  else{
    corr->SetTitle("Data");
    c2->Print("plots/corr_muon.png");
    c2->Print("plots/corr_muon.pdf");
  }

  for (int j = 0; j<nthinslices; ++j){

    double correction = corr->GetPointY(j);
    double correction_err = corr->GetErrorY(j);
    plotdaughter_michel_score(hdaughter_michel_scoreSlice[j], j, correction, correction_err, fakedata);
  }

  for (int i = 0; i<4; ++i){
    cout<<(*fitresults)[i]<<endl;
  }

  plotmediandEdx(hmediandEdx, -1, (*fitresults)[0], (*fitresults)[1], fakedata);

  plotdaughter_michel_score(hdaughter_michel_score, -1, (*fitresults)[2], (*fitresults)[3], fakedata);
  
}
