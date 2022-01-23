#include "TFile.h"
#include "TEfficiency.h"
#include "TH1D.h"
#include "TH2D.h"

TFile *file;

void ploteff(string name, string xtitle, string ytitle, string title, double max = 1.){

  TCanvas *c1 = new TCanvas("c1","c1");

  TEfficiency *eff = (TEfficiency*)file->Get(Form("h_Eff_%s",name.c_str()));
  eff->SetTitle(Form("%s;%s;%s", title.c_str(), xtitle.c_str(), ytitle.c_str()));
  eff->SetMarkerStyle(20);
  eff->SetMarkerSize(0.5);
  eff->Draw();
  gPad->Update();
  auto graph = eff->GetPaintedGraph();
  graph->SetMaximum(max);

  c1->Print(Form("plots/eff_%s.png", name.c_str()));
  c1->Print(Form("plots/eff_%s.pdf", name.c_str()));
  delete c1;

}

void plotres(string name){
  
  TCanvas *c1 = new TCanvas("c1","c1");

  TH2D *h2 = (TH2D*)file->Get(Form("h_reco_true_%s",name.c_str()));
  h2->Draw("colz");

  c1->Print(Form("plots/res2d_%s.png", name.c_str()));
  c1->Print(Form("plots/res2d_%s.pdf", name.c_str()));

  TH1D *h1 = (TH1D*)file->Get(Form("h_res_%s",name.c_str()));
  h1->Draw();
  h1->Fit("gaus");

  gPad->Update();
  TPaveStats *st = (TPaveStats*)h1->FindObject("stats");
  st->SetX1NDC(0.6); //new x start position
  //st->SetX2NDC(newx2); //new x end position

  c1->Print(Form("plots/res1d_%s.png", name.c_str()));
  c1->Print(Form("plots/res1d_%s.pdf", name.c_str()));

  delete c1;

}

void plot(){

  gSystem->Exec("rm -rf plots");
  gSystem->Exec("mkdir plots");

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(11);
  gErrorIgnoreLevel = kWarning;

  file = TFile::Open("../../job/effval.root");
  ploteff("Ppi","Momentum (MeV/c)", "Efficiency", "Changed pion");
  ploteff("Ppi_michel","Momentum (MeV/c)", "Efficiency", "Changed pion with Michel electron", 0.2);
  ploteff("thetapi","#theta (deg)", "Efficiency", "Changed pion");
  ploteff("thetapi_michel","#theta (deg)", "Efficiency", "Changed pion with Michel electron", 0.2);
  ploteff("phipi","#phi (deg)", "Efficiency", "Changed pion");
  ploteff("phipi_michel","#phi (deg)", "Efficiency", "Changed pion with Michel electron", 0.2);

  ploteff("Pp","Momentum (MeV/c)", "Efficiency", "Proton");
  ploteff("thetap","#theta (deg)", "Efficiency", "Proton", 0.6);
  ploteff("phip","#phi (deg)", "Efficiency", "Proton", 0.4);

  plotres("Ppi_sel");
  plotres("Ppi_michel");
  plotres("thetapi_sel");
  plotres("thetapi_michel");
  plotres("phipi_sel");
  plotres("phipi_michel");

  plotres("Pp_sel");
  plotres("thetap_sel");
  plotres("phip_sel");

}
