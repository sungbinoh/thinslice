#include "../../SliceParams.h"

void plotXS(){

  gStyle->SetOptStat(0);

  TFile *file = TFile::Open("../../build/XSMC.root");
  
  TH1D *h_sel_data = (TH1D*)file->Get("hdata");
  TH1D *h_sel_sig_inc = (TH1D*)file->Get("hsiginc");
  TH1D *hval_sel_sig_inc = (TH1D*)file->Get("hval_siginc_reco");
  TH1D *h_sel_sig_inc_uf = (TH1D*)file->Get("hsiginc_uf");
  TH1D *hval_sel_sig_inc_uf = (TH1D*)file->Get("hval_trueinc");
  TH1D *h_sel_sig_int = (TH1D*)file->Get("hsignal");
  TH1D *hval_sel_sig_int = (TH1D*)file->Get("hval_signal_reco");
  TH1D *h_sel_sig_int_uf = (TH1D*)file->Get("hsignal_uf");
  TH1D *hval_sel_sig_int_uf = (TH1D*)file->Get("hval_trueint");
  TGraphErrors *gr_inc = (TGraphErrors*)file->Get("gr_inc");
  TGraphErrors *gr_int = (TGraphErrors*)file->Get("gr_int");
  
  
  TCanvas *c1 = new TCanvas("c1","c1");
  h_sel_data->SetLineColor(3);
  h_sel_data->SetMarkerColor(3);
  h_sel_data->SetTitle("All Pions;Reco SliceID;Events");
  h_sel_data->DrawCopy();
  h_sel_sig_inc->DrawCopy("same");
  hval_sel_sig_inc->SetLineColor(2);
  hval_sel_sig_inc->Draw("same hist");
  TLegend *leg1 = new TLegend(0.5,0.6,0.8,0.9);
  leg1->SetFillStyle(0);
  leg1->AddEntry(h_sel_data,"Selected","ple");
  leg1->AddEntry(h_sel_sig_inc, "Selected #times purity","ple");
  leg1->AddEntry(hval_sel_sig_inc,"Selected true pions","l");
  leg1->Draw();

  TCanvas *c2 = new TCanvas("c2","c2");
  h_sel_sig_inc_uf->SetLineColor(4);
  h_sel_sig_inc_uf->SetMarkerColor(4);
  h_sel_sig_inc_uf->SetTitle("All Pions; True SliceID; Events");
  h_sel_sig_inc_uf->Draw();
  h_sel_sig_inc->DrawCopy("same");
  hval_sel_sig_inc_uf->SetLineColor(2);
  hval_sel_sig_inc_uf->SetMarkerColor(2);
  hval_sel_sig_inc_uf->Draw("same hist");
  TLegend *leg2 = new TLegend(0.5,0.6,0.8,0.9);
  leg2->SetFillStyle(0);
  leg2->AddEntry(h_sel_sig_inc, "Selected #times purity","ple");
  leg2->AddEntry(h_sel_sig_inc_uf,"Unfolded pions","ple");
  leg2->AddEntry(hval_sel_sig_inc_uf,"True pions","l");
  leg2->Draw();

  TCanvas *c3 = new TCanvas("c3","c3"); // compared to c1, only miss a few pion elastic events
  h_sel_data->SetLineColor(3);
  h_sel_data->SetMarkerColor(3);
  h_sel_data->SetTitle("Pion Inelastic Scatterings;Reco SliceID;Events");
  h_sel_data->Draw();
  h_sel_sig_int->DrawCopy("same");
  hval_sel_sig_int->SetLineColor(2);
  hval_sel_sig_int->Draw("same hist");
  TLegend *leg3 = new TLegend(0.5,0.6,0.8,0.9);
  leg3->SetFillStyle(0);
  leg3->AddEntry(h_sel_data,"Selected data","ple");
  leg3->AddEntry(h_sel_sig_int, "After bkg subtraction","ple");
  leg3->AddEntry(hval_sel_sig_int,"Selected true pions","l");
  leg3->Draw();

  TCanvas *c4 = new TCanvas("c4","c4");
  h_sel_sig_int_uf->SetLineColor(4);
  h_sel_sig_int_uf->SetMarkerColor(4);
  h_sel_sig_int_uf->SetTitle("Pion Inelastic Scatterings; True SliceID; Events");
  h_sel_sig_int_uf->Draw();
  h_sel_sig_int->DrawCopy("same");
  hval_sel_sig_int_uf->SetLineColor(2);
  hval_sel_sig_int_uf->SetMarkerColor(2);
  hval_sel_sig_int_uf->Draw("same hist");
  TLegend *leg4 = new TLegend(0.5,0.6,0.8,0.9);
  leg4->SetFillStyle(0);
  leg4->AddEntry(h_sel_sig_int, "Selected signal","ple");
  leg4->AddEntry(h_sel_sig_int_uf,"After unfolding","ple");
  leg4->AddEntry(hval_sel_sig_int_uf,"True pions","l");
  leg4->Draw();

  TCanvas *c6 = new TCanvas("c6","c6");
  gr_inc->SetTitle("");
  gr_inc->SetLineWidth(2);
  gr_inc->SetLineColor(4);
  gr_inc->SetMarkerColor(4);
  gr_inc->GetXaxis()->SetTitle("Slice ID");
  gr_inc->GetYaxis()->SetTitle("N_{Inc}");
  gr_inc->GetXaxis()->SetRangeUser(0, 23);
  gr_inc->Draw("ape");

  TCanvas *c7 = new TCanvas("c7","c7");
  gr_int->SetTitle("");
  gr_int->SetLineWidth(2);
  gr_int->SetLineColor(4);
  gr_int->SetMarkerColor(4);
  gr_int->GetXaxis()->SetTitle("Slice ID");
  gr_int->GetYaxis()->SetTitle("N_{Int}");
  gr_int->GetXaxis()->SetRangeUser(0, 23);
  gr_int->Draw("ape");


  double xs[pi::nthinslices] = {0};
  double err_xs[pi::nthinslices] = {0};
  double incE[pi::nthinslices] = {0};
  double err_incE[pi::nthinslices] = {0};
  
  TGraph *xs_new = (TGraph*)file->Get("gr_recoxs");
  TGraphErrors *gr_trueincE = (TGraphErrors*)file->Get("gr_trueincE");
  TGraphErrors *gr_truexs = (TGraphErrors*)file->Get("gr_truexs");
  for (int i = 0; i<pi::nthinslices; ++i){
    xs[i] = xs_new->GetPointY(i);
    err_xs[i] = xs_new->GetErrorY(i);
    incE[i] = gr_trueincE->GetPointY(i);
    err_incE[i] = gr_trueincE->GetErrorY(i);
  }
  TFile f2("../../files/exclusive_xsec.root");
  TGraph *total_inel_KE = (TGraph*)f2.Get("total_inel_KE");

  TGraphErrors *gr_recoxs = new TGraphErrors(pi::nthinslices, incE, xs, 0, err_xs);
  TCanvas *c5 = new TCanvas("c5", "c5");
  gr_recoxs->SetTitle("Pion Inelastic Cross Section");
  gr_recoxs->GetXaxis()->SetTitle("Pion Kinetic Energy (MeV)");
  gr_recoxs->GetXaxis()->SetRangeUser(360, 900);
  gr_recoxs->GetYaxis()->SetTitle("#sigma_{inelastic} (mb)");
  gr_recoxs->GetYaxis()->SetRangeUser(400, 900);
  gr_recoxs->SetLineWidth(2);
  gr_recoxs->Draw("ape");
  gr_truexs->SetMarkerColor(3);
  gr_truexs->SetLineColor(3);
  gr_truexs->Draw("pe");
  total_inel_KE->SetLineColor(2);
  total_inel_KE->Draw("c");
  TLegend *leg5 = new TLegend(0.3,0.65,0.8,0.9);
  leg5->SetFillStyle(0);
  leg5->AddEntry(gr_recoxs, "MC with reconstruction", "pe");
  leg5->AddEntry(gr_truexs, "MC truth", "pe");
  leg5->AddEntry(total_inel_KE, "Geant4 (theory prediction)", "l");
  leg5->Draw();

  c1->Print("plots/xs_sliceidinc_reco.pdf");
  c2->Print("plots/xs_sliceidinc_true.pdf");
  c3->Print("plots/xs_sliceidint_reco.pdf");
  c4->Print("plots/xs_sliceidint_true.pdf");
  c5->Print("plots/xs_pi+inel.pdf");
  c6->Print("plots/xs_Ninc.pdf");
  c7->Print("plots/xs_Nint.pdf");

  c1->Print("plots/xs_sliceidinc_reco.png");
  c2->Print("plots/xs_sliceidinc_true.png");
  c3->Print("plots/xs_sliceidint_reco.png");
  c4->Print("plots/xs_sliceidint_true.png");
  c5->Print("plots/xs_pi+inel.png");
  c6->Print("plots/xs_Ninc.png");
  c7->Print("plots/xs_Nint.png");

}
  
