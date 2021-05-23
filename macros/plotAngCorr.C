{

  TFile *file = TFile::Open("../install/bin/mcprod4a.root");

  TGraphErrors *gr_trueAngCorr = (TGraphErrors*)file->Get("gr_trueAngCorr");
  TGraphErrors *gr_recoAngCorr = (TGraphErrors*)file->Get("gr_recoAngCorr");
  TGraphErrors *gr_reco_trueAngCorr = (TGraphErrors*)file->Get("gr_reco_trueAngCorr");

  TCanvas *c1 = new TCanvas("c1","c1");
  gr_trueAngCorr->Draw("ape");
  gr_trueAngCorr->SetTitle("");
  gr_trueAngCorr->GetXaxis()->SetTitle("Slice ID");
  gr_trueAngCorr->GetXaxis()->SetRangeUser(-1, 23);
  gr_trueAngCorr->GetYaxis()->SetRangeUser(0.7, 1);
  gr_trueAngCorr->GetYaxis()->SetTitle("Angular correction");
  gr_recoAngCorr->SetLineColor(2);
  gr_recoAngCorr->SetMarkerColor(2);
  gr_recoAngCorr->Draw("pe");
  TLegend *leg = new TLegend(0.2,0.3,0.6,0.5);
  leg->SetFillStyle(0);
  leg->AddEntry(gr_trueAngCorr, "True angular correction", "pe");
  leg->AddEntry(gr_recoAngCorr, "Reco angular correction", "pe");
  leg->Draw();

  TCanvas *c2 = new TCanvas("c2","c2");
  gr_reco_trueAngCorr->Draw("ape");
  gr_reco_trueAngCorr->SetTitle("");
  gr_reco_trueAngCorr->GetXaxis()->SetTitle("Slice ID");
  gr_reco_trueAngCorr->GetXaxis()->SetRangeUser(-1, 23);
  gr_reco_trueAngCorr->GetYaxis()->SetTitle("Reco - True angular correction (MeV)");

  c1->Print("plots/AngCorr.pdf");
  c2->Print("plots/AngCorrdiff.pdf");

}
