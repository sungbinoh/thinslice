{

  TFile *file = TFile::Open("../install/bin/mcprod4a.root");

  TGraphErrors *gr_trueincE = (TGraphErrors*)file->Get("gr_trueincE");
  TGraphErrors *gr_recoincE = (TGraphErrors*)file->Get("gr_recoincE");
  TGraphErrors *gr_reco_trueincE = (TGraphErrors*)file->Get("gr_reco_trueincE");

  TCanvas *c1 = new TCanvas("c1","c1");
  gr_trueincE->Draw("ape");
  gr_trueincE->SetTitle("");
  gr_trueincE->GetXaxis()->SetTitle("Slice ID");
  gr_trueincE->GetXaxis()->SetRangeUser(-1, 23);
  gr_trueincE->GetYaxis()->SetTitle("Pion KE (MeV)");
  gr_recoincE->SetLineColor(2);
  gr_recoincE->SetMarkerColor(2);
  gr_recoincE->Draw("pe");
  TLegend *leg = new TLegend(0.5,0.65,0.9,0.9);
  leg->SetFillStyle(0);
  leg->AddEntry(gr_trueincE, "True KE", "pe");
  leg->AddEntry(gr_recoincE, "Reco KE", "pe");
  leg->Draw();

  TCanvas *c2 = new TCanvas("c2","c2");
  gr_reco_trueincE->Draw("ape");
  gr_reco_trueincE->SetTitle("");
  gr_reco_trueincE->GetXaxis()->SetTitle("Slice ID");
  gr_reco_trueincE->GetXaxis()->SetRangeUser(-1, 23);
  gr_reco_trueincE->GetYaxis()->SetTitle("Reco KE - True KE (MeV)");

  c1->Print("plots/incE.pdf");
  c2->Print("plots/incEdiff.pdf");

}
