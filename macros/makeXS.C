{

  gStyle->SetOptStat(0);

  TFile *file = TFile::Open("../install/bin/mc.root");

  TH1D *eff_num_Int = (TH1D*)file->Get("eff_num_Int");
  TH1D *eff_den_Int = (TH1D*)file->Get("eff_den_Int");
  TH1D *eff_num_Inc = (TH1D*)file->Get("eff_num_Inc");
  TH1D *eff_den_Inc = (TH1D*)file->Get("eff_den_Inc");
  TH1D *pur_num_Int = (TH1D*)file->Get("pur_num_Int");
  TH1D *pur_num_Inc = (TH1D*)file->Get("pur_num_Inc");
  TH1D *pur_den = (TH1D*)file->Get("pur_den");

  TH1D *h_truesliceid_pion_uf = (TH1D*)file->Get("h_truesliceid_pion_uf");
  TH1D *h_truesliceid_pioninelastic_uf = (TH1D*)file->Get("h_truesliceid_pioninelastic_uf");

  TH1D *eff_Int = (TH1D*)eff_num_Int->Clone("eff_Int");
  eff_Int->Divide(eff_den_Int);

  TH1D *eff_Inc = (TH1D*)eff_num_Inc->Clone("eff_Inc");
  eff_Inc->Divide(eff_den_Inc);

  TH1D *pur_Int = (TH1D*)pur_num_Int->Clone("pur_Int");
  pur_Int->Divide(pur_den);

  TH1D *pur_Inc = (TH1D*)pur_num_Inc->Clone("pur_Inc");
  pur_Inc->Divide(pur_den);

  TH1D *h_truesliceid_pion_all = (TH1D*)file->Get("h_truesliceid_pion_all");
  TH1D *h_truesliceid_pion_cuts = (TH1D*)file->Get("h_truesliceid_pion_cuts");
  TH1D *h_truesliceid_pioninelastic_all = (TH1D*)file->Get("h_truesliceid_pioninelastic_all");
  TH1D *h_truesliceid_pioninelastic_cuts = (TH1D*)file->Get("h_truesliceid_pioninelastic_cuts");
  TH1D *h_recosliceid_allevts_cuts = (TH1D*)file->Get("h_recosliceid_allevts_cuts");
  TH1D *h_recosliceid_pion_cuts = (TH1D*)file->Get("h_recosliceid_pion_cuts");
  TH1D *h_recosliceid_pioninelastic_cuts = (TH1D*)file->Get("h_recosliceid_pioninelastic_cuts");

  TH1D *hinc = (TH1D*)h_recosliceid_allevts_cuts->Clone("hinc");
  TH1D *hint = (TH1D*)h_recosliceid_allevts_cuts->Clone("hint");

  hinc->Multiply(pur_Inc);
  TCanvas *c1 = new TCanvas("c1","c1");
  h_recosliceid_allevts_cuts->SetLineColor(3);
  h_recosliceid_allevts_cuts->SetMarkerColor(3);
  h_recosliceid_allevts_cuts->SetTitle("All Pions;Reco SliceID;Events");
  h_recosliceid_allevts_cuts->Draw();
  hinc->DrawCopy("same");
  h_recosliceid_pion_cuts->SetLineColor(2);
  h_recosliceid_pion_cuts->Draw("same hist");
  TLegend *leg1 = new TLegend(0.5,0.6,0.8,0.9);
  leg1->SetFillStyle(0);
  leg1->AddEntry(h_recosliceid_allevts_cuts,"Selected","ple");
  leg1->AddEntry(hinc, "Selected #times purity","ple");
  leg1->AddEntry(h_recosliceid_pion_cuts,"Selected true pions","l");
  leg1->Draw();

  TCanvas *c2 = new TCanvas("c2","c2");
  h_truesliceid_pion_uf->SetLineColor(4);
  h_truesliceid_pion_uf->SetMarkerColor(4);
  h_truesliceid_pion_uf->SetTitle("All Pions; True SliceID; Events");
  h_truesliceid_pion_uf->Draw();
  //hinc->SetLineColor(3);
  //hinc->SetMarkerColor(3);
  hinc->DrawCopy("same");
  h_truesliceid_pion_all->SetLineColor(2);
  h_truesliceid_pion_all->SetMarkerColor(2);
  h_truesliceid_pion_all->Draw("same hist");
  TLegend *leg2 = new TLegend(0.5,0.6,0.8,0.9);
  leg2->SetFillStyle(0);
  leg2->AddEntry(hinc, "Selected #times purity","ple");
  leg2->AddEntry(h_truesliceid_pion_uf,"Unfolded pions","ple");
  leg2->AddEntry(h_truesliceid_pion_all,"True pions","l");
  leg2->Draw();

  hint->Multiply(pur_Int);
  TCanvas *c3 = new TCanvas("c3","c3");
  h_recosliceid_allevts_cuts->SetLineColor(3);
  h_recosliceid_allevts_cuts->SetMarkerColor(3);
  h_recosliceid_allevts_cuts->SetTitle("Pion Inelastic Scatterings;Reco SliceID;Events");
  h_recosliceid_allevts_cuts->Draw();
  hint->DrawCopy("same");
  h_recosliceid_pion_cuts->SetLineColor(2);
  h_recosliceid_pion_cuts->Draw("same hist");
  TLegend *leg3 = new TLegend(0.5,0.6,0.8,0.9);
  leg3->SetFillStyle(0);
  leg3->AddEntry(h_recosliceid_allevts_cuts,"Selected","ple");
  leg3->AddEntry(hint, "Selected #times purity","ple");
  leg3->AddEntry(h_recosliceid_pion_cuts,"Selected true pions","l");
  leg3->Draw();

  TCanvas *c4 = new TCanvas("c4","c4");
  h_truesliceid_pion_uf->SetLineColor(4);
  h_truesliceid_pion_uf->SetMarkerColor(4);
  h_truesliceid_pion_uf->SetTitle("Pion Inelastic Scatterings; True SliceID; Events");
  h_truesliceid_pion_uf->Draw();
  //hint->SetLineColor(3);
  //hint->SetMarkerColor(3);
  hint->DrawCopy("same");
  h_truesliceid_pion_all->SetLineColor(2);
  h_truesliceid_pion_all->SetMarkerColor(2);
  h_truesliceid_pion_all->Draw("same hist");
  TLegend *leg4 = new TLegend(0.5,0.6,0.8,0.9);
  leg4->SetFillStyle(0);
  leg4->AddEntry(hint, "Selected #times purity","ple");
  leg4->AddEntry(h_truesliceid_pion_uf,"Unfolded pions","ple");
  leg4->AddEntry(h_truesliceid_pion_all,"True pions","l");
  leg4->Draw();


}
  
