{

  TFile *file = TFile::Open("../install/bin/hadana.root");

  TH1D *h_truesliceid_pion_all = (TH1D*)file->Get("h_truesliceid_pion_all");
  TH1D *h_truesliceid_pion_cuts = (TH1D*)file->Get("h_truesliceid_pion_cuts");
  TH1D *h_truesliceid_pioninelastic_all = (TH1D*)file->Get("h_truesliceid_pioninelastic_all");
  TH1D *h_truesliceid_pioninelastic_cuts = (TH1D*)file->Get("h_truesliceid_pioninelastic_cuts");
  TH1D *h_recosliceid_allevts_cuts = (TH1D*)file->Get("h_recosliceid_allevts_cuts");
  TH1D *h_recosliceid_pion_cuts = (TH1D*)file->Get("h_recosliceid_pion_cuts");
  TH1D *h_recosliceid_pioninelastic_cuts = (TH1D*)file->Get("h_recosliceid_pioninelastic_cuts");

  TCanvas *c1 = new TCanvas("c1","c1");
  TEfficiency *eff_pion = 0;
  if (TEfficiency::CheckConsistency(*h_truesliceid_pion_cuts, *h_truesliceid_pion_all)){
    eff_pion = new TEfficiency(*h_truesliceid_pion_cuts, *h_truesliceid_pion_all);
  }
  eff_pion->SetTitle(";True slice ID;Pion efficiency");
  eff_pion->Draw();

  TCanvas *c2 = new TCanvas("c2","c2");
  TEfficiency *eff_pioninel = 0;
  if (TEfficiency::CheckConsistency(*h_truesliceid_pioninelastic_cuts, *h_truesliceid_pioninelastic_all)){
    eff_pioninel = new TEfficiency(*h_truesliceid_pioninelastic_cuts, *h_truesliceid_pioninelastic_all);
  }
  eff_pioninel->SetTitle(";True slice ID;Pion inelastic scattering efficiency");
  eff_pioninel->Draw();

  TCanvas *c3 = new TCanvas("c3","c3");
  TEfficiency *pur_pion = 0;
  if (TEfficiency::CheckConsistency(*h_recosliceid_pion_cuts, *h_recosliceid_allevts_cuts)){
    pur_pion = new TEfficiency(*h_recosliceid_pion_cuts, *h_recosliceid_allevts_cuts);
  }
  pur_pion->SetTitle(";Reco slice ID;Pion purity");
  pur_pion->Draw();

  TCanvas *c4 = new TCanvas("c4","c4");
  TEfficiency *pur_pioninel = 0;
  if (TEfficiency::CheckConsistency(*h_recosliceid_pioninelastic_cuts, *h_recosliceid_allevts_cuts)){
    pur_pioninel = new TEfficiency(*h_recosliceid_pioninelastic_cuts, *h_recosliceid_allevts_cuts);
  }
  pur_pioninel->SetTitle(";Reco slice ID;Pion inelastic scattering purity");
  pur_pioninel->Draw();

}
