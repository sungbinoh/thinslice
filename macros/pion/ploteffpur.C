{
  gStyle->SetOptStat(0);

  TFile *file = TFile::Open("../../build/mcprod4a.root");

  TH1D *eff_num_Int = (TH1D*)file->Get("eff_num_Int");
  TH1D *eff_den_Int = (TH1D*)file->Get("eff_den_Int");
  TH1D *eff_num_Inc = (TH1D*)file->Get("eff_num_Inc");
  TH1D *eff_den_Inc = (TH1D*)file->Get("eff_den_Inc");
  TH1D *pur_num_Int = (TH1D*)file->Get("pur_num_Int");
  TH1D *pur_num_Inc = (TH1D*)file->Get("pur_num_Inc");
  TH1D *pur_den = (TH1D*)file->Get("pur_den");
  TH2D *response_SliceID_Int = (TH2D*)file->Get("response_SliceID_Int");
  TH2D *response_SliceID_Inc = (TH2D*)file->Get("response_SliceID_Inc");

  
  TCanvas *c1 = new TCanvas("c1","c1");
  TEfficiency *eff_pion = 0;
  if (TEfficiency::CheckConsistency(*eff_num_Inc, *eff_den_Inc)){
    eff_pion = new TEfficiency(*eff_num_Inc, *eff_den_Inc);
  }
  eff_pion->SetTitle("All Pions;True slice ID;Efficiency");
  eff_pion->Draw();

  /*for (int i=1; i<=eff_den_Int->GetNbinsX(); ++i){
    if (eff_den_Int->GetBinContent(i)<eff_num_Int->GetBinContent(i)){
      eff_num_Int->SetBinContent(i, eff_den_Int->GetBinContent(i));
      cout<<"$$$"<<i<<endl;
    }
  }*/
  TCanvas *c2 = new TCanvas("c2","c2");
  TEfficiency *eff_pioninel = 0;
  if (TEfficiency::CheckConsistency(*eff_num_Int, *eff_den_Int)){
    eff_pioninel = new TEfficiency(*eff_num_Int, *eff_den_Int);
  }
  eff_pioninel->SetTitle("Pion Inelastic Scatterings;True slice ID;Efficiency");
  eff_pioninel->Draw();

  /*TCanvas *c3 = new TCanvas("c3","c3");
  TEfficiency *pur_pion = 0;
  if (TEfficiency::CheckConsistency(*pur_num_Inc, *pur_den)){
    pur_pion = new TEfficiency(*pur_num_Inc, *pur_den);
  }
  pur_pion->SetTitle("All Pions;Reco slice ID;Purity");
  pur_pion->Draw();

  TCanvas *c4 = new TCanvas("c4","c4");
  TEfficiency *pur_pioninel = 0;
  if (TEfficiency::CheckConsistency(*pur_num_Int, *pur_den)){
    pur_pioninel = new TEfficiency(*pur_num_Int, *pur_den);
  }
  pur_pioninel->SetTitle("Pion Inelastic Scatterings;Reco slice ID;Purity");
  pur_pioninel->Draw();

  TCanvas *c5 = new TCanvas("c5","c5");
  response_SliceID_Int->Draw("colz");

  TCanvas *c6 = new TCanvas("c6","c6");
  response_SliceID_Inc->Draw("colz");

  c1->Print("plots/pioneff.png");
  c2->Print("plots/pionineleff.png");
  c3->Print("plots/pionpur.png");
  c4->Print("plots/pioninelpur.png");
  c5->Print("plots/pionres.png");
  c6->Print("plots/pioninelres.png");

  c1->Print("plots/pioneff.pdf");
  c2->Print("plots/pionineleff.pdf");
  c3->Print("plots/pionpur.pdf");
  c4->Print("plots/pioninelpur.pdf");
  c5->Print("plots/pionres.pdf");
  c6->Print("plots/pioninelres.pdf");*/

}
