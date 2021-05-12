{

#include "../SliceParams.h"

  gStyle->SetOptStat(0);

  TFile *file = TFile::Open("../install/bin/mcprod4.root");

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
  h_recosliceid_pioninelastic_cuts->SetLineColor(2);
  h_recosliceid_pioninelastic_cuts->Draw("same hist");
  TLegend *leg3 = new TLegend(0.5,0.6,0.8,0.9);
  leg3->SetFillStyle(0);
  leg3->AddEntry(h_recosliceid_allevts_cuts,"Selected","ple");
  leg3->AddEntry(hint, "Selected #times purity","ple");
  leg3->AddEntry(h_recosliceid_pioninelastic_cuts,"Selected true pions","l");
  leg3->Draw();

  TCanvas *c4 = new TCanvas("c4","c4");
  h_truesliceid_pioninelastic_uf->SetLineColor(4);
  h_truesliceid_pioninelastic_uf->SetMarkerColor(4);
  h_truesliceid_pioninelastic_uf->SetTitle("Pion Inelastic Scatterings; True SliceID; Events");
  h_truesliceid_pioninelastic_uf->Draw();
  //hint->SetLineColor(3);
  //hint->SetMarkerColor(3);
  hint->DrawCopy("same");
  h_truesliceid_pioninelastic_all->SetLineColor(2);
  h_truesliceid_pioninelastic_all->SetMarkerColor(2);
  h_truesliceid_pioninelastic_all->Draw("same hist");
  TLegend *leg4 = new TLegend(0.5,0.6,0.8,0.9);
  leg4->SetFillStyle(0);
  leg4->AddEntry(hint, "Selected #times purity","ple");
  leg4->AddEntry(h_truesliceid_pioninelastic_uf,"Unfolded pions","ple");
  leg4->AddEntry(h_truesliceid_pioninelastic_all,"True pions","l");
  leg4->Draw();

  double Ninc[nthinslices] = {0};
  double Nint[nthinslices] = {0};
  double err_inc[nthinslices] = {0};
  double err_int[nthinslices] = {0};

  for (int i = 0; i<nthinslices; ++i){
    Nint[i] = h_truesliceid_pioninelastic_uf->GetBinContent(i+2);
    err_int[i] = h_truesliceid_pioninelastic_uf->GetBinError(i+2);
    for (int j = i; j<=nthinslices; ++j){
      Ninc[i] += h_truesliceid_pion_uf->GetBinContent(j+2);
      err_inc[i] += pow(h_truesliceid_pion_uf->GetBinError(j+2),2);
      //cout<<i<<" "<<j<<" "<<h_truesliceid_pion_uf->GetBinContent(j+2)<<" "<<Ninc[i]<<endl;
    }
    err_inc[i] = sqrt(err_inc[i]);
  }

  double NA=6.02214076e23;
  double MAr=39.95; //gmol
  double Density = 1.39; // g/cm^3

  double xs[nthinslices] = {0};
  double err_xs[nthinslices] = {0};
  double incE[nthinslices] = {0};

  TGraphErrors *gr_trueincE = (TGraphErrors*)file->Get("gr_trueincE");
  TGraphErrors *gr_truexs = (TGraphErrors*)file->Get("gr_truexs");
  for (int i = 0; i<nthinslices; ++i){
    xs[i] = MAr/(Density*NA*thinslicewidth)*log(Ninc[i]/(Ninc[i]-Nint[i]))*1e27;
    //err_xs[i] = MAr/(Density*NA*thinslicewidth)*1e27*sqrt(N_int[i]+pow(N_int[i],2)/N_inc[i])/N_incidents[i];
    err_xs[i] = MAr/(Density*NA*thinslicewidth)*1e27*sqrt(pow(Nint[i]*err_inc[i]/Ninc[i]/(Ninc[i]-Nint[i]),2)+pow(err_int[i]/(Ninc[i]-Nint[i]),2));
    incE[i] = gr_trueincE->GetPointY(i);
    //std::cout<<i<<" "<<Ninc[i]<<" "<<Nint[i]<<" "<<xs[i]<<" "<<incE[i]<<std::endl;
  }

  TFile f2("../files/exclusive_xsec.root");
  TGraph *total_inel_KE = (TGraph*)f2.Get("total_inel_KE");
//  TGraph *abs_KE = (TGraph*)f2.Get("abs_KE");
//  TGraph *cex_KE = (TGraph*)f2.Get("cex_KE");

  TGraphErrors *gr_recoxs = new TGraphErrors(nthinslices, incE, xs, 0, err_xs);
  TCanvas *c5 = new TCanvas("c5");
  gr_recoxs->SetTitle("Pion Inelastic Cross Section");
  gr_recoxs->GetXaxis()->SetTitle("Pion Kinetic Energy (MeV)");
  gr_recoxs->GetYaxis()->SetTitle("#sigma_{inelastic} (mb)");
  gr_recoxs->Draw("ape");
  gr_truexs->SetMarkerColor(3);
  gr_truexs->SetLineColor(3);
  gr_truexs->Draw("pe");
  total_inel_KE->SetLineColor(2);
  total_inel_KE->Draw("c");
  TLegend *leg5 = new TLegend(0.3,0.6,0.8,0.9);
  leg5->SetFillStyle(0);
  leg5->AddEntry(gr_recoxs, "MC with reco information", "pe");
  leg5->AddEntry(gr_truexs, "MC with truth information", "pe");
  leg5->AddEntry(total_inel_KE, "Geant4 v4_10_6_p01c", "l");
  leg5->Draw();

  c5->Print("pi+inel.png");
}
  
