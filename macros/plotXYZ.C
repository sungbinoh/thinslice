{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TFile *fmc = TFile::Open("../install/bin/mcprod4a.root");
  TFile *fdata = TFile::Open("../install/bin/data.root");

  TH1D *hreco_beam_startX_SCE_data = (TH1D*)fdata->Get("hreco_beam_startX_SCE_2_0");
  TH1D *hreco_beam_startX_SCE_mc = (TH1D*)fmc->Get("hreco_beam_startX_SCE_2_0");
  TH1D *hreco_beam_startY_SCE_data = (TH1D*)fdata->Get("hreco_beam_startY_SCE_2_0");
  TH1D *hreco_beam_startY_SCE_mc = (TH1D*)fmc->Get("hreco_beam_startY_SCE_2_0");
  TH1D *hreco_beam_startZ_SCE_data = (TH1D*)fdata->Get("hreco_beam_startZ_SCE_2_0");
  TH1D *hreco_beam_startZ_SCE_mc = (TH1D*)fmc->Get("hreco_beam_startZ_SCE_2_0");

  hreco_beam_startX_SCE_data->Scale(1./hreco_beam_startX_SCE_data->Integral());
  hreco_beam_startX_SCE_mc->Scale(1./hreco_beam_startX_SCE_mc->Integral());
  hreco_beam_startY_SCE_data->Scale(1./hreco_beam_startY_SCE_data->Integral());
  hreco_beam_startY_SCE_mc->Scale(1./hreco_beam_startY_SCE_mc->Integral());
  hreco_beam_startZ_SCE_data->Scale(1./hreco_beam_startZ_SCE_data->Integral());
  hreco_beam_startZ_SCE_mc->Scale(1./hreco_beam_startZ_SCE_mc->Integral());

  TCanvas *c1 = new TCanvas("c1","c1");
  hreco_beam_startX_SCE_data->SetMaximum(1.1*TMath::Max(hreco_beam_startX_SCE_data->GetMaximum(), hreco_beam_startX_SCE_mc->GetMaximum()));
  hreco_beam_startX_SCE_data->Draw();
  hreco_beam_startX_SCE_data->Fit("gaus","RQ","",-60,10);
  TF1 *fun = (TF1*)hreco_beam_startX_SCE_data->FindObject("gaus");
  fun->SetLineColor(1);
  std::cout<<"startX_data "<<fun->GetParameter(1)<<"+/-"<<fun->GetParameter(2)<<std::endl;
  hreco_beam_startX_SCE_mc->SetLineColor(2);
  hreco_beam_startX_SCE_mc->SetMarkerColor(2);
  hreco_beam_startX_SCE_mc->Draw("same");
  hreco_beam_startX_SCE_mc->Fit("gaus","RQ","",-60,10);
  fun = (TF1*)hreco_beam_startX_SCE_mc->FindObject("gaus");
  std::cout<<"startX_mc "<<fun->GetParameter(1)<<"+/-"<<fun->GetParameter(2)<<std::endl;
  TLegend *leg = new TLegend(0.6,0.7,0.9,0.9);
  leg->SetFillStyle(0);
  leg->AddEntry(hreco_beam_startX_SCE_data, "Data","ple");
  leg->AddEntry(hreco_beam_startX_SCE_mc, "MC","ple");
  leg->Draw();

  TCanvas *c2 = new TCanvas("c2","c2");
  hreco_beam_startY_SCE_data->SetMaximum(1.1*TMath::Max(hreco_beam_startY_SCE_data->GetMaximum(), hreco_beam_startY_SCE_mc->GetMaximum()));
  hreco_beam_startY_SCE_data->Draw();
  hreco_beam_startY_SCE_data->Fit("gaus","RQ","",400,480);
  fun = (TF1*)hreco_beam_startY_SCE_data->FindObject("gaus");
  fun->SetLineColor(1);
  std::cout<<"startY_data "<<fun->GetParameter(1)<<"+/-"<<fun->GetParameter(2)<<std::endl;
  hreco_beam_startY_SCE_mc->SetLineColor(2);
  hreco_beam_startY_SCE_mc->SetMarkerColor(2);
  hreco_beam_startY_SCE_mc->Draw("same");
  hreco_beam_startY_SCE_mc->Fit("gaus","RQ","",400,480);
  fun = (TF1*)hreco_beam_startY_SCE_mc->FindObject("gaus");
  std::cout<<"startY_mc "<<fun->GetParameter(1)<<"+/-"<<fun->GetParameter(2)<<std::endl;
  leg->Draw();

  TCanvas *c3 = new TCanvas("c3","c3");
  hreco_beam_startZ_SCE_data->SetMaximum(1.1*TMath::Max(hreco_beam_startZ_SCE_data->GetMaximum(), hreco_beam_startZ_SCE_mc->GetMaximum()));
  hreco_beam_startZ_SCE_data->Draw();
  hreco_beam_startZ_SCE_data->Fit("gaus","RQ","",0,8);
  fun = (TF1*)hreco_beam_startZ_SCE_data->FindObject("gaus");
  fun->SetLineColor(1);
  std::cout<<"startZ_data "<<fun->GetParameter(1)<<"+/-"<<fun->GetParameter(2)<<std::endl;
  hreco_beam_startZ_SCE_mc->SetLineColor(2);
  hreco_beam_startZ_SCE_mc->SetMarkerColor(2);
  hreco_beam_startZ_SCE_mc->Draw("same");
  hreco_beam_startZ_SCE_mc->Fit("gaus","RQ","",-2,2);
  fun = (TF1*)hreco_beam_startZ_SCE_mc->FindObject("gaus");
  std::cout<<"startZ_mc "<<fun->GetParameter(1)<<"+/-"<<fun->GetParameter(2)<<std::endl;
  leg->Draw();

  c1->Print("plots/startX.pdf");
  c2->Print("plots/startY.pdf");
  c3->Print("plots/startZ.pdf");

}
