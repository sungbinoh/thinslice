{

  #include "../SliceParams.h"

  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.05);
  gStyle->SetPadRightMargin(0.05);


  TCanvas *c1 = new TCanvas("c1","c1",800,500);

  TH1D *h = new TH1D("h",";z (cm)", (nthinslices + 2), -10, (nthinslices+1)*10);
  for (int i = 1; i<=h->GetNbinsX(); ++i){
    h->SetBinContent(i,1);
  }
  h->SetFillColor(38);
  h->SetBarWidth(.9);
  h->Draw("b");
  h->GetYaxis()->SetRangeUser(0, 1.);
  h->GetYaxis()->SetLabelSize(0);
  h->GetYaxis()->SetTickLength(0);
  h->GetXaxis()->CenterTitle();
  gPad->RedrawAxis();
  
  TLatex l;
  l.SetTextAngle(90);
  l.SetTextColor(0);
  for (int i = 1; i<=h->GetNbinsX(); ++i){
    l.DrawLatex(h->GetBinCenter(i)+2, 0.4, Form("Slice %d",i-2));
  }

  c1->Print("plots/slices.pdf");
}
