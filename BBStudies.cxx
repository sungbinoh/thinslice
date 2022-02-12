#include "BetheBloch.h"
#include "util.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include <iostream>

using namespace std;

int main(){

  SetProtoDUNEStyle();
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);

  BetheBloch bb(13);

  //cout<<bb.meandEdx(1000)<<" "<<bb.MPVdEdx(1000, 0.5)<<endl;

  const int np = 13;
  double spline_KE[np] = {10, 14, 20, 30, 40, 80, 100, 140, 200, 300, 400, 800, 1000};
  double spline_Range[np] = {0.70437, 1.27937, 2.37894, 4.72636, 7.5788, 22.0917, 30.4441, 48.2235, 76.1461, 123.567, 170.845, 353.438, 441.476};
  double bbRange[np];

  for (int i = 0; i<np; ++i){
    //cout<<spline_KE[i]<<" "<<bb.RangeFromKE(spline_KE[i])<<endl;
    bbRange[i] = bb.RangeFromKE(spline_KE[i]);
  }

  TGraph *grmuKERangePDG = new TGraph(np, spline_Range, spline_KE);
  TGraph *grmuKERangeBB = new TGraph(np, bbRange, spline_KE);

  TCanvas *c1 = new TCanvas("c1","c1");
  grmuKERangePDG->SetMarkerStyle(24);
  grmuKERangePDG->SetMarkerColor(2);
  grmuKERangePDG->SetLineColor(2);
  grmuKERangePDG->Draw("apc");
  grmuKERangePDG->SetTitle("Muon");
  grmuKERangePDG->GetXaxis()->SetTitle("Range (cm)");
  grmuKERangePDG->GetYaxis()->SetTitle("KE (MeV)");
  grmuKERangeBB->SetMarkerStyle(24);
  grmuKERangeBB->SetMarkerColor(4);
  grmuKERangeBB->SetLineColor(4);
  grmuKERangeBB->Draw("pc");
  gPad->SetLogx();
  gPad->SetLogy();
  TLegend *leg1 = new TLegend(0.2,0.7,0.5,0.9);
  leg1->SetFillStyle(0);
  leg1->AddEntry(grmuKERangePDG, "PDG","lp");
  leg1->AddEntry(grmuKERangeBB, "Bethe Bloch","lp");
  leg1->Draw();
  c1->Print("muKERange.png");
  c1->Print("muKERange.pdf");

  return 0;

}
