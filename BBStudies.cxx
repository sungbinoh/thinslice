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
  //https://pdg.lbl.gov/2012/AtomicNuclearProperties/MUON_ELOSS_TABLES/muonloss_289.pdf
  const int np = 29;
  double spline_KE[np] = {
    10,    14,    20,    30,    40,     80,     100,    140,    200,   300,
    400,   800,   1000,  1400,  2000,   3000,   4000,   8000,   10000, 14000,
    20000, 30000, 40000, 80000, 100000, 140000, 200000, 300000, 400000};

  double spline_Range[np] = {
    9.833E-1, 1.786E0, 3.321E0, 6.598E0, 1.058E1, 3.084E1, 4.250E1, 6.732E1,
    1.063E2,  1.725E2, 2.385E2, 4.934E2, 6.163E2, 8.552E2, 1.202E3, 1.758E3,
    2.297E3,  4.359E3, 5.354E3, 7.298E3, 1.013E4, 1.469E4, 1.910E4, 3.558E4,
    4.326E4,  5.768E4, 7.734E4, 1.060E5, 1.307E5};

  for (int i = 0; i<np; ++i) spline_Range[i] /= 1.396;

  double bbRange[np];
  double bbRangeSpline[np];

  for (int i = 0; i<np; ++i){
    //cout<<spline_KE[i]<<" "<<bb.RangeFromKE(spline_KE[i])<<endl;
    bbRange[i] = bb.RangeFromKE(spline_KE[i]);
    bbRangeSpline[i] = bb.RangeFromKESpline(spline_KE[i]);
    //cout<<spline_KE[i]<<" "<<spline_Range[i]<<" "<<bbRange[i]<<" "<<bbRangeSpline[i]<<endl;
  }

  TGraph *grmuKERangePDG = new TGraph(np, spline_Range, spline_KE);
  TGraph *grmuKERangeBB = new TGraph(np, bbRange, spline_KE);
  TGraph *grmuKERangeBBSpline = new TGraph(np, bbRangeSpline, spline_KE);

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
  grmuKERangeBBSpline->SetMarkerStyle(24);
  grmuKERangeBBSpline->SetMarkerColor(3);
  grmuKERangeBBSpline->SetLineColor(3);
  grmuKERangeBBSpline->Draw("pc");
  gPad->SetLogx();
  gPad->SetLogy();
  TLegend *leg1 = new TLegend(0.2,0.65,0.5,0.9);
  leg1->SetFillStyle(0);
  leg1->AddEntry(grmuKERangePDG, "PDG","lp");
  leg1->AddEntry(grmuKERangeBB, "Bethe Bloch","lp");
  leg1->AddEntry(grmuKERangeBBSpline, "Bethe Bloch spline","lp");
  leg1->Draw();
  c1->Print("muKERange.png");
  c1->Print("muKERange.pdf");

  return 0;

}
