#include "util.h"
#include "TVector3.h"
#include "TMath.h"
#include "TStyle.h"
#include "TGaxis.h"
#include <iostream>

using namespace std;

double CalWeight(const anavar & evt, const int &partype){
  double weight = 1.;
  
  if (partype == 3) // kMuon
    weight *= 1.;//58;
  
  /*if (partype == 0) { // fake data
    if (evt.true_beam_PDG == -13)
      weight = 1.6;
  }*/
  
  return weight;
}

void FillHistVec1D(TH1D *hist[pi::nIntTypes+1], const double &value, const int &partype, double weight, bool fill_underflow, bool fill_overflow){
  //hist[0]->Fill(value);
  if (partype>=0 && partype < pi::nIntTypes+1){
    if (value < hist[partype]->GetXaxis()->GetXmin()){ // underflow values
      if (fill_underflow) hist[partype]->Fill(hist[partype]->GetXaxis()->GetXmin(), weight);
    }
    else if (value < hist[partype]->GetXaxis()->GetXmax()){
      hist[partype]->Fill(value, weight);
    }
    else{ // overflow values
      if (fill_overflow) hist[partype]->Fill(hist[partype]->GetXaxis()->GetXmax()-0.000001, weight);
    }
  }
}

void FillHistVec2D(TH2D *hist[pi::nIntTypes+1], const double &value1, const double &value2, const int &partype, double weight){
  //hist[0]->Fill(value1, value2);
  if (partype>=0 && partype < pi::nIntTypes+1){
    hist[partype]->Fill(value1, value2, weight);
  }
}

void FillHist1D(TH1D *hist, const double &value, const double &wei){
  if (value<hist->GetXaxis()->GetXmin()){
    hist->Fill(hist->GetXaxis()->GetXmin(), wei);
  }
  else if (value<hist->GetXaxis()->GetXmax()){
    hist->Fill(value, wei);
  }
  else{
    hist->Fill(hist->GetXaxis()->GetXmax()-0.000001, wei);
  }
}

double GetPionKE(double length){
  if (length>0){
    // https://arxiv.org/abs/1306.1712
    return (8/(1-0.37))*pow(length,1-0.37);
  }
  else{
    cout<<"Invalid length "<<length<<", returning -1"<<endl;
    return -1;
  }
}

double GetTheta(double x, double y, double z){
  TVector3 v(x,y,z);
  if (!v.Mag()){
    cout<<"Empty vector, returning -1"<<endl;
    return -1;
  }
  else{
    return v.Theta()*180/TMath::Pi();
  }
}

double GetPhi(double x, double y, double z){
  TVector3 v(x,y,z);
  if (!v.Mag()){
    cout<<"Empty vector, returning -1"<<endl;
    return -1;
  }
  else{
    return v.Phi()*180/TMath::Pi();
  }
}

double GetThetaxz(double x, double y, double z){
  TVector3 v(x,y,z);
  if (!v.Mag()){
    cout<<"Empty vector, returning -1"<<endl;
    return -1;
  }
  else{
    return atan2(v.X(), v.Z())*180/TMath::Pi();
  }
}

double GetThetayz(double x, double y, double z){
  TVector3 v(x,y,z);
  if (!v.Mag()){
    cout<<"Empty vector, returning -1"<<endl;
    return -1;
  }
  else{
    return atan2(v.Y(), v.Z())*180/TMath::Pi();
  }
}

void SetProtoDUNEStyle(){

  // This is the protoDUNE style file

  TStyle* protoDUNEStyle = new  TStyle("protoDUNEStyle", "ProtoDUNE Style");

  // Colors

  //set the background color to white
  protoDUNEStyle->SetFillColor(10);
  protoDUNEStyle->SetFrameFillColor(10);
  protoDUNEStyle->SetCanvasColor(10);
  protoDUNEStyle->SetPadColor(10);
  protoDUNEStyle->SetTitleFillColor(10);
  protoDUNEStyle->SetStatColor(10);

  //dont put a colored frame around the plots
  protoDUNEStyle->SetFrameBorderMode(0);
  protoDUNEStyle->SetCanvasBorderMode(0);
  protoDUNEStyle->SetPadBorderMode(0);

  //use the primary color palette
  protoDUNEStyle->SetPalette(1,0);

  //set the default line color for a histogram to be black
  protoDUNEStyle->SetHistLineColor(kBlack);

  //set the default line color for a fit function to be red
  protoDUNEStyle->SetFuncColor(kRed);

  //make the axis labels black
  protoDUNEStyle->SetLabelColor(kBlack,"xyz");

  //set the default title color to be black
  protoDUNEStyle->SetTitleColor(kBlack);

  // Sizes

  //set the margins
  protoDUNEStyle->SetPadBottomMargin(0.125);
  protoDUNEStyle->SetPadTopMargin(0.075);
  protoDUNEStyle->SetPadLeftMargin(0.12);
  protoDUNEStyle->SetPadRightMargin(0.1);

  //set axis label and title text sizes
  protoDUNEStyle->SetLabelSize(0.045,"xy");
  protoDUNEStyle->SetLabelSize(0.035,"z");
  protoDUNEStyle->SetLabelOffset(0.005,"xy");
  //protoDUNEStyle->SetLabelOffset(0.005,"z");
  protoDUNEStyle->SetTitleSize(0.05,"xyz");
  protoDUNEStyle->SetTitleOffset(1.15,"x");
  protoDUNEStyle->SetTitleOffset(1.15,"yz");
  protoDUNEStyle->SetStatFontSize(0.05);
  protoDUNEStyle->SetTextSize(0.05);
  protoDUNEStyle->SetTitleBorderSize(0);
  protoDUNEStyle->SetStatBorderSize(0);

  //set line widths
  protoDUNEStyle->SetHistLineWidth(3);
  protoDUNEStyle->SetFrameLineWidth(2);
  protoDUNEStyle->SetFuncWidth(2);

  // Misc

  //align the titles to be centered
  //protoDUNEStyle->SetTitleAlign(22);

  //set the number of divisions to show
  protoDUNEStyle->SetNdivisions(506, "xy");

  //turn off xy grids
  protoDUNEStyle->SetPadGridX(0);
  protoDUNEStyle->SetPadGridY(0);

  //set the tick mark style
  protoDUNEStyle->SetPadTickX(1);
  protoDUNEStyle->SetPadTickY(1);

  //show the fit parameters in a box
  protoDUNEStyle->SetOptFit(1111);

  //turn off all other stats
  //protoDUNEStyle->SetOptStat(0000000);
  protoDUNEStyle->SetStatBorderSize(1);
  protoDUNEStyle->SetStatFont      (62);
  protoDUNEStyle->SetOptStat       (111111);
  protoDUNEStyle->SetStatColor     (0);
  protoDUNEStyle->SetStatX         (0.93);
  protoDUNEStyle->SetStatY         (0.90);
  protoDUNEStyle->SetLegendBorderSize(0);

  //marker settings
  protoDUNEStyle->SetMarkerStyle(20);
  protoDUNEStyle->SetMarkerSize(0.9);

  // Fonts
  const int kProtoDUNEFont = 62;

  protoDUNEStyle->SetStatFont(kProtoDUNEFont);
  protoDUNEStyle->SetLabelFont(kProtoDUNEFont,"xyz");
  protoDUNEStyle->SetTitleFont(kProtoDUNEFont,"xyz");
  protoDUNEStyle->SetTitleFont(kProtoDUNEFont,"");
  protoDUNEStyle->SetTextFont(kProtoDUNEFont);

  // ---------------------------------------------------------
  // Additions from George 26/06/2019
  protoDUNEStyle->SetCanvasBorderSize(0);
  protoDUNEStyle->SetFrameBorderSize(0);
  protoDUNEStyle->SetDrawBorder(0);
  protoDUNEStyle->SetTitleBorderSize(0);

  // Set the size (in pixels) of the small lines drawn at the end of the error bars
  protoDUNEStyle->SetEndErrorSize(4);

  // Set option to strip decimals when drawing axis labels.
  protoDUNEStyle->SetStripDecimals(kFALSE);

  protoDUNEStyle->SetLegendBorderSize(0);
  protoDUNEStyle->SetLegendFont(kProtoDUNEFont);

  protoDUNEStyle->SetLabelOffset(0.015, "x");
  protoDUNEStyle->SetLabelOffset(0.015, "y");
  protoDUNEStyle->SetLabelOffset(0.01, "z");

  protoDUNEStyle->SetTitleStyle(0);
  protoDUNEStyle->SetTitleFont(kProtoDUNEFont, "pad");
  protoDUNEStyle->SetTitleX(0.1f);
  protoDUNEStyle->SetTitleY(.98f);
  protoDUNEStyle->SetTitleW(0.8f);
  protoDUNEStyle->SetLineStyleString(2, "[12 12]"); // postscript dashes

  // get rid of X error bars and y error bar caps
  protoDUNEStyle->SetErrorX(0.001);

  protoDUNEStyle->SetNumberContours(255);
  protoDUNEStyle->SetPalette(kBird);


  // ---------------------------------------------------------
  //done
  protoDUNEStyle->cd();

  gROOT->ForceStyle();
  gStyle->ls();

/*
  if (gROOT->GetVersionInt()>51200) {
    TColor::InitializeColors();
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
  }
*/


  // Avoid too many decimal places in the axis labels
  //TGaxis::SetMaxDigits(4);
}
