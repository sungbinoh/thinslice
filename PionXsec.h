#ifndef PIONXSEC_H
#define PIONXSEC_H

#include "TFile.h"
#include "HadAna.h"
#include "FillHist_Helper.h"
#include "RecoDaughter.h"

class anavar;
class HadAna;
class FillHist_Helper;
//class RecoDaughter;
class TH1D;
class TH2D;

class PionXsec {

 public:

  PionXsec();

  HadAna hadana;
  FillHist_Helper Hist;
  
  std::string fOutputFileName;
  void SetOutputFileName(std::string name){fOutputFileName = name;};

  // == Call Daughter Collections
  vector<RecoDaughter> GetAllRecoDaughters(const anavar & evt);
  vector<RecoDaughter> GetPions(const vector<RecoDaughter> in);
  vector<RecoDaughter> GetProtons(const vector<RecoDaughter> in);
  vector<RecoDaughter> GetTruePions(const vector<RecoDaughter> in);
  vector<RecoDaughter> GetTrueProtons(const vector<RecoDaughter> in);

  // == Draw Histograms
  void BookHistograms();
  void FillHistBeam(const anavar & evt, double weight, TString suffix);
  void FillHistDaughterTrue(const vector<RecoDaughter> daughters, const anavar & evt, double weight, TString suffix);
  void FillHistDaughterPurity(const vector<RecoDaughter> daughters, const anavar & evt, double weight, TString suffix);
  void FillHistQE_MCstudy(const vector<RecoDaughter> daughters, const anavar & evt, double weight, TString suffix);
  void SaveHistograms();

  // == Run
  void ProcessEvent(const anavar & evt);
  void Run(anavar & evt, Long64_t nentries);

 private:
  
  TFile *outputFile;

  // == Beam P scale
  double beamP_scale = 1.0;

  // == Beam P_{true} reweight
  double Gaussian_Reweight(double mu1, double sigma1, double mu2, double sigma2, double x_min, double x_max, double x);

  // == Event Topology
  double Get_EQE(double P_pion, double cos_theta);
  double Get_EQE_NC_Pion(double P_pion, double cos_theta, double E_binding, int which_sol);

};

#endif
