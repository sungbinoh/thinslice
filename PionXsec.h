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
  void Fill_Eslice_Histograms(TString histname, TString suffix,double KE_init, double KE_end, double KE_int, double weight, const vector<double> binning, int N_bin, bool fill_int);
  void Fill_Eslice_Study(const anavar & evt, double weight, TString suffix);
  void FillHistQE_MCstudy(const anavar & evt, double weight, TString suffix);
  void FillHistQE_Reco(const vector<RecoDaughter> daughters, const anavar & evt, double weight, TString suffix);
  void SaveHistograms();

  // == Run
  bool Pass_Beam_PID(const anavar & evt, int PID);
  void ProcessEvent(const anavar & evt);
  void Run_Beam(const anavar & evt, double weight, TString suffix);
  void Run_Daughter(const vector<RecoDaughter> daughters, const anavar & evt, double weight, TString suffix);
  void Run(anavar & evt, Long64_t nentries);

  //int N_binning_100MeV = 10;
  //vector<double> binning_100MeV = {0., 100., 200., 300., 400., 500., 600., 700., 800., 900., 1000.};

 private:
  
  TFile *outputFile;

  // == Beam P scale
  TString P_beam_str = "1.0";
  double P_beam_inst_scale = 1.0;
  double P_beam_inst_cut_upper = 9999.;
  double P_beam_inst_cut_lower = 0.;

  // == Beam P to P_ff
  double Convert_P_Spectrometer_to_P_ff(const anavar & evt, double P_beam_inst, TString suffix, TString key, int syst);

  // == Beam P_{true} reweight
  double Gaussian_Reweight(double mu1, double sigma1, double mu2, double sigma2, double x_min, double x_max, double x);

  // == Event Topology
  double Get_EQE(double P_pion, double cos_theta);
  double Get_EQE_NC_Pion(double P_pion, double cos_theta, double E_binding, int which_sol);

  // == Binnings for cross sections
  int N_binning_100MeV = 10;
  vector<double> binning_100MeV = {0., 100., 200., 300., 400., 500., 600., 700., 800., 900., 1000.};
  

};

#endif
