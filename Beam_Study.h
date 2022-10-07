#ifndef BEAM_STUDY_H
#define BEAM_STUDY_H

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

class Beam_Study {

 public:

  Beam_Study();

  HadAna hadana;
  FillHist_Helper Hist;
  
  std::string fOutputFileName;
  void SetOutputFileName(std::string name){fOutputFileName = name;};

  // == Draw Histograms
  void BookHistograms();
  void FillHistBeam(const anavar & evt, double weight, TString suffix, double KE_fit_gaussian);
  void SaveHistograms();

  // == Run
  bool Pass_Beam_PID(const anavar & evt, int PID);
  void ProcessEvent(const anavar & evt);
  void Run_Beam(const anavar & evt, double weight, TString suffix, int beam_PID);
  void Run(anavar & evt, Long64_t nentries);

 private:
  
  TFile *outputFile;

  // == Beam P scale
  TString P_beam_str = "1.0";
  double P_beam_inst_scale = 1.0;
  double P_beam_inst_cut_upper = 9999.;
  double P_beam_inst_cut_lower = 0.;

  // == Beam P_{true} reweight
  double Gaussian_Reweight(double mu1, double sigma1, double mu2, double sigma2, double x_min, double x_max, double x);
  double Beam_TrueP_Reweight(const anavar & evt);
};

#endif
