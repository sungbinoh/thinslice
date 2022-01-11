#ifndef EFFEVAL_H
#define EFFEVAL_H

#include "TFile.h"
#include "SliceParams.h"
#include "HadAna.h"

class anavar;
class HadAna;
class TEfficiency;
class TH1D;
class TH2D;

class EffEval {

 public:

  EffEval();

  HadAna hadana;
  
  std::string fOutputFileName;
  void SetOutputFileName(std::string name){fOutputFileName = name;};
  void BookHistograms();
  void FillHistograms(const anavar & evt);
  void SaveHistograms();

  void ProcessEvent(const anavar & evt);

  void Run(anavar & evt, Long64_t nentries);

 private:
  
  TFile *outputFile;
  TEfficiency *h_Eff_Ppi;
  TEfficiency *h_Eff_Ppi_michel;
  TH1D *h_true_Ppi_all;
  TH1D *h_true_Ppi_sel;
  TH2D *h_reco_true_Ppi_sel;
  TH1D *h_res_Ppi_sel;
  TH2D *h_res_Ppi_michelscore;
  TH1D *h_true_Ppi_michel;
  TH2D *h_reco_true_Ppi_michel;
  TH1D *h_res_Ppi_michel;

};

#endif
