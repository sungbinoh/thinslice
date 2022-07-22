#ifndef BEAMSAMPLEANA_H
#define BEAMSAMPLEANA_H

#include "TH1D.h"
#include "util.h"
#include "TFile.h"
#include "BeamNtuple.h"
#include "FillHist_Helper.h"

class BeamNtuple;
class FillHist_Helper;
//class HadAna;
class TH1D;
class TH2D;

class BeamSampleAna {

 public:

  BeamSampleAna();

  FillHist_Helper Hist;
  //HadAna hadana;
  
  std::string fOutputFileName;
  void SetOutputFileName(std::string name){fOutputFileName = name;};
  void BookHistograms();
  void FillHistograms(const BeamNtuple & evt);
  void SaveHistograms();
  void ProcessEvent(const BeamNtuple & evt);
  void Run(BeamNtuple & evt, Long64_t nentries);

 private:
  
  TFile *outputFile;

  // ==== Histograms for normalization
  TH1D *h_cutflow;

  // ==== Histograms for beam
  
};

#endif
