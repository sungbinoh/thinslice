#ifndef BEAMTREEREPRODUCER_H
#define BEAMTREEREPRODUCER_H

#include "TH1D.h"
#include "util.h"
#include "TFile.h"
#include "BeamNtuple.h"
#include "BeamVirtualDetector.h"
#include "FillHist_Helper.h"

class BeamNtuple;
class BeamVirtualDetector;
class FillHist_Helper;
//class HadAna;
class TH1D;
class TH2D;

class BeamTreeReproducer {

 public:

  BeamTreeReproducer();

  FillHist_Helper Hist;
  //HadAna hadana;
  
  std::string fOutputFileName;
  void SetOutputFileName(std::string name){fOutputFileName = name;};
  void BookHistograms();
  void FillHistograms(const BeamNtuple & evt);
  void FillHistograms(const BeamVirtualDetector & evt, TString detector_str);
  void SaveHistograms();
  void ProcessEvent(const BeamNtuple & evt);
  void Run(BeamNtuple & evt, Long64_t nentries);
  void Run(BeamVirtualDetector & evt, Long64_t nentries, TString detector_str);

 private:
  
  TFile *outputFile;

  // ==== Histograms for normalization
  TH1D *h_cutflow;

  // ==== Histograms for beam
  
};

#endif
