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
  //void Run(vector<BeamVirtualDetector> & branch_vector, const vector<TString> branch_name_vector, Long64_t nentries);
  void Find_Matched_Event(BeamVirtualDetector & detector_tree, int Evt_ID, int Trk_ID, int PDG_ID, TString detector_str);
  void Run(BeamVirtualDetector & evt_AfterTarget, BeamVirtualDetector & evt_TOF1, BeamVirtualDetector & evt_COLL1, BeamVirtualDetector & evt_BPROF1,
	   BeamVirtualDetector & evt_BPROF2, BeamVirtualDetector & evt_BPROF3, BeamVirtualDetector & evt_TRIG1, BeamVirtualDetector & evt_BPROFEXT,
	   BeamVirtualDetector & evt_BPROF4, BeamVirtualDetector & evt_TRIG2, Long64_t nentries);
  

 private:
  
  TFile *outputFile;

  // ==== Histograms for normalization
  TH1D *h_cutflow;

  // ==== Histograms for beam
  
};

#endif
