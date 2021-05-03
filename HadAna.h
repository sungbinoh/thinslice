#ifndef HADANA_H
#define HADANA_H

#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "EventType.h"
#include "EventSelection.h"
#include "anavar.h"

class HadAna : public anavar{
 public: 
  //Selected true pdg list
  std::vector<int> truepdglist;
  void AddTruePDG(int pdg){
    truepdglist.push_back(pdg);
  };

  //Check is the current particle is selected
  bool isTrueSelectedPart();

  int GetParType();

  //Pandora slice pdg
  int pandora_slice_pdg;
  void SetPandoraSlicePDG(int pdg){
    pandora_slice_pdg = pdg;
  };

  bool PassPandoraSliceCut();
  bool PassBeamQualityCut();
  bool PassAPA3Cut();
  bool PassCaloSizeCut();

  TH1D *htrue_beam_endZ[nCuts][nParTypes+1];
  TH1D *hreco_beam_endZ[nCuts][nParTypes+1];
  TH1D *hreco_true_beam_endZ[nCuts][nParTypes+1];
  TH2D *hreco_vs_true_beam_endZ[nCuts][nParTypes+1];

  std::string fOutputFileName;
  TFile *outputFile;
  void SetOutputFileName(std::string name){fOutputFileName = name;};
  void BookHistograms();
  void FillHistograms(int cut);
  void SaveHistograms();

  using anavar::anavar;

};

#endif
