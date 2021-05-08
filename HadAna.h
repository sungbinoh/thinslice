#ifndef HADANA_H
#define HADANA_H

#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "EventType.h"
#include "EventSelection.h"
#include "SliceParams.h"
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

  // Pandora slice pdg
  int pandora_slice_pdg;
  void SetPandoraSlicePDG(int pdg){
    pandora_slice_pdg = pdg;
  };

  bool PassPandoraSliceCut() const;
  bool PassBeamQualityCut() const;
  bool PassAPA3Cut() const;
  bool PassCaloSizeCut() const;
  bool PassMichelScoreCut() const;
  bool PassMediandEdxCut() const;
  bool PassAllCuts() const;

  // Event information
  void ProcessEvent();
  int partype;
  double median_dEdx;
  double daughter_michel_score;

  using anavar::anavar;

};

#endif
