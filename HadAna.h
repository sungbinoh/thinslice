#ifndef HADANA_H
#define HADANA_H

#include "EventType.h"
#include "EventSelection.h"
#include "SliceParams.h"
#include "anavar.h"

class HadAna : public anavar{
 public: 

  HadAna(TTree *tree=0);

  void AddTruePDG(int pdg);

  //Check if the desired particle is selected
  bool isSelectedPart() const;

  // Set beam cut values
  void SetBeamQualityCuts(double dx_min = -3, double dx_max = 3,
                          double dy_min = -3, double dy_max = 3,
                          double dz_min = -3, double dz_max = 3,
                          double costh_min = 0.95, double costh_max = 2);

  int GetParType();

  // Pandora slice pdg
  int pandora_slice_pdg;
  void SetPandoraSlicePDG(int pdg);

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
  double beam_dx, beam_dy, beam_dz, beam_costh;

  using anavar::anavar;

 private:
  
  //Selected true pdg list
  std::vector<int> truepdglist;

  double beamcut_dx_min, beamcut_dx_max;
  double beamcut_dy_min, beamcut_dy_max;
  double beamcut_dz_min, beamcut_dz_max;
  double beamcut_costh_min, beamcut_costh_max;

};

#endif
