#ifndef PROTONENERGY_H
#define PROTONENERGY_H

#include "TFile.h"
#include "HadAna.h"

class anavar;
class HadAna;
class TH1D;
class TH2D;

class ProtonEnergy {

 public:

  ProtonEnergy();

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
  // ratio of the reconstructed track length and the predicted track length based on the beam momentum
  TH1D *htrack_length_ratio[p::nIntTypes+1];
  // same as above but consider a 46 MeV upstream energy loss
  TH1D *htrack_length_ratio_eloss46[p::nIntTypes+1];
  // Proton energy at the track end point based on the beam momentum and track length
  TH1D *hend_energy[p::nIntTypes+1];
};

#endif
