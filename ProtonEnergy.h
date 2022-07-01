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
  TH1D *htrack_length_ratio[p::nIntTypes+1];

};

#endif
