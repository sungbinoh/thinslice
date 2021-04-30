#ifndef EVENTTYPE_H
#define EVENTTYPE_H

#include "anavar.h"

const unsigned int nParTypes = 13;

enum parType{
  kPrimMuP = 0,
  kPrimPiPInEl,
  kPrimPiPEl,
  kPrimKaPInEl,
  kPrimKaPEl,
  kPrimProInEl,
  kPrimProEl,
  kPrimMuM,
  kPrimPiMInEl,
  kPrimPiMEl,
  kPrimKaMInEl,
  kMisID,
  kUnknown
};

int GetParType(const anavar &t){

  if (!t.reco_beam_true_byE_matched){
    return kMisID;
  }
  else if (t.true_beam_PDG == -13){
    return kPrimMuP;
  }
  else if (t.true_beam_PDG == 13){
    return kPrimMuM;
  }
  else if (t.true_beam_PDG == 211){
    if ((*t.true_beam_endProcess) == "pi+Inelastic"){
      return kPrimPiPInEl;
    }
    else return kPrimPiPEl;
  }
  else if (t.true_beam_PDG == 2212){
    if ((*t.true_beam_endProcess) == "protonInelastic"){
      return kPrimProInEl;
    }
    else return kPrimProEl;
  }
  
  return kUnknown;
}

#endif
