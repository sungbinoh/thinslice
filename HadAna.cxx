#include "HadAna.h"

void HadAna::AddTruePDG(int pdg){
  truepdglist.push_back(pdg);
}

bool HadAna::isTrueSelectedPart(){
  for (size_t i = 0; i<truepdglist.size(); ++i){
    if (true_beam_PDG == truepdglist[i]) return true;
  }
  return false;
}
