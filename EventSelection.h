#ifndef EVENTSELECTION_H
#define EVENTSELECTION_H

#include "anavar.h"

bool isTrueSelected(const anavar &t, const std::vector<int> & pdglist){

  for (size_t i = 0; i<pdglist.size(); ++i){
    if (t.true_beam_PDG == pdglist[i]) return true;
  }
  return false;
}

#endif
