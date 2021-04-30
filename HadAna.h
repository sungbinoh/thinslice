#ifndef HADANA_H
#define HADANA_H

#include "anavar.h"

class HadAna : public anavar{
 public: 
  //Selected true pdg list
  std::vector<int> truepdglist;
  void AddTruePDG(int pdg);

  //Check is the current particle is selected
  bool isTrueSelectedPart();

  using anavar::anavar;

};

#endif
