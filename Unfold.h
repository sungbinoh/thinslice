#ifndef UNFOLD_H
#define UNFOLD_H

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"

class Unfold {
  
 public:

  Unfold(int nb, double xlo, double xhi);

  RooUnfoldResponse response_SliceID_Int;  //Interaction
  RooUnfoldResponse response_SliceID_Inc;  //Incident

};

#endif
