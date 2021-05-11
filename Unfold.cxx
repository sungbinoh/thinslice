#include "Unfold.h"

Unfold::Unfold(int nb, double xlo, double xhi){

  response_SliceID_Int.Reset();
  response_SliceID_Int.Setup(nb, xlo, xhi);

  response_SliceID_Inc.Reset();
  response_SliceID_Inc.Setup(nb, xlo, xhi);

}  
