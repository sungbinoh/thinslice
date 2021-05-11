#include "Unfold.h"

Unfold::Unfold(int nb, double xlo, double xhi){

  response_SliceID_Int.Reset();
  response_SliceID_Int.Setup(nb, xlo, xhi);

  response_SliceID_Inc.Reset();
  response_SliceID_Inc.Setup(nb, xlo, xhi);

  eff_num_Int = new TH1D("eff_num_Int", "eff_num_Int", nb, xlo, xhi);
  eff_den_Int = new TH1D("eff_den_Int", "eff_den_Int", nb, xlo, xhi);
  eff_num_Inc = new TH1D("eff_num_Inc", "eff_num_Inc", nb, xlo, xhi);
  eff_den_Inc = new TH1D("eff_den_Inc", "eff_den_Inc", nb, xlo, xhi);
  pur_num_Int = new TH1D("pur_num_Int", "pur_num_Int", nb, xlo, xhi);
  pur_num_Inc = new TH1D("pur_num_Inc", "pur_num_Inc", nb, xlo, xhi);
  pur_den     = new TH1D("pur_den",     "pur_den",     nb, xlo, xhi);

}  

void Unfold::SaveHistograms(){

  eff_num_Int->Write("eff_num_Int");
  eff_den_Int->Write("eff_den_Int");
  eff_num_Inc->Write("eff_num_Inc");
  eff_den_Inc->Write("eff_den_Inc");
  pur_num_Int->Write("pur_num_Int");
  pur_num_Inc->Write("pur_num_Inc");
  pur_den->Write("pur_den");

}
