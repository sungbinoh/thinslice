#include "../../SliceParams.h"

void plotXS_data(){

  gStyle->SetOptStat(0);

  TFile *file = TFile::Open("/dune/app/users/yinrui/thinslice/build/XS.root");
  
  TGraphErrors *gr_truexs = (TGraphErrors*)file->Get("gr_truexs");
  TGraphErrors *gr_recoxs = (TGraphErrors*)file->Get("gr_recoxs");
  
  TFile f2("/dune/app/users/yinrui/thinslice/files/exclusive_xsec.root");
  TGraph *total_inel_KE = (TGraph*)f2.Get("total_inel_KE");
  
  /*for (int i: {0, 12, 12, 12, 12, 12, 12, 12}) {
    gr_truexs->RemovePoint(i);
    gr_recoxs->RemovePoint(i);
  }*/
  ofstream myfile;
  myfile.open ("toys/example.txt", ios::app);
  
  double MCg4rw[20] = {
    1,1,1,1,1,1,1,1,1,1,
    1,1,1,1,1,1,1,1,1,1,
  };
  for (int i=0; i<pi::nthinslices; ++i) {
    double KE = gr_recoxs->GetPointX(i);
    double xs = gr_recoxs->GetPointY(i);
    double xs_err = gr_recoxs->GetErrorY(i);
    double xs_true = total_inel_KE->Eval(KE);
    bool cover = abs(xs - xs_true) < xs_err;
    myfile<<cover<<", ";
  }
  myfile<<endl;
  myfile.close();
}
  
